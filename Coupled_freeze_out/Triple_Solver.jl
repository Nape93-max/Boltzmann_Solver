using Plots
using LaTeXStrings

include("../FreezeOut.jl") #useful functions of the 1D solver

# Useful constants
const SUNN = 3 #No of colours of the group
const SUNN2 = 2*SUNN
const CF = (SUNN*SUNN-1)/SUNN2
const kappa = (SUNN*SUNN-1)*(SUNN*SUNN-2)/(16*SUNN^3)

const g_quark = 4*SUNN #degeneracy of the Dirac quark
const g_Diquark = SUNN2*(SUNN2-1) #Degeneracy of the Diquark (check all possible combinations of colour, spin and particle-antiparticle)
const g_Baryon = 4*SUNN/factorial(SUNN-2) #Degeneracy of the baryon

function YeqQuark(g_deg, gS, x) #Calculates the equilibrium yield of species with degeneracy g_deg at time x with gS relativistic entropic d.o.f.
    return 90/(2*pi)^3.5*g_deg/gS*x*sqrt(x)*exp(-x)
end

function YeqDiQuark(g_deg, gS, x, Ebin) #Calculates the equilibrium yield of species with degeneracy g_deg at time x with gS relativistic entropic d.o.f.
    y = (2-Ebin)*x
    return 90/(2*pi)^3.5*g_deg/gS*y*sqrt(y)*exp(-y)
end

function YeqBaryon(g_deg, gS, x, Ebin) #Calculates the equilibrium yield of species with degeneracy g_deg at time x with gS relativistic entropic d.o.f.
    y = (3-Ebin)*x
    return 90/(2*pi)^3.5*g_deg/gS*y*sqrt(y)*exp(-y)
end

function Multi_dim_Newton_Raphson_step(x, W_old, cs, g_deg, gS, Ebin, procs) #Method to calculate the next 3 W=log(Y) point in the implicit backward Euler method.
    # x is the x_(n+1) time step, W_old is W_n (a 3 vector!) and cs is Delta_t*lambda(t_(n+1))*geff(t_(n+1)).
    # cs is an array of length corresponding to the number of different cross sections (currently 3). g_deg is the degeneracy of the annihilating particles (also a 3-vector!).
    # gS are the effective entropy dofs. Ebin is a 3-vector of binding energies of the 3 species (Ebin[1] = 0)
    # procs is a boolean array of length 14 stating which processes are active and which are not  
    W_try_initial = W_old
    W_new = Multi_dim_Newton_Raphson_iteration(x, W_old, cs, W_try_initial, g_deg, gS, Ebin, procs)
    diff = sum(abs.(log.(abs.(W_new./W_try_initial))))

    while diff > 1E-2 
        W_try = deepcopy(W_new);
        W_new = Multi_dim_Newton_Raphson_iteration(x, W_old, cs, W_try, g_deg, gS, Ebin, procs)
        diff = sum(abs.(log.(abs.(W_new./W_try))))
    end

    if any(isnan.(W_new))
        println("ALARM: W_new in function Newton_Raphson is NaN. Maybe relax bound on diff in the while loop as a fix.")
    end
    
    return W_new
end

function Multi_dim_Newton_Raphson_iteration(x, W_old, cs, W_previous, g_deg, gS, Ebin, procs) #Does one NR-step to calculate a trial W^(i+1)_(n+1). W_previous = W^(i)_(n+1)
    # All W are 3-vectors!
    J = JMatrix(x, cs, W_previous, g_deg, gS, Ebin, procs) #Set the Jacobian matrix

    J_inverse = Matrix_Inversion(J) #Set the inverse of the Jacobian matrix
    
    W_next = W_previous .- (J_inverse*Multi_dim_NR_func_array(W_previous, W_old, cs, g_deg, gS, x, Ebin, procs))

    return W_next
end

function Multi_dim_NR_func_array(W_prev, W_old, lam, g_deg, gS, x, Ebin, procs) 
    #auxiliary function that returns the 3-vector of equations that have to be solved in the multi-dim Euler backward step

    final_func = zeros(3) #Initialise final output vector

    final_func[1] = W_prev[1] - W_old[1] - rhsQuark(W_prev, lam, g_deg, gS, x, Ebin, procs)
    final_func[2] = W_prev[2] - W_old[2] - rhsDiQuark(W_prev, lam, g_deg, gS, x, Ebin, procs)
    final_func[3] = W_prev[3] - W_old[3] - rhsBaryon(W_prev, lam, g_deg, gS, x, Ebin, procs)

    return final_func
end

function rearrangement_cross_section(x, m, alpha)
    #return kappa*pert_acs(alpha,m)
    return pi/(m^2*alpha*CF*sqrt(SUNN))
end

function Coupled_Freeze_Out(x_init, x_final, Delta_t, mq, Lambda, procs) #Function that solves the coupled system of Boltzmann equations between x_init and x_final
    # mq is the quark mass in GeV, Lambda is the confinement scale in GeV. Delta_t is the logarithmic step size between x_init and x_final
    # procs is an array of length 14 indicating which processes are active.  
    Alpha_ann = running_coupling_from_pole(2*mq, Lambda, (11*SUNN-2)/3) #coupling at quark-antiquark annihilation scale

    #Binding energies are defined as E_B/m_quark, since this cancels in the definition of Yeq(x). Values taken from "Weakly coupled dark baryon"
    E_binding_meson = 0.25*(Alpha_ann*CF)^2
    E_binding_diquark = 1/16*(Alpha_ann*CF)^2
    E_binding_baryon = 0.26*(Alpha_ann*CF)^2

    sigma_ann = kappa*pert_acs(Alpha_ann, mq) 
    sigma_cap = deepcopy(sigma_ann)

    t_initial = log(x_init)
    t_final = log(x_final)

    #First define initial conditions:
    tvec = collect(t_initial:Delta_t:t_final)
    xvec = exp.(tvec)
    Npoints = length(tvec)

    g_star_eff_vec = 106.75*ones(Npoints) 
    entropy_eff_dof = 106.75*ones(Npoints) # h_eff_dof.(m_quark./xvec) #set effective entropic degrees of freedom
    #g_star_eff_vec = eff_dof_sqrt.(m_quark./xvec) #Effective degrees of freedom
    #entropy_eff_dof = h_eff_dof.(m_quark./xvec) #set effective entropic degrees of freedom

    #Initialise the different cross sections
    sigma_v_ann = ones(Npoints) #Initialise array for interpolated annihaltion xs
    #sigma_v_averaged_coeffs = sigma_v_interpolation(sigma0, sigma_v_x_values, sigma_v_y_values) #Cubic spline fit coefficient vector (beta, gamma, delta)
    sigma_v_cap = ones(Npoints) #Initialise array for interpolated capture xs
    sigma_v_RA = ones(Npoints) #Initialise array for interpolated rearrangement xs

    #Set equilibrium yields and cross sections
    EquilibriumYieldQuark = zeros(Npoints)
    EquilibriumYieldDiQuark = zeros(Npoints)
    EquilibriumYieldBaryon = zeros(Npoints)

    for i = 1:Npoints
        EquilibriumYieldQuark[i] = YeqQuark(g_quark, entropy_eff_dof[i], xvec[i])
        EquilibriumYieldDiQuark[i] = YeqDiQuark(g_Diquark, entropy_eff_dof[i], xvec[i], E_binding_diquark)
        EquilibriumYieldBaryon[i] = YeqBaryon(g_Baryon, entropy_eff_dof[i], xvec[i], E_binding_baryon)
        
        #Assign the cross sections
        #sigma_v_averaged[i] = sigma_v_cspline(xvec[i], sigma_v_x_values, sigma_v_y_values, sigma0, sigma_v_averaged_coeffs)
        #sigma_v_averaged[i] = sigma_ann
        sigma_v_ann[i] = sigma_ann
        sigma_v_cap[i] = sigma_cap
        sigma_v_RA[i] = rearrangement_cross_section(xvec[i], mq, Alpha_ann)
    end

    #Set the initial conditions
    WxQuark = zeros(Npoints)
    YxQuark = zeros(Npoints)
    YxQuark[1] = EquilibriumYieldQuark[1]
    WxQuark[1] = log(YxQuark[1])

    WxDiQuark = zeros(Npoints)
    YxDiQuark = zeros(Npoints)  
    YxDiQuark[1] = EquilibriumYieldDiQuark[1]
    WxDiQuark[1] = log(YxDiQuark[1])

    WxBaryon = zeros(Npoints)
    YxBaryon = zeros(Npoints)
    YxBaryon[1] = EquilibriumYieldBaryon[1]
    WxBaryon[1] = log(YxBaryon[1])

    #Define sets of quantities needed for the solution
    Wx = [WxQuark, WxDiQuark, WxBaryon] #Complete solution vector
    sigma_v_averaged = [sigma_v_ann, sigma_v_cap, sigma_v_RA] #Complete vector of all the distinct cross sections
    degeneracy_vec = [g_quark, g_Diquark, g_Baryon] #Complete vector of all 3 distinct dofs
    E_bin_vec = zeros(3) #Complete vector of binding energies
    E_bin_vec[2] = E_binding_diquark
    E_bin_vec[3] = E_binding_baryon

    BigConstant = 0.5*Delta_t*bc_constant(mq) #This constant is a UNIVERSAL part of the prefactor of all the Boltzmann equations, because x = m_quark/T_SM. 

    #Solution to the Boltzmann equations for the decoupled freeze-out
    for i = 2:Npoints
        Wx_new = zeros(3)
        Wx_old = zeros(3)
        lambda = zeros(3) # equals BigConstant/xvec[i]*sqrt(g_star)*sigma_v

        h1 = BigConstant*g_star_eff_vec[i]/xvec[i]
        for j = 1:3
            Wx_old[j] = Wx[j][i-1]
            lambda[j] = h1*sigma_v_averaged[j][i]
        end
    
        Wx_new = Multi_dim_Newton_Raphson_step(xvec[i], Wx_old, lambda, degeneracy_vec, entropy_eff_dof[i], E_bin_vec, procs)

        for j = 1:3
            Wx[j][i] = Wx_new[j]
        end
    end

    return exp.(Wx[1]), exp.(Wx[2]), exp.(Wx[3])
end

function f4(x, gS, Ebin) #This function returns the equilibrium ratio Y_quark^2/Y_Diquark. Appears in process No. 4
    pre_fac = 90/(2*pi)^3.5*g_quark^2/(g_Diquark*gS)
    x_fac = (x/(2 - Ebin))^1.5
    exp_fac = exp(-Ebin*x)
    return pre_fac*x_fac*exp_fac
end

function f5(x, gS, Ebin_D, Ebin_B) #This function returns the rquilibrium ratio Y_quark*X_Diquark/Y_Baryon. Appears in process No. 4
    pre_fac = 90/(2*pi)^3.5*g_quark*g_Diquark/(g_Baryon*gS)
    x_fac = (x*(2-Ebin_D)/(3 - Ebin_B))^1.5
    exp_fac = exp((Ebin_D - Ebin_B)*x)
    return pre_fac*x_fac*exp_fac
end

function rhsQuark(W_prev, lam, g_deg, gS, x, Ebin, procs) #!!! UNDER CONSTRUCTION !!! 
    #First define relevant quantities
    sigma_ann = lam[1]
    sigma_cap = lam[2]
    sigma_RA = lam[3]
    YQuark = exp(W_prev[1])
    YQuarkEquilibrium = YeqQuark(g_deg[1], gS, x)
    YDiQuark = exp(W_prev[2])
    YDiQuarkEquilibrium = YeqDiQuark(g_deg[2], gS, x, Ebin[2])
    YBaryon = exp(W_prev[3])
    YBaryonEquilibrium = YeqBaryon(g_deg[3], gS, x, Ebin[3])

    rhs = 0 #Initialise final output

    if procs[1] #Process No. 1 (quark-antiquark annihilation)
        rhs -= generic_process(sigma_ann, YQuark, YQuark, 1., 1., YQuarkEquilibrium^2)
    end

    if procs[4] #Process No. 4 (quark-capture into diquark)
        eq_ratio4 = f4(x, gS, Ebin[2])
        rhs -= 2*generic_process(sigma_cap, YQuark, YQuark, YDiQuark, 1., eq_ratio4)
    end

    if procs[5] #Process No. 5 (quark-diquark capture into baryons)
        eq_ratio5 = f5(x, gS, Ebin[2], Ebin[3])
        rhs -= generic_process(sigma_cap, YQuark, YDiQuark, YBaryon, 1., eq_ratio5)
    end

    return rhs
end

function rhsDiQuark(W_prev, lam, g_deg, gS, x, Ebin, procs) #!!! UNDER CONSTRUCTION !!! 
    #First define relevant quantities
    sigma_ann = lam[1]
    sigma_cap = lam[2]
    sigma_RA = lam[3]
    YQuark = exp(W_prev[1])
    YQuarkEquilibrium = YeqQuark(g_deg[1], gS, x)
    YDiQuark = exp(W_prev[2])
    YDiQuarkEquilibrium = YeqDiQuark(g_deg[2], gS, x, Ebin[2])
    YBaryon = exp(W_prev[3])
    YBaryonEquilibrium = YeqBaryon(g_deg[3], gS, x, Ebin[3])

    rhs = 0 #Initialise final output

    if procs[2]
        rhs -= generic_process(sigma_RA, YDiQuark, YDiQuark, 1., 1., YDiQuarkEquilibrium^2)
    end

    if procs[4] #Process No. 4 (quark-capture into diquark)
        eq_ratio4 = f4(x, gS, Ebin[2])
        rhs += generic_process(sigma_cap, YQuark, YQuark, YDiQuark, 1, eq_ratio4)
    end

    if procs[5] #Process No. 5 (quark-diquark capture into baryons)
        eq_ratio5 = f5(x, gS, Ebin[2], Ebin[3])
        rhs -= generic_process(sigma_cap, YQuark, YDiQuark, YBaryon, 1., eq_ratio5)
    end

    return rhs
end

function rhsBaryon(W_prev, lam, g_deg, gS, x, Ebin, procs) #!!! UNDER CONSTRUCTION !!! 
    #First define relevant quantities
    sigma_ann = lam[1]
    sigma_cap = lam[2]
    sigma_RA = lam[3]
    YQuark = exp(W_prev[1])
    YQuarkEquilibrium = YeqQuark(g_deg[1], gS, x)
    YDiQuark = exp(W_prev[2])
    YDiQuarkEquilibrium = YeqDiQuark(g_deg[2], gS, x, Ebin[2])
    YBaryon = exp(W_prev[3])
    YBaryonEquilibrium = YeqBaryon(g_deg[3], gS, x, Ebin[3])

    rhs = 0 #Initialise final output

    if procs[2]
        rhs -= generic_process(sigma_RA, YBaryon, YBaryon, 1., 1., YBaryonEquilibrium^2)
    end

    if procs[5] #Process No. 5 (quark-diquark capture into baryons)
        eq_ratio5 = f5(x, gS, Ebin[2], Ebin[3])
        rhs += generic_process(sigma_cap, YQuark, YDiQuark, YBaryon, 1., eq_ratio5)
    end

    return rhs 
end

function JMatrix(x, cs, W_previous, g_deg, gS, Ebin, procs) #This function calculates the Jacobian matrix in each Euler backward timestep
    #Argument definitions as in function "Multi_dim_Newton_Raphson_iteration"
    J = zeros(3,3)
    J[1,1] = J11(x, cs, W_previous, g_deg, gS, Ebin, procs)
    J[1,2] = J12(x, cs, W_previous, g_deg, gS, Ebin, procs)
    J[1,3] = J13(x, cs, W_previous, g_deg, gS, Ebin, procs)
    J[2,1] = J21(x, cs, W_previous, g_deg, gS, Ebin, procs)
    J[2,2] = J22(x, cs, W_previous, g_deg, gS, Ebin, procs)
    J[2,3] = J23(x, cs, W_previous, g_deg, gS, Ebin, procs)
    J[3,1] = J31(x, cs, W_previous, g_deg, gS, Ebin, procs)
    J[3,2] = J32(x, cs, W_previous, g_deg, gS, Ebin, procs)
    J[3,3] = J33(x, cs, W_previous, g_deg, gS, Ebin, procs)
    return J
end

function Matrix_Inversion(J) #Returns the inverse of the 3x3 J input matrix.  
    #For a 3x3 matrix, we use A^(-1) = 1/det(A)*adj(A)
    det = Determinant_Sarrus(J)

    temp_mat = Adjoint_Matrix(J)

    Jinv = det*temp_mat
    
    return Jinv
end

function Determinant_Sarrus(J) #Calculates the inverse determinant of any 3x3 matrix
    det = J[1,1]*J[2,2]*J[3,3] + J[2,1]*J[3,2]*J[1,3] + J[3,1]*J[1,2]*J[2,3] - J[2,3]*J[3,2]*J[1,1] - J[3,3]*J[1,2]*J[2,1] - J[3,1]*J[2,2]*J[1,3]
    det_inv = 1/det
    if isnan(1/det)
        println("HORROR: Jacobian non invertible due to vanishing determinant")
    end
    return det_inv
end

function Adjoint_Matrix(J) #Calculates the adjoint matrix of a 3x3 matrix J 
    M_adj = zeros(3,3)

    for i = 1:3
        for j = 1:3
            M_adj[j,i] = Cofactor(J, i, j) #Transposition via M(j,i)
        end
    end

    return M_adj 
end

function Cofactor(M, i, j) #Calculates the minor of the (i,j) component of a 3x3 matrix M
    #the sign of the cofactor is automatically included in the way the minor is calculated!
    cofac = M[1 + mod(i,3), 1 + mod(j,3)]*M[1 + mod(i+1,3), 1 + mod(j+1,3)] - M[1 + mod(i,3), 1 + mod(j+1,3)]*M[1 + mod(i+1,3), 1 + mod(j,3)]
    return cofac
end

function generic_process(lambda, Y1, Y2, Y3, Y4, equi_ratio) # Generic cross section of process No. 1 defined according to eq. (3.9) of the long squeezeout paper.
    #Definition in the double logartithmic scheme, i.e. Y_a == exp(W_a). The entire prefactor outside the brackets is encoded in lambda.
    #equi_ratio is any combination of equilibrium yields appearing in the expression
    #The correct sign needs to be implemented by the user!
    return lambda*(Y2 - Y3*Y4*equi_ratio/Y1)
end

function generic_derivative_of_process(lambda, Y1, Y2, Y3, Y4, equi_ratio, bool_array) #Derivative of a generic process needed in the Jacobian
    #The correct sign needs to be implemented by the user!
    #Definition in the double logartithmic scheme, i.e. Y == exp(W). The entire prefactor outside the brackets is encoded in lambda
    #bool_array is an array of length 4 encoding which of the four yields the expression is derived.

    der = 0 #initialise final result 

    h = Y3*Y4*equi_ratio/Y1 #auxiliary quantity appearing frequently
    if bool_array[1] 
        der += h
    end

    if bool_array[3]
        der -= h
    end

    if bool_array[4]
        der -= h
    end

    if bool_array[2]
        der += Y2
    end

    return lambda*der
end

function D_rhs_Quark_Quark(x, cs, W_previous, g_deg, gS, Ebin, procs) # !!! UNDER CONSTRUCTION !!! Returns the derivative of the rhs of the quark equation w.r.t. the quark field
    #First define relevant quantities
    sigma_ann = cs[1]
    sigma_cap = cs[2]
    sigma_RA = cs[3]
    YQuark = exp(W_previous[1])
    YQuarkEquilibrium = YeqQuark(g_deg[1], gS, x)
    YDiQuark = exp(W_previous[2])
    YDiQuarkEquilibrium = YeqDiQuark(g_deg[2], gS, x, Ebin[2])
    YBaryon = exp(W_previous[3])
    YBaryonEquilibrium = YeqBaryon(g_deg[3], gS, x, Ebin[3])

    der = 0 #Initialise final output

    if procs[1] #Process No. 1 (quark-antiquark annihilation)
        der -= generic_derivative_of_process(sigma_ann, YQuark, YQuark, 1, 1, YQuarkEquilibrium^2, [true, true, false, false])
    end

    if procs[4] #Process No. 4 (quark-capture into diquark)
        eq_ratio4 = f4(x, gS, Ebin[2])
        der -= generic_derivative_of_process(2*sigma_cap, YQuark, YQuark, YDiQuark, 1, eq_ratio4, [true, true, false, false])
    end

    if procs[5] #Process No. 5 (quark-diquark capture into baryons)
        eq_ratio5 = f5(x, gS, Ebin[2], Ebin[3])
        der -= generic_derivative_of_process(sigma_cap, YQuark, YDiQuark, YBaryon, 1., eq_ratio5, [true, false, false, false])
    end

    return der
end

function D_rhs_Quark_DiQuark(x, cs, W_previous, g_deg, gS, Ebin, procs) # !!! UNDER CONSTRUCTION !!! Returns the derivative of the rhs of the quark equation w.r.t. the quark field
    #First define relevant quantities
    sigma_ann = cs[1]
    sigma_cap = cs[2]
    sigma_RA = cs[3]
    YQuark = exp(W_previous[1])
    YQuarkEquilibrium = YeqQuark(g_deg[1], gS, x)
    YDiQuark = exp(W_previous[2])
    YDiQuarkEquilibrium = YeqDiQuark(g_deg[2], gS, x, Ebin[2])
    YBaryon = exp(W_previous[3])
    YBaryonEquilibrium = YeqBaryon(g_deg[3], gS, x, Ebin[3])

    der = 0 #Initialise final output

    if procs[4] #Process No. 4 (quark-capture into diquark)
        eq_ratio4 = f4(x, gS, Ebin[2])
        der -= generic_derivative_of_process(2*sigma_cap, YQuark, YQuark, YDiQuark, 1, eq_ratio4, [false, false, true, false])
    end

    if procs[5] #Process No. 5 (quark-diquark capture into baryons)
        eq_ratio5 = f5(x, gS, Ebin[2], Ebin[3])
        der -= generic_derivative_of_process(sigma_cap, YQuark, YDiQuark, YBaryon, 1., eq_ratio5, [false, true, false, false])
    end

    return der
end

function D_rhs_Quark_Baryon(x, cs, W_previous, g_deg, gS, Ebin, procs) # !!! UNDER CONSTRUCTION !!! Returns the derivative of the rhs of the quark equation w.r.t. the quark field
    #First define relevant quantities
    sigma_ann = cs[1]
    sigma_cap = cs[2]
    sigma_RA = cs[3]
    YQuark = exp(W_previous[1])
    YQuarkEquilibrium = YeqQuark(g_deg[1], gS, x)
    YDiQuark = exp(W_previous[2])
    YDiQuarkEquilibrium = YeqDiQuark(g_deg[2], gS, x, Ebin[2])
    YBaryon = exp(W_previous[3])
    YBaryonEquilibrium = YeqBaryon(g_deg[3], gS, x, Ebin[3])

    der = 0 #Initialise final output

    if procs[5] #Process No. 5 (quark-diquark capture into baryons)
        eq_ratio5 = f5(x, gS, Ebin[2], Ebin[3])
        der -= generic_derivative_of_process(sigma_cap, YQuark, YDiQuark, YBaryon, 1., eq_ratio5, [false, false, true, false])
    end

    return der
end

function D_rhs_DiQuark_Quark(x, cs, W_previous, g_deg, gS, Ebin, procs) # !!! UNDER CONSTRUCTION !!! Returns the derivative of the rhs of the quark equation w.r.t. the quark field
    #First define relevant quantities
    sigma_ann = cs[1]
    sigma_cap = cs[2]
    sigma_RA = cs[3]
    YQuark = exp(W_previous[1])
    YQuarkEquilibrium = YeqQuark(g_deg[1], gS, x)
    YDiQuark = exp(W_previous[2])
    YDiQuarkEquilibrium = YeqDiQuark(g_deg[2], gS, x, Ebin[2])
    YBaryon = exp(W_previous[3])
    YBaryonEquilibrium = YeqBaryon(g_deg[3], gS, x, Ebin[3])

    der = 0 #Initialise final output

    if procs[4] #Process No. 4 (quark-capture into diquark)
        eq_ratio4 = f4(x, gS, Ebin[2])
        der += generic_derivative_of_process(sigma_cap, YQuark, YQuark, YDiQuark, 1, eq_ratio4, [true, true, false, false])
    end

    if procs[5] #Process No. 5 (quark-diquark capture into baryons)
        eq_ratio5 = f5(x, gS, Ebin[2], Ebin[3])
        der -= generic_derivative_of_process(sigma_cap, YQuark, YDiQuark, YBaryon, 1., eq_ratio5, [true, false, false, false])
    end
    
    return der
end

function D_rhs_DiQuark_DiQuark(x, cs, W_previous, g_deg, gS, Ebin, procs) # !!! UNDER CONSTRUCTION !!! Returns the derivative of the rhs of the quark equation w.r.t. the quark field
    #First define relevant quantities
    sigma_ann = cs[1]
    sigma_cap = cs[2]
    sigma_RA = cs[3]
    YQuark = exp(W_previous[1])
    YQuarkEquilibrium = YeqQuark(g_deg[1], gS, x)
    YDiQuark = exp(W_previous[2])
    YDiQuarkEquilibrium = YeqDiQuark(g_deg[2], gS, x, Ebin[2])
    YBaryon = exp(W_previous[3])
    YBaryonEquilibrium = YeqBaryon(g_deg[3], gS, x, Ebin[3])

    der = 0 #Initialise final output

    if procs[2] #Process No. 1
        der -= generic_derivative_of_process(sigma_RA, YDiQuark, YDiQuark, 1., 1., YDiQuarkEquilibrium^2, [true, true, false, false])
    end

    if procs[4] #Process No. 4 (quark-capture into diquark)
        eq_ratio4 = f4(x, gS, Ebin[2])
        der += generic_derivative_of_process(sigma_cap, YQuark, YQuark, YDiQuark, 1, eq_ratio4, [false, false, true, false])
    end

    if procs[5] #Process No. 5 (quark-diquark capture into baryons)
        eq_ratio5 = f5(x, gS, Ebin[2], Ebin[3])
        der -= generic_derivative_of_process(sigma_cap, YQuark, YDiQuark, YBaryon, 1., eq_ratio5, [false, true, false, false])
    end
    
    return der
end

function D_rhs_DiQuark_Baryon(x, cs, W_previous, g_deg, gS, Ebin, procs) # !!! UNDER CONSTRUCTION !!! Returns the derivative of the rhs of the quark equation w.r.t. the quark field
    #First define relevant quantities
    sigma_ann = cs[1]
    sigma_cap = cs[2]
    sigma_RA = cs[3]
    YQuark = exp(W_previous[1])
    YQuarkEquilibrium = YeqQuark(g_deg[1], gS, x)
    YDiQuark = exp(W_previous[2])
    YDiQuarkEquilibrium = YeqDiQuark(g_deg[2], gS, x, Ebin[2])
    YBaryon = exp(W_previous[3])
    YBaryonEquilibrium = YeqBaryon(g_deg[3], gS, x, Ebin[3])

    der = 0 #Initialise final output

    if procs[5] #Process No. 5 (quark-diquark capture into baryons)
        eq_ratio5 = f5(x, gS, Ebin[2], Ebin[3])
        der -= generic_derivative_of_process(sigma_cap, YQuark, YDiQuark, YBaryon, 1., eq_ratio5, [false, false, true, false])
    end
    
    return der
end

function D_rhs_Baryon_Quark(x, cs, W_previous, g_deg, gS, Ebin, procs) # !!! UNDER CONSTRUCTION !!! Returns the derivative of the rhs of the quark equation w.r.t. the quark field
    #First define relevant quantities
    sigma_ann = cs[1]
    sigma_cap = cs[2]
    sigma_RA = cs[3]
    YQuark = exp(W_previous[1])
    YQuarkEquilibrium = YeqQuark(g_deg[1], gS, x)
    YDiQuark = exp(W_previous[2])
    YDiQuarkEquilibrium = YeqDiQuark(g_deg[2], gS, x, Ebin[2])
    YBaryon = exp(W_previous[3])
    YBaryonEquilibrium = YeqBaryon(g_deg[3], gS, x, Ebin[3])

    der = 0 #Initialise final output

    if procs[5] #Process No. 5 (quark-diquark capture into baryons)
        eq_ratio5 = f5(x, gS, Ebin[2], Ebin[3])
        der += generic_derivative_of_process(sigma_cap, YQuark, YDiQuark, YBaryon, 1., eq_ratio5, [true, false, false, false])
    end
    
    return der
end

function D_rhs_Baryon_DiQuark(x, cs, W_previous, g_deg, gS, Ebin, procs) # !!! UNDER CONSTRUCTION !!! Returns the derivative of the rhs of the quark equation w.r.t. the quark field
    #First define relevant quantities
    sigma_ann = cs[1]
    sigma_cap = cs[2]
    sigma_RA = cs[3]
    YQuark = exp(W_previous[1])
    YQuarkEquilibrium = YeqQuark(g_deg[1], gS, x)
    YDiQuark = exp(W_previous[2])
    YDiQuarkEquilibrium = YeqDiQuark(g_deg[2], gS, x, Ebin[2])
    YBaryon = exp(W_previous[3])
    YBaryonEquilibrium = YeqBaryon(g_deg[3], gS, x, Ebin[3])

    der = 0 #Initialise final output

    if procs[5] #Process No. 5 (quark-diquark capture into baryons)
        eq_ratio5 = f5(x, gS, Ebin[2], Ebin[3])
        der += generic_derivative_of_process(sigma_cap, YQuark, YDiQuark, YBaryon, 1., eq_ratio5, [false, true, false, false])
    end
    
    return der
end

function D_rhs_Baryon_Baryon(x, cs, W_previous, g_deg, gS, Ebin, procs) # !!! UNDER CONSTRUCTION !!! Returns the derivative of the rhs of the quark equation w.r.t. the quark field
    #First define relevant quantities
    sigma_ann = cs[1]
    sigma_cap = cs[2]
    sigma_RA = cs[3]
    YQuark = exp(W_previous[1])
    YQuarkEquilibrium = YeqQuark(g_deg[1], gS, x)
    YDiQuark = exp(W_previous[2])
    YDiQuarkEquilibrium = YeqDiQuark(g_deg[2], gS, x, Ebin[2])
    YBaryon = exp(W_previous[3])
    YBaryonEquilibrium = YeqBaryon(g_deg[3], gS, x, Ebin[3])

    der = 0 #Initialise final output

    if procs[3] #Process No. 1
        der -= generic_derivative_of_process(sigma_RA, YBaryon, YBaryon, 1., 1., YBaryonEquilibrium^2, [true, true, false, false])
    end
    
    if procs[5] #Process No. 5 (quark-diquark capture into baryons)
        eq_ratio5 = f5(x, gS, Ebin[2], Ebin[3])
        der += generic_derivative_of_process(sigma_cap, YQuark, YDiQuark, YBaryon, 1., eq_ratio5, [false, false, true, false])
    end

    return der
end

function J11(x, cs, W_previous, g_deg, gS, Ebin, procs) # This function calculates the Jacobian matrix element J11 in each Euler backward timestep
    #argument definitions as in function "JMatrix"
    return 1 - D_rhs_Quark_Quark(x, cs, W_previous, g_deg, gS, Ebin, procs)
end

function J12(x, cs, W_previous, g_deg, gS, Ebin, procs) # This function calculates the Jacobian matrix element J12 in each Euler backward timestep
    #argument definitions as in function "JMatrix"
    return -D_rhs_Quark_DiQuark(x, cs, W_previous, g_deg, gS, Ebin, procs)
end

function J13(x, cs, W_previous, g_deg, gS, Ebin, procs) # This function calculates the Jacobian matrix element J13 in each Euler backward timestep
    #argument definitions as in function "JMatrix"
    return -D_rhs_Quark_Baryon(x, cs, W_previous, g_deg, gS, Ebin, procs)
end

function J21(x, cs, W_previous, g_deg, gS, Ebin, procs) # This function calculates the Jacobian matrix element J21 in each Euler backward timestep
    #argument definitions as in function "JMatrix"
    return -D_rhs_DiQuark_Quark(x, cs, W_previous, g_deg, gS, Ebin, procs)
end

function J22(x, cs, W_previous, g_deg, gS, Ebin, procs) # This function calculates the Jacobian matrix element J22 in each Euler backward timestep
    #argument definitions as in function "JMatrix"
    return 1 - D_rhs_DiQuark_DiQuark(x, cs, W_previous, g_deg, gS, Ebin, procs)
end

function J23(x, cs, W_previous, g_deg, gS, Ebin, procs) # This function calculates the Jacobian matrix element J23 in each Euler backward timestep
    #argument definitions as in function "JMatrix"
    return -D_rhs_DiQuark_Baryon(x, cs, W_previous, g_deg, gS, Ebin, procs)
end

function J31(x, cs, W_previous, g_deg, gS, Ebin, procs) # This function calculates the Jacobian matrix element J31 in each Euler backward timestep
    #argument definitions as in function "JMatrix"
    return -D_rhs_Baryon_Quark(x, cs, W_previous, g_deg, gS, Ebin, procs)
end

function J32(x, cs, W_previous, g_deg, gS, Ebin, procs) # This function calculates the Jacobian matrix element J32 in each Euler backward timestep
    #argument definitions as in function "JMatrix"
    return -D_rhs_Baryon_DiQuark(x, cs, W_previous, g_deg, gS, Ebin, procs)
end

function J33(x, cs, W_previous, g_deg, gS, Ebin, procs) # This function calculates the Jacobian matrix element J33 in each Euler backward timestep
    #argument definitions as in function "JMatrix"
    return 1 - D_rhs_Baryon_Baryon(x, cs, W_previous, g_deg, gS, Ebin, procs)
end