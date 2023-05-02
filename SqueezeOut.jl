using QuadGK

function Landau_pole(mu, alpha, beta0) #This calculates the Landau pole for a theory with coefficient beta0 given a coupling alpha(mu) at another scale mu. 
    Lambda = mu*exp(-2*pi/(beta0*alpha))
end

function running_coupling_from_pole(Q, Lambda, beta0) #Running coupling at the scale Q given a Landau pole Lambda with a theory with a beta0 coefficient
    return 2*pi/(beta0*log(Q/Lambda))
end

function running_coupling_from_scale(Q, mu, alphamu, beta0) #Running coupling at the scale Q given the running alphamu at another scale mu with a beta0 coefficient. 
    return alphamu/(1 + alphamu*beta0*log(Q/mu)/(2*pi))
end

function T_convert_u(u, ui, Ti) #Converts the variable u used in the g-average into temperatures
    T = Ti*(ui/u)^(2/3)
end

function g_average(ui, Ti, uf) #Calculates the averaged eff dofs of entropy needed in the entropy dilution
    integral_denom, err_denom = quadgk(u->u^(2/3)*exp(-u), ui, uf, rtol=1e-8)
    integral_num, err_num = quadgk(u->u^(2/3)*exp(-u)*h_eff_dof(T_convert_u(u, ui, Ti))^(1/3), ui, uf, rtol=1e-8)
    return integral_num/integral_denom
end

function rolled_index(big_ind, num_array) #This function converts a rolled index ind into the subindices small_inds = (n1, n2, ...). The maximally allowed values for the ni are stored in num_array
    # The formula that holds is 
    # ind = n1 + (n2-1)*L1 + (n3-1)*L1*L2 + ... + (nd-1)*L1*L2*...*L_(d-1)
    # Via basic arithmetics, this relation can be inverted for the ni. 
    d = length(num_array)
    small_inds = ones(d)
    
    new_array = deepcopy(num_array)
    new_ind = deepcopy(big_ind)
    
    h1 = prod(new_array) #It works, believe me.
    for j in d:-1:1 
        if new_ind == h1 || new_ind == 0
            for k in 1:j
                small_inds[k] = num_array[k]
            end
            break
        elseif j == 1
            small_inds[j] = new_ind
        else
            pop!(new_array)
            h1 = prod(new_array)
            h2 = ceil(new_ind/h1) 
            if new_ind > h1
               small_inds[j] = h2
               new_ind -= h1*(h2-1)
            end
        end
    end
    return Int.(small_inds)
end

function pocket_radius(Lambda) # returns pocket radius in units 1/Lambda
    R0 = 1E-6*(Lambda/Mpl)^(-0.9)
    R1 = (Mpl/(1E4*Lambda))^(2/3)
    return max(R0, R1)
end

function aux_func_GB_decay(xfo, m_GB, R)
    sigma32 = (4*pi)^3/Ndark^6/m_GB
    lhs = log(xfo)*(5/2) + 2*xfo - log(h_eff_dof(m_GB/xfo)*R/(180*pi)*(Mpl*sigma32/sqrt(4/5*pi^3*g_eff_dof(m_GB/xfo)))^(3/2))
    return lhs
end

function GB_freeze_out_estimate(xft, m_GB, R) #Numerical estimate of the freeze-out x with the secant method. xft is the initial guess, m is the DM mass and R is the entropy ratio. acs is the constant annihilation xs. cs_xvalues are the x_values of the fully dressed acs. cs_yvalues are the corresponding acs values. coeffs are the spline interpolation coefficients.
    xf0 = xft #first try
    xf1 = xf0 - aux_func_GB_decay(xf0, m_GB, R)*2*0.001/(aux_func_GB_decay(xf0 + 0.001, m_GB, R)-aux_func_GB_decay(xf0 - 0.001, m_GB, R))
    diff = abs(xf1 - xf0)
    if diff < 1E-4
        xf2 = xf1
    else
        while diff > 1E-4
            xf2 = xf1 - aux_func_GB_decay(xf1, m_GB, R)*(xf1 - xf0)/(aux_func_GB_decay(xf1, m_GB, R) - aux_func_GB_decay(xf0, m_GB, R))
            diff = abs(xf2 - xf1)
            xf0 = copy(xf1)
            xf1 = copy(xf2)
        end
    end
    return xf2
end