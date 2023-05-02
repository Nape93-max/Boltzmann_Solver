include("FreezeOut.jl")

g_quark = 4 #degeneracy of the Dirac quark

m_quark = 1000
Alpha_DM = 0.01

BigConstant = bc_constant(m_quark)
sigma0 = pert_acs(Alpha_DM, m_quark)

Y_final = quark_freeze_out(1000, m_quark, sigma0, BigConstant, g_quark)

Omega_relic = reduced_Hubble_squared*Y_final*s0*m_quark/rho_crit
println(Omega_relic, "\t", Y_final)