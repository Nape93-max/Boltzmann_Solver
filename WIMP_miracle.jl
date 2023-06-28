include("FreezeOut.jl")
using Plots
using GR
using LaTeXStrings
#using PlotlyJS

g_quark = 4 #degeneracy of the Dirac quark
#=
Lambda = 1E2
ratio = 1000
m_quark = ratio*Lambda
Alpha_DM = running_coupling_from_pole(m_quark, Lambda, 11)
=#
### BLOCK FOR HEATMAP GENERATION
len_scales = 10
len_masses = 10
array_scales = 10.0.^collect(range(-3, 7, len_scales)) 
array_masses = 10.0.^collect(range(1.4, 7, len_masses)) 
xi_hm = ones(len_scales, len_masses)
for i = 1:len_scales
    for j = 1:len_masses
        Lambda = array_masses[i]
        m_quark = array_masses[j]*Lambda
        Alpha_DM = running_coupling_from_pole(m_quark, Lambda, 11)

        BigConstant = bc_constant(m_quark)
        sigma0 = pert_acs(Alpha_DM, m_quark)

        Y_final = quark_freeze_out(Lambda, m_quark, sigma0, BigConstant, g_quark)
        xi_hm[i,j] = Alpha_DM*m_quark/(Lambda*cbrt(Y_final*2*pi*pi/45*h_eff_dof(Lambda)))
        println(xi_hm[i,j])
    end
end
Plots.heatmap(array_masses, array_scales, xi_hm, 
c=cgrad([:blue, :white,:red, :yellow], scale=:log), ylabel=L"\Lambda / GeV", xlabel=L"m_Q/\Lambda",
title="Inter quark spacing "*L"\xi", xaxis=:log10, yaxis=:log10)
Plots.savefig("Inter_quark_spacing.png")


#=
BigConstant = bc_constant(m_quark)
sigma0 = pert_acs(Alpha_DM, m_quark)

Y_final = quark_freeze_out(Lambda, m_quark, sigma0, BigConstant, g_quark)

Omega_relic = reduced_Hubble_squared*Y_final*s0*m_quark/rho_crit
inter_quark_spacing = 1/cbrt(Y_final*2*pi*pi/45*h_eff_dof(Lambda))
println(Omega_relic, "\t", Y_final)
=#


#=
### Here the plotting business starts ###

ytics = 10.0.^collect(range(-30, -1, step=1))
xtics = 10.0.^collect(range(log10(x_initial), log10(x_today), step=1))

x_vec_squeezeout = range(x_PT, length=100, stop=x_dilution)
y_vec_squeezeout = Yx_squeezeout*ones(100)
x_vec_dilution = range(x_dilution, length=100, stop=x_today)
y_vec_dilution = Yx_dilution*ones(100)

#plot(sigma_v_x_values, sigma_v_y_values, xaxis=:log, yaxis=:log)
#savefig("sigma_v_data.png")

plot(xvec, [EquilibriumYield, Yx], title="WIMP freeze-out", label=[L"Y_{eq}(x)" L"Y(x)"], yticks = ytics, xticks = xtics, minorticks = 10, minorgrid = true, xlabel="x = m/T", ylabel="Y(x)", xaxis=:log, yaxis=:log, xlims = (x_initial, x_final*1000), ylims = (1E-30, 1E-1))
plot!(x_vec_squeezeout, y_vec_squeezeout, label = L"Y_S")
plot!(x_vec_dilution, y_vec_dilution, label = L"Y_\infty")
plot!([x_MR], seriestype = :vline, label = L"x_{MR}")
plot!([x_PT], seriestype = :vline, label = L"x_\Lambda")
plot!([x_dilution], seriestype = :vline, label = L"x_{GB}")
savefig("FreezeOut.png")

plot(xvec, Yx, title="Relic Yield", minorticks = 10, minorgrid = true, xlabel=L"x = m/T", ylabel=L"Y(x)", label = L"Y(x)", yticks = ytics, xaxis=:log, yaxis=:log, xticks = xtics, xlims = (x_initial, x_final), ylims = (1E-16, 1E-1))
savefig("Y_plot.png")
=#