using Plots
using LaTeXStrings
using CSV, DataFrames

include("../FreezeOut.jl") #Include important functions regarding freeze-out

function line_finder(x_data, y_data, error, minormax)
    xline = sort(unique(x_data))
    ylineMIN = ones(length(xline))
    ylineMAX = ones(length(xline))

    for i in 1:length(xline)
       ylineMIN[i] = minimum(y_data[abs.(xline[i] .- x_data) .< error])
       ylineMAX[i] = maximum(y_data[abs.(xline[i] .- x_data) .< error])
    end

    if minormax == 0
        return xline, ylineMIN
    else
        return xline, ylineMAX
    end
end

function BBN_time(Lambda_data, GB_data, LTLimit) #function that returns the indices of those measurements
    # with Lambda_data and GB_data which have a lifetime ABOVE LTLimit and are thereby EXCLUDED
    tau_GB = 1 ./ GB_data
    g_eff_sqrt_vec = eff_dof_sqrt.(Lambda_data) 
    t_Lambda = 0.75/pi*sqrt(5/pi)*Mpl ./ (Lambda_data .^ 2) ./ g_eff_sqrt_vec
    t_total = tau_GB + t_Lambda
    return findall(t_total .> LTLimit)
end

#Define the limits according to physics
Ndark = 3
CF = 4/3
Relic_abundance_limit = 0.12
BBN_lifetime = 1.52*1E24 #Upper bound on glueball decay rate in Gev^-1.
Bullet_Cluster_velocity = 0.015 #vrel = 4500 km s^-1 in units v/c
Bullet_Cluster_Limit = 4638 #Limit on sigma/m \sim 1 cm^2 g^-1 in GeV^-3
Positron_mass_limit = 350000 #Upper mass bound found from AMS-02 positron spectra for 100% of DM 

feff = 0.2 #CMB feff parameter for annihilation into W bosons (fig. 4 of 1506.03811)
CMB_bound_W_bosons = 3.495*1E-11 #Upper CMB bound on feff*<sigma v>/m_B in GeV^(-3) 

test_df = DataFrame(CSV.File("V_model_scan.csv")) #Read in data as a dataframe

#Initialise the data in respective arrays
y_values = test_df[!, 2] #values for Lambda
x_values = test_df[!, 1]./test_df[!, 2] #values for m/Lambda
mass_values = test_df[!, 1] #values for the DM mass
coupling_values = test_df[!, 4] #values for the DM coupling alpha(m) (should be replaced by alpha_B, but who cares at this point?)
GB_rate_values = test_df[!, 11]
relic_abundance_values = test_df[!, 13] 

Yfo_values = test_df[!, 7] #Yields after quark freeze-out 
Ysqo_values = test_df[!, 8] #Yields after 1st order phase transition
Dil_fac_values = test_df[!, 12] #Entropy dilution factor

DM_self_scattering_values = 4*pi*sqrt(2)./(CF*Bullet_Cluster_velocity*sqrt(3)*coupling_values.*mass_values.^3)

CMB_baryon_values = feff*sqrt(8/Ndark)*pi./(CF*coupling_values.*mass_values.^3) #Values for feff*<sigma v>/m_B in GeV^(-3) 

# sigma/m for the rearrangement cross section in GeV^-3

#=

relic_abundance_fo_dil_values = Yfo_values.*relic_abundance_values./Ysqo_values
relic_abundance_fo_only_values = Yfo_values .* relic_abundance_values ./ Dil_fac_values ./ Ysqo_values
relic_abundance_fo_sq_values = relic_abundance_values./Dil_fac_values


#Data for scatter plot of excluded regions: Freeze-Out only
excluded_relic_indices_fo_only = findall(relic_abundance_fo_only_values.>Relic_abundance_limit) #Find indices of data which are EXCLUDED by observations
x_relic_fo_only = x_values[excluded_relic_indices_fo_only] #Filter the data to obtain only the EXCLUDED regions
y_relic_fo_only = y_values[excluded_relic_indices_fo_only]

#Data for scatter plot of excluded regions: Freeze-Out + Squeezeout
excluded_relic_indices_fo_sq= findall(relic_abundance_fo_sq_values.>Relic_abundance_limit) #Find indices of data which are EXCLUDED by observations
x_relic_fo_sq = x_values[excluded_relic_indices_fo_sq] #Filter the data to obtain only the EXCLUDED regions
y_relic_fo_sq = y_values[excluded_relic_indices_fo_sq]
=#
#Data for scatter plot of excluded regions: Freeze-Out + Entropy dilution

excluded_Bullet_cluster = findall(DM_self_scattering_values.>Bullet_Cluster_Limit) #Find indices of data which are EXCLUDED 
# by Bullet cluster observations.
x_Bullet = x_values[excluded_Bullet_cluster] #Filter the data to obtain only the EXCLUDED regions
y_Bullet = y_values[excluded_Bullet_cluster]

excluded_Positron_spectra = findall((mass_values*Relic_abundance_limit./relic_abundance_values.<Positron_mass_limit)
 .& (Relic_abundance_limit.>relic_abundance_values)) #Find indices of data which are EXCLUDED by AMS-02 observations
x_Positron = x_values[excluded_Positron_spectra]
y_Positron = y_values[excluded_Positron_spectra]

excluded_CMB_values = findall(CMB_baryon_values .> CMB_bound_W_bosons)
x_CMB = x_values[excluded_CMB_values]
y_CMB = y_values[excluded_CMB_values]

#=
excluded_relic_indices_fo_dil = findall(relic_abundance_fo_dil_values.>Relic_abundance_limit) #Find indices of data which are EXCLUDED by observations
x_relic_fo_dil = x_values[excluded_relic_indices_fo_dil] #Filter the data to obtain only the EXCLUDED regions
y_relic_fo_dil = y_values[excluded_relic_indices_fo_dil]
=#

#Data for scatter plot of excluded regions: All processes
excluded_relic_indices = findall(relic_abundance_values.>Relic_abundance_limit) #Find indices of data which are EXCLUDED by observations
x_relic = x_values[excluded_relic_indices] #Filter the data to obtain only the EXCLUDED regions
y_relic = y_values[excluded_relic_indices]
excluded_BBN_indices = BBN_time(y_values, GB_rate_values, BBN_lifetime)
x_BBN = x_values[excluded_BBN_indices]
y_BBN = y_values[excluded_BBN_indices]

#Scatter plot to get a feeling for the result
scatter(x_relic, y_relic, title = "V model", minorgrid = true, minorticks = 10, legend=:bottomleft,
 xscale = :log10, yscale = :log10, ylabel = L"\Lambda/GeV", xlabel = L"m_Q/\Lambda", labels = L"\Omega_{DM} h^2>0.12",
 xlims = (2, 1E4), ylims = (1E-3, 1E7))
scatter!(x_BBN, y_BBN, labels = L"\tau_{GB} > 1s")
scatter!(x_CMB, y_CMB, labels = "CMB")
#scatter!(x_Bullet, y_Bullet, labels = L"\sigma/m < 1cm^2 g^{-1}")
#scatter!(x_Positron, y_Positron, labels = L"e^+ AMS-02")
savefig("Exclusion_plot_Vmodel_scatter.png")
#=
x_relic_line_fo_only, y_relic_line_fo_only = line_finder(x_relic_fo_only, y_relic_fo_only, 10, 0)
x_relic_line_fo_sq, y_relic_line_fo_sq = line_finder(x_relic_fo_sq, y_relic_fo_sq, 10, 0)
x_relic_line_fo_dil, y_relic_line_fo_dil = line_finder(x_relic_fo_dil, y_relic_fo_dil, 10, 0)
=#

x_relic_line, y_relic_line = line_finder(x_relic, y_relic, 0.1, 0)
x_BBN_line, y_BBN_line = line_finder(x_BBN, y_BBN, 0.1, 1)
x_Bullet_line, y_Bullet_line = line_finder(x_Bullet, y_Bullet, 0.1, 1)
x_Positron_line, y_Positron_line = line_finder(x_Positron, y_Positron, 0.1, 1)
x_CMB_line, y_CMB_line = line_finder(x_CMB, y_CMB, 0.1, 1)

x_mass1_line = 10.0.^collect(range(0.3, 4, length = 100))
y_mass1_line = (1*10^0)./x_mass1_line

x_mass2_line = 10.0.^collect(range(0.3, 4, length = 100))
y_mass2_line = (1*10^3)./x_mass2_line

x_mass3_line = 10.0.^collect(range(0.3, 4, length = 100))
y_mass3_line = (5*10^6)./x_mass3_line

x_mass4_line = 10.0.^collect(range(0.3, 4, length = 100))
y_mass4_line = (1*10^9)./x_mass4_line

x_grav_decay_mass_line = 10.0.^collect(range(0.3, 4, length = 100))
y_grav_decay_mass_line = (1.7*10^10)./x_grav_decay_mass_line

plot(x_relic_line, y_relic_line, fillrange = 1E7*ones(size(x_relic_line)), fillalpha = 0.35, c = 1, xlims = (2, 1E4), ylims = (1E-3, 1E7),
 legend=false, ylabel = L"\Lambda/GeV", xlabel = L"m_Q/\Lambda", xscale = :log10, yscale = :log10,
  minorgrid = true, xticks=5, yticks = 5, minorticks=10)
#plot!(x_relic_line_fo_dil, y_relic_line_fo_dil, linestyle=:dash, seriescolor=:blue, linewidth = 3)
#plot!(x_relic_line_fo_sq, y_relic_line_fo_sq, linestyle=:dash, seriescolor=:blue)
plot!(x_BBN_line, y_BBN_line, fillrange = 1E-3*ones(size(x_BBN_line)), fillalpha = 0.35, c = 2)
plot!(x_CMB_line, y_CMB_line, fillrange = 1E-3*ones(size(x_CMB_line)), fillalpha = 0.35, c = 3)
#plot!(x_Bullet_line, y_Bullet_line, fillrange = 1E-3*ones(size(x_Bullet_line)), fillalpha = 0.35, c = 3)
plot!(x_grav_decay_mass_line, y_grav_decay_mass_line, fillrange = 1E7*ones(size(x_grav_decay_mass_line)), fillalpha=1.0, c=:black)
plot!(x_mass1_line, y_mass1_line, c=:black, linestyle=:dash)
annotate!(3, 1, text(L"m_Q = 1"*" GeV", :black, :left, 12))
annotate!(500, 0.7, text(L"m_Q = 1"*" TeV", :black, :right, 12))
annotate!(1000, 2000, text(L"m_Q = 10^3"*" TeV", :black, :right, 12))
annotate!(3000, 8*10^4, text(L"m_Q = 10^6"*" TeV", :black, :right, 12))
#annotate!(8000, 5, text("No squeezeout", :blue, :right, 12))
plot!(x_mass2_line, y_mass2_line, c=:black, linestyle=:dash)
plot!(x_mass3_line, y_mass3_line, c=:black, linestyle=:dash)
plot!(x_mass4_line, y_mass4_line, c=:black, linestyle=:dash)
annotate!(7000, 10, text(L"\tau_{GB} > 1s", :red, :right, 12))
annotate!(15, 1E5, text(L"\Omega_{DM} h^2>0.12", :blue, :right, 12))
annotate!(30, 0.01, text("CMB: "*L"f_{eff}\frac{<\sigma v>}{m_B}", :green, :right, 12))
#annotate!(3000, 1E5, text("w/o entropy dilution", :black, :right, 12))
#annotate!(800, 30, text("FO only (+- dilution)", :black, :right, 12))
savefig("Exclusion_plot_Vmodel_line.png")

#=
#Scatter plot to get a feeling for the result 
scatter(x_relic, y_relic, title = "V model", minorgrid = true, minorticks = 10, legend=:bottomleft, xscale = :log10, yscale = :log10, ylabel = L"\Lambda/GeV", xlabel = L"m_Q/\Lambda", labels = L"\Omega_{DM} h^2>0.12")
scatter!(x_BBN, y_BBN, labels = L"\tau_{GB} > 1s")
savefig("Exclusion_plot_Vmodel_scatter_plot.png")

plot(x_BBN_line, y_BBN_line, fillrange = 1*ones(size(x_BBN_line)), fillalpha = 0.35, c = 2, xlims = (1E2, 1E4), ylims = (1, 1E7), title="V model", minorticks = 10,  minorgrid = true, legend=false, ylabel = L"\Lambda/GeV", xlabel = L"m_Q/\Lambda", xscale = :log10, yscale = :log10)
plot!(x_relic_line, y_relic_line, fillrange = 1E7*ones(size(x_relic_line)), fillalpha = 0.35, c = 1)
annotate!(7000,5, text(L"\tau_{GB} > 1s", :black, :right, 12))
annotate!(300, 1E6, text(L"\Omega_{DM} h^2>0.12", :black, :right, 12))
savefig("Exclusion_plot_Vmodel_line_plot.png")
=#