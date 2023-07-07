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
const Relic_abundance_limit = 0.12 
const BBN_lifetime = 1.52*1E24 #Upper bound on glueball decay rate in Gev^-1.

test_df = DataFrame(CSV.File("V_model_scan.csv")) #Read in data as a dataframe

#Initialise the data in respective arrays
y_values = test_df[!, 2] #values for Lambda
x_values = test_df[!, 1]./test_df[!, 2] #values for m/Lambda
GB_rate_values = test_df[!, 11]
relic_abundance_values = test_df[!, 13] 

Yfo_values = test_df[!, 7] #Yields after quark freeze-out 
Ysqo_values = test_df[!, 8] #Yields after 1st order phase transition
Dil_fac_values = test_df[!, 12] #Entropy dilution factor

relic_abundance_fo_only_values = Yfo_values .* relic_abundance_values ./ Dil_fac_values ./ Ysqo_values
relic_abundance_fo_sq_values = relic_abundance_values./Dil_fac_values
relic_abundance_fo_dil_values = Yfo_values.*relic_abundance_values./Ysqo_values

#Data for scatter plot of excluded regions: Freeze-Out only
excluded_relic_indices_fo_only = findall(relic_abundance_fo_only_values.>Relic_abundance_limit) #Find indices of data which are EXCLUDED by observations
x_relic_fo_only = x_values[excluded_relic_indices_fo_only] #Filter the data to obtain only the EXCLUDED regions
y_relic_fo_only = y_values[excluded_relic_indices_fo_only]

#Data for scatter plot of excluded regions: Freeze-Out + Squeezeout
excluded_relic_indices_fo_sq= findall(relic_abundance_fo_sq_values.>Relic_abundance_limit) #Find indices of data which are EXCLUDED by observations
x_relic_fo_sq = x_values[excluded_relic_indices_fo_sq] #Filter the data to obtain only the EXCLUDED regions
y_relic_fo_sq = y_values[excluded_relic_indices_fo_sq]

#Data for scatter plot of excluded regions: Freeze-Out + Entropy dilution
excluded_relic_indices_fo_dil = findall(relic_abundance_fo_dil_values.>Relic_abundance_limit) #Find indices of data which are EXCLUDED by observations
x_relic_fo_dil = x_values[excluded_relic_indices_fo_dil] #Filter the data to obtain only the EXCLUDED regions
y_relic_fo_dil = y_values[excluded_relic_indices_fo_dil]

#Data for scatter plot of excluded regions: All processes
excluded_relic_indices = findall(relic_abundance_values.>Relic_abundance_limit) #Find indices of data which are EXCLUDED by observations
x_relic = x_values[excluded_relic_indices] #Filter the data to obtain only the EXCLUDED regions
y_relic = y_values[excluded_relic_indices]
excluded_BBN_indices = BBN_time(y_values, GB_rate_values, BBN_lifetime)
x_BBN = x_values[excluded_BBN_indices]
y_BBN = y_values[excluded_BBN_indices]

#Scatter plot to get a feeling for the result
scatter(x_relic_fo_only, y_relic_fo_only, title = "V model", minorgrid = true, minorticks = 10, legend=:bottomleft,
 xscale = :log10, yscale = :log10, ylabel = L"\Lambda/GeV", xlabel = L"m_Q/\Lambda", labels = L"\Omega_{DM} h^2>0.12",
 xlims = (25, 1E4), ylims = (1E-3, 1E7))
scatter!(x_BBN, y_BBN, labels = L"\tau_{GB} > 1s")
savefig("Exclusion_plot_Vmodel_scatter_relic_only.png")

x_relic_line_fo_only, y_relic_line_fo_only = line_finder(x_relic_fo_only, y_relic_fo_only, 10, 0)
x_relic_line_fo_sq, y_relic_line_fo_sq = line_finder(x_relic_fo_sq, y_relic_fo_sq, 10, 0)
x_relic_line_fo_dil, y_relic_line_fo_dil = line_finder(x_relic_fo_dil, y_relic_fo_dil, 10, 0)
x_relic_line, y_relic_line = line_finder(x_relic, y_relic, 10, 0)
x_BBN_line, y_BBN_line = line_finder(x_BBN, y_BBN, 10, 1)

plot(x_relic_line, y_relic_line, fillrange = 1E7*ones(size(x_relic_line)), fillalpha = 0.35, c = 1, xlims = (25, 1E4), ylims = (1E-3, 1E7), title="V model",
minorticks = 10,  minorgrid = true, legend=false, ylabel = L"\Lambda/GeV", xlabel = L"m_Q/\Lambda", xscale = :log10, yscale = :log10)
plot!(x_relic_line_fo_only, y_relic_line_fo_only, linestyle=:dash, seriescolor=:blue, fillrange = 1E7*ones(size(x_relic_line_fo_only)), fillalpha = 0.35, c = 1)
plot!(x_relic_line_fo_sq, y_relic_line_fo_sq, linestyle=:dash, seriescolor=:blue)
plot!(x_BBN_line, y_BBN_line, fillrange = 1E-3*ones(size(x_BBN_line)), fillalpha = 0.35, c = 2)
annotate!(7000, 0.005, text(L"\tau_{GB} > 1s", :black, :right, 12))
annotate!(300, 1E6, text(L"\Omega_{DM} h^2>0.12", :black, :right, 12))
annotate!(3000, 1E5, text("w/o entropy dilution", :black, :right, 12))
annotate!(800, 30, text("FO only (+- dilution)", :black, :right, 12))
savefig("Exclusion_plot_Vmodel_line_Testplot.png")

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