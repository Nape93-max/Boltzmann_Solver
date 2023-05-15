using Plots
using LaTeXStrings
using CSV, DataFrames

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

function excluded_relic_inds(relic_values, delta_values, delta, Yukawa_values, Yukawa)
    findall((relic_values.>Relic_abundance_limit) .& (abs.(delta_values .- delta) .< 1E-1) .& ((abs.(Yukawa_values .- Yukawa) .< 1E-8))) #Find indices of data which are EXCLUDED by observations
end

function excluded_BBN_inds(Gamma_values, delta_values, delta, Yukawa_values, Yukawa)
    findall((Gamma_values .< BBN_lifetime) .& (abs.(delta_values .- delta) .< 1E-1) .& ((abs.(Yukawa_values .- Yukawa) .< 1E-8))) #Find indices of data which are EXCLUDED by observations
end

function relic_line(x_data, y_data, relic_abundance_values, delta_values, delta, Yukawa_values, Yukawa)
    relic_inds = excluded_relic_inds(relic_abundance_values, delta_values, delta, Yukawa_values, Yukawa)
    x_relic = x_data[relic_inds]
    y_relic = y_data[relic_inds]
    x_relic_line, y_relic_line = line_finder(x_relic, y_relic, 10, 0)
    return x_relic_line, y_relic_line
end

function BBN_line(x_data, y_data, BBN_values, delta_values, delta, Yukawa_values, Yukawa)
    BBN_inds = excluded_BBN_inds(BBN_values, delta_values, delta, Yukawa_values, Yukawa)
    x_BBN = x_data[BBN_inds]
    y_BBN = y_data[BBN_inds]
    x_BBN_line, y_BBN_line = line_finder(x_BBN, y_BBN, 10, 1)
    return x_BBN_line, y_BBN_line
end

#Define the limits according to physics
const Relic_abundance_limit = 0.12 
const BBN_lifetime = 6.58*1E-25 #Lower bound on glueball decay rate.

test_df = DataFrame(CSV.File("NplusL_model_scan.csv")) #Read in data as a dataframe

Yukawa_coupling = 1E-1 #For simple plots
#Mass_delta = 12.9154966

#Initialise the data in respective arrays
y_values = test_df[!, 2] #values for Lambda
x_values = test_df[!, 1] ./ test_df[!, 2] #values for m/Lambda
delta_values = test_df[!, 5] #values for delta
Yukawa_values = test_df[!, 7] #values for yd
GB_rate_values = test_df[!, 13]
relic_abundance_values = test_df[!, 15] 

#BLOCK for colour maps


#BLOCK for the simple plots
#=
correct_relic_indices = excluded_relic_inds(relic_abundance_values, delta_values, Mass_delta, Yukawa_values, Yukawa_coupling) #Find indices of data which are EXCLUDED by observations
correct_BBN_indices = excluded_BBN_inds(GB_rate_values, delta_values, Mass_delta,  Yukawa_values, Yukawa_coupling)

x_relic = x_values[correct_relic_indices] #Filter the data to obtain only the EXCLUDED regions
y_relic = y_values[correct_relic_indices]
x_BBN = x_values[correct_BBN_indices]
y_BBN = y_values[correct_BBN_indices]

#Scatter plot to get a feeling for the result 
scatter(x_relic, y_relic, title = "N + L model. yd = $Yukawa_coupling, delta = $Mass_delta", minorgrid = true, minorticks = 10, xlims = (100, 1E4), ylims = (1, 1E7), legend=:bottomleft, xscale = :log10, yscale = :log10, xlabel = L"\Lambda/GeV", ylabel = L"m_N/\Lambda", labels="Overclosure")
scatter!(x_BBN, y_BBN, labels="BBN constraint on GB lifetime")
savefig("Plots/Exclusion_plot_NplusLmodel_scatter_plot_No9.png")

x_relic_line, y_relic_line = relic_line(x_values, y_values, relic_abundance_values, delta_values, Mass_delta, Yukawa_values,  Yukawa_coupling)
x_BBN_line, y_BBN_line = BBN_line(x_values, y_values, GB_rate_values, delta_values, Mass_delta, Yukawa_values,  Yukawa_coupling)

plot(x_BBN_line, y_BBN_line, fillrange = 1*ones(size(x_BBN_line)), fillalpha = 0.35, c = 2, xlims = (1E2, 1E4), ylims = (1, 1E7), title="N + L model. yd = $Yukawa_coupling, delta = $Mass_delta", minorticks = 10,  minorgrid = true, legend=false, ylabel = L"\Lambda/GeV", xlabel = L"m_N/\Lambda", xscale = :log10, yscale = :log10)
plot!(x_relic_line, y_relic_line, fillrange = 1E7*ones(size(x_relic_line)), fillalpha = 0.35, c = 1)
if !isempty(x_BBN_line)
    annotate!(7000,5, text(L"\tau_{GB} > 1s", :black, :right, 12))
end
if !isempty(x_relic_line)
    annotate!(275, 5E6, text(L"\Omega_{DM} h^2>0.12", :black, :right, 12))
end
savefig("Plots/Exclusion_plot_NplusLmodel_line_plot_No9.png")
=#
