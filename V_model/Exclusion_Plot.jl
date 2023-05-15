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

#Define the limits according to physics
const Relic_abundance_limit = 0.12 
const BBN_lifetime = 6.58*1E-25 #Lower bound on glueball decay rate.

test_df = DataFrame(CSV.File("V_model_scan.csv")) #Read in data as a dataframe

#Initialise the data in respective arrays
y_values = test_df[!, 2] #values for Lambda
x_values = test_df[!, 1]./test_df[!, 2] #values for m/Lambda
GB_rate_values = test_df[!, 10]
relic_abundance_values = test_df[!, 12] 

#Data for scatter plot of excluded regions
excluded_relic_indices = findall(relic_abundance_values.>Relic_abundance_limit) #Find indices of data which are EXCLUDED by observations
excluded_BBN_indices = findall(GB_rate_values.<BBN_lifetime)
x_relic = x_values[excluded_relic_indices] #Filter the data to obtain only the EXCLUDED regions
y_relic = y_values[excluded_relic_indices]
x_BBN = x_values[excluded_BBN_indices]
y_BBN = y_values[excluded_BBN_indices]

#Scatter plot to get a feeling for the result 
scatter(x_relic, y_relic, title = "V model", minorgrid = true, minorticks = 10, legend=:bottomleft, xscale = :log10, yscale = :log10, ylabel = L"\Lambda/GeV", xlabel = L"m_Q/\Lambda", labels = L"\Omega_{DM} h^2>0.12")
scatter!(x_BBN, y_BBN, labels = L"\tau_{GB} > 1s")
savefig("Exclusion_plot_Vmodel_scatter_plot.png")

x_relic_line, y_relic_line = line_finder(x_relic, y_relic, 10, 0)
x_BBN_line, y_BBN_line = line_finder(x_BBN, y_BBN, 10, 1)

plot(x_BBN_line, y_BBN_line, fillrange = 1*ones(size(x_BBN_line)), fillalpha = 0.35, c = 2, xlims = (1E2, 1E4), ylims = (1, 1E7), title="V model", minorticks = 10,  minorgrid = true, legend=false, ylabel = L"\Lambda/GeV", xlabel = L"m_Q/\Lambda", xscale = :log10, yscale = :log10)
plot!(x_relic_line, y_relic_line, fillrange = 1E7*ones(size(x_relic_line)), fillalpha = 0.35, c = 1)
annotate!(7000,5, text(L"\tau_{GB} > 1s", :black, :right, 12))
annotate!(300, 1E6, text(L"\Omega_{DM} h^2>0.12", :black, :right, 12))
savefig("Exclusion_plot_Vmodel_line_plot.png")