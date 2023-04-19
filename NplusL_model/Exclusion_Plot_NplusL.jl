using Plots
using LaTeXStrings
using CSV, DataFrames

#Define the limits according to physics
const Relic_abundance_limit = 0.12 
const BBN_lifetime = 6.58*1E-25 #Lower bound on glueball decay rate.

test_df = DataFrame(CSV.File("NplusL_model_Full_scan.csv")) #Read in data as a dataframe

Yukawa_coupling = 0.0483293023857175
Mass_delta = 29.7635144163132

#Initialise the data in respective arrays
x_values = test_df[!, 2] #values for Lambda
y_values = test_df[!, 3] #values for m/Lambda
delta_values = test_df[!, 5] #values for delta
Yukawa_values = test_df[!, 7] #values for yd
GB_rate_values = test_df[!, 13]
relic_abundance_values = test_df[!, 15] 
correct_relic_indices = findall((relic_abundance_values.>Relic_abundance_limit) .& (abs.(delta_values .- Mass_delta) .< 1E-8) .& ((abs.(Yukawa_values .- Yukawa_coupling) .< 1E-8))) #Find indices of data which are EXCLUDED by observations
correct_BBN_indices = findall((GB_rate_values.<BBN_lifetime) .& (abs.(delta_values .- Mass_delta) .< 1E-8) .& ((abs.(Yukawa_values .- Yukawa_coupling) .< 1E-8)))

x_relic = x_values[correct_relic_indices] #Filter the data to obtain only the EXCLUDED regions
y_relic = y_values[correct_relic_indices]
x_BBN = x_values[correct_BBN_indices]
y_BBN = y_values[correct_BBN_indices]

#Scatter plot to get a feeling for the result 
scatter(x_relic, y_relic, title = "N + L model. yd = $Yukawa_coupling, delta = $Mass_delta", minorgrid = true, minorticks = 10, legend=:bottomleft, xscale = :log, yscale = :log, xlabel = L"\Lambda/GeV", ylabel = L"m_Q/\Lambda", labels="Overclosure")
scatter!(x_BBN, y_BBN, labels="BBN constraint on GB lifetime")
savefig("Plots/Exclusion_plot_NplusLmodel_scatter_plot_No9.png")