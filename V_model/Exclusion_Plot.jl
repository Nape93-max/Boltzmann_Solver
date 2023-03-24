using Plots
using LaTeXStrings
using CSV, DataFrames

#Define the limits according to physics
const Relic_abundance_limit = 0.12 
const BBN_lifetime = 1/1.52*1E-22 #Lower bound on glueball decay rate.

test_df = DataFrame(CSV.File("V_model_scan.csv")) #Read in data as a dataframe

#Initialise the data in respective arrays
x_values = test_df[!, 2] #values for Lambda
y_values = test_df[!, 3] #values for m/Lambda
GB_rate_values = test_df[!, 10]
relic_abundance_values = test_df[!, 12] 
correct_relic_indices = findall(relic_abundance_values.>Relic_abundance_limit) #Find indices of data which are EXCLUDED by observations
correct_BBN_indices = findall(GB_rate_values.<BBN_lifetime)

x_relic = x_values[correct_relic_indices] #Filter the data to obtain only the EXCLUDED regions
y_relic = y_values[correct_relic_indices]
x_BBN = x_values[correct_BBN_indices]
y_BBN = y_values[correct_BBN_indices]

#Scatter plot to get a feeling for the result 
scatter(x_relic, y_relic, title = "V model", minorgrid = true, minorticks = 10, legend=:bottomleft, xscale = :log, yscale = :log, xlabel = L"\Lambda/GeV", ylabel = L"m_Q/\Lambda", labels="Overclosure")
scatter!(x_BBN, y_BBN, labels="BBN constraint on GB lifetime")
savefig("Exclusion_plot_Vmodel_scatter_plot.png")