using Plots
using LaTeXStrings
using CSV, DataFrames

include("../FreezeOut.jl") #Include important functions

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

function excluded_BBN_inds(Gamma_values, Lambda_data, delta_values, delta, Yukawa_values, Yukawa)
    tau_GB = 1 ./ Gamma_values
    g_eff_sqrt_vec = eff_dof_sqrt.(Lambda_data) 
    t_Lambda = 0.75/pi*sqrt(5/pi)*Mpl ./ (Lambda_data .^ 2) ./ g_eff_sqrt_vec
    t_total = tau_GB + t_Lambda
    
    findall((t_total .> BBN_lifetime) .& (abs.(delta_values .- delta) .< 1E-1) .& ((abs.(Yukawa_values .- Yukawa) .< 1E-8))) #Find indices of data which are EXCLUDED by observations
end

function relic_line(x_data, y_data, relic_abundance_values, delta_values, delta, Yukawa_values, Yukawa)
    relic_inds = excluded_relic_inds(relic_abundance_values, delta_values, delta, Yukawa_values, Yukawa)
    x_relic = x_data[relic_inds]
    y_relic = y_data[relic_inds]
    x_relic_line, y_relic_line = line_finder(x_relic, y_relic, 10, 0)
    return x_relic_line, y_relic_line
end

function BBN_line(x_data, y_data, BBN_values, delta_values, delta, Yukawa_values, Yukawa) 
    BBN_inds = excluded_BBN_inds(BBN_values, y_data, delta_values, delta, Yukawa_values, Yukawa)
    x_BBN = x_data[BBN_inds]
    y_BBN = y_data[BBN_inds]
    x_BBN_line, y_BBN_line = line_finder(x_BBN, y_BBN, 10, 1)
    return x_BBN_line, y_BBN_line
end

#Define the limits according to physics
const Relic_abundance_limit = 0.12 
const BBN_lifetime = 1.52*1E24 #Upper bound on glueball decay rate in Gev^-1.

test_df = DataFrame(CSV.File("NplusL_model_scan.csv")) #Read in data as a dataframe

#Initialise the data in respective arrays
y_values = test_df[!, 2] #values for Lambda
x_values = test_df[!, 1] ./ test_df[!, 2] #values for m/Lambda
delta_values = test_df[!, 5] #values for delta
Yukawa_values = test_df[!, 7] #values for yd
GB_rate_values = test_df[!, 13]
relic_abundance_values = test_df[!, 15] 

#BLOCK for colour maps w/ fixed delta

Delta_fixed_value = 4.6
Yukawa1 = 1E-1
Yukawa2 = 1E-2
Yukawa3 = 1E-3

x_relic1,  y_relic1 = relic_line(x_values, y_values, relic_abundance_values, delta_values, Delta_fixed_value, Yukawa_values, Yukawa1)
x_BBN1, y_BBN1 = BBN_line(x_values, y_values, GB_rate_values, delta_values, Delta_fixed_value, Yukawa_values, Yukawa1) 
x_relic2,  y_relic2 = relic_line(x_values, y_values, relic_abundance_values, delta_values, Delta_fixed_value, Yukawa_values, Yukawa2)
x_BBN2, y_BBN2 = BBN_line(x_values, y_values, GB_rate_values, delta_values, Delta_fixed_value, Yukawa_values, Yukawa2) 
x_relic3,  y_relic3 = relic_line(x_values, y_values, relic_abundance_values, delta_values, Delta_fixed_value, Yukawa_values, Yukawa3)
x_BBN3, y_BBN3 = BBN_line(x_values, y_values, GB_rate_values, delta_values, Delta_fixed_value, Yukawa_values, Yukawa3) 

x_mass1_line = 10.0.^collect(range(2, 4, length = 100))
y_mass1_line = (1*10^0)./x_mass1_line

x_mass2_line = 10.0.^collect(range(2, 4, length = 100))
y_mass2_line = (1*10^3)./x_mass2_line

x_mass3_line = 10.0.^collect(range(2, 4, length = 100))
y_mass3_line = (1*10^6)./x_mass3_line

x_mass4_line = 10.0.^collect(range(2, 4, length = 100))
y_mass4_line = (1*10^9)./x_mass4_line

plot(x_relic1, y_relic1, xlims = (1E2, 1E4), ylims = (1, 1E7), title="N + L model: "*L"m_L = "*"$Delta_fixed_value "*L"m_N",
 minorticks = 10,  minorgrid = true, ylabel = L"\Lambda/GeV",
  xlabel = L"m_N/\Lambda", xscale = :log10, yscale = :log10, labels=L"y_d = 0.1",legend=false,
   seriescolor=:blue, fillrange = 1E7*ones(size(x_relic1)), fillalpha = 0.35, c = 1,)
plot!(x_BBN1, y_BBN1, seriescolor=:red, labels=L"y_d = 0.1", linestyle=:dash)
plot!(x_relic2, y_relic2, seriescolor=:blue, labels=L"y_d = 0.01", linestyle=:dash)
plot!(x_BBN2, y_BBN2, seriescolor=:red, labels=L"y_d = 0.01", linestyle=:dash)
plot!(x_relic3, y_relic3, seriescolor=:blue, labels=L"y_d = 0.001", linestyle=:dash)
plot!(x_BBN3, y_BBN3, seriescolor=:red, labels=L"y_d = 0.001", fillrange = ones(size(x_BBN3)), fillalpha = 0.35, c = 2)
plot!(x_mass1_line, y_mass1_line, c=:black, linestyle=:dash)
plot!(x_mass2_line, y_mass2_line, c=:black, linestyle=:dash)
plot!(x_mass3_line, y_mass3_line, c=:black, linestyle=:dash)
plot!(x_mass4_line, y_mass4_line, c=:black, linestyle=:dash)
annotate!(500, 0.01, text(L"m_Q = 1"*" GeV", :black, :right, 12))
annotate!(300, 15, text(L"m_Q = 1"*" TeV", :black, :right, 12))
annotate!(800, 850, text(L"m_Q = 10^3"*" TeV", :black, :right, 12))
annotate!(2000, 3*10^6, text(L"m_Q = 10^6"*" TeV", :black, :right, 12))
annotate!(7000, 15, text(L"\tau_{GB} > 1s", :black, :right, 12))
annotate!(900, 2E5, text(L"\Omega_{DM} h^2>0.12", :black, :right, 12))
annotate!(200, 3E6, text(L"y_d = 0.001", :blue, :right, 9))
annotate!(200, 1E5, text(L"y_d = 0.1", :blue, :right, 9))
#annotate!(2000, 5, text(L"y_d = 0.1", :red, :right, 9))
annotate!(2500, 3000, text(L"y_d = 0.001", :red, :right, 9))
#annotate!(4000, 800000, text(L"y_d = 0.01", :red, :right, 9))
savefig("Exclusion_Plot_NplusL_variable_Yukawa.png")

#BLOCK for colour maps w/ fixed yd
#=
Yukawa_coupling = 1E-3
delta1 = 1
delta2 = 4.6
delta3 = 35.9

x_relic1,  y_relic1 = relic_line(x_values, y_values, relic_abundance_values, delta_values, delta1, Yukawa_values, Yukawa_coupling)
x_BBN1, y_BBN1 = BBN_line(x_values, y_values, GB_rate_values, delta_values, delta1, Yukawa_values, Yukawa_coupling) 
x_relic2,  y_relic2 = relic_line(x_values, y_values, relic_abundance_values, delta_values, delta2, Yukawa_values, Yukawa_coupling)
x_BBN2, y_BBN2 = BBN_line(x_values, y_values, GB_rate_values, delta_values, delta2, Yukawa_values, Yukawa_coupling) 
x_relic3,  y_relic3 = relic_line(x_values, y_values, relic_abundance_values, delta_values, delta3, Yukawa_values, Yukawa_coupling)
x_BBN3, y_BBN3 = BBN_line(x_values, y_values, GB_rate_values, delta_values, delta3, Yukawa_values, Yukawa_coupling) 

x_mass1_line = 10.0.^collect(range(2, 4, length = 100))
y_mass1_line = (1*10^0)./x_mass1_line

x_mass2_line = 10.0.^collect(range(2, 4, length = 100))
y_mass2_line = (1*10^3)./x_mass2_line

x_mass3_line = 10.0.^collect(range(2, 4, length = 100))
y_mass3_line = (1*10^6)./x_mass3_line

x_mass4_line = 10.0.^collect(range(2, 4, length = 100))
y_mass4_line = (1*10^9)./x_mass4_line

plot(x_relic1, y_relic1, xlims = (1E2, 1E4), ylims = (1, 1E7), title="N + L model: yd = $Yukawa_coupling",
 minorticks = 10,  minorgrid = true, ylabel = L"\Lambda/GeV",
  xlabel = L"m_N/\Lambda", xscale = :log10, yscale = :log10, labels=L"m_L = m_N",legend=false,
   seriescolor=:blue, fillrange = 1E7*ones(size(x_relic1)), fillalpha = 0.35, c = 1,)
plot!(x_BBN1, y_BBN1, seriescolor=:red, labels=L"m_L =m_N", linestyle=:dash)
plot!(x_relic2, y_relic2, seriescolor=:blue, labels=L"m_L = 5*m_N", linestyle=:dash)
plot!(x_BBN2, y_BBN2, seriescolor=:red, labels=L"m_L = 5*m_N", linestyle=:dash)
plot!(x_relic3, y_relic3, seriescolor=:blue, labels=L"m_L = 36*m_N", linestyle=:dash)
plot!(x_BBN3, y_BBN3, seriescolor=:red, labels=L"m_L = 36*m_N", fillrange = ones(size(x_BBN3)), fillalpha = 0.35, c = 2)
plot!(x_mass1_line, y_mass1_line, c=:black, linestyle=:dash)
plot!(x_mass2_line, y_mass2_line, c=:black, linestyle=:dash)
plot!(x_mass3_line, y_mass3_line, c=:black, linestyle=:dash)
plot!(x_mass4_line, y_mass4_line, c=:black, linestyle=:dash)
annotate!(500, 0.01, text(L"m_Q = 1"*" GeV", :black, :right, 12))
annotate!(300, 15, text(L"m_Q = 1"*" TeV", :black, :right, 12))
annotate!(800, 10000, text(L"m_Q = 10^3"*" TeV", :black, :right, 12))
annotate!(2000, 3*10^6, text(L"m_Q = 10^6"*" TeV", :black, :right, 12))
annotate!(7000, 15, text(L"\tau_{GB} > 1s", :black, :right, 12))
annotate!(500, 8E6, text(L"\Omega_{DM} h^2>0.12", :black, :right, 9))
annotate!(200, 1E5, text(L"m_L = m_N", :blue, :right, 9))
annotate!(200, 2E6, text(L"m_L = 5 m_N", :blue, :right, 9))
annotate!(2000, 5, text(L"m_L = m_N", :red, :right, 9))
annotate!(2000, 1000, text(L"m_L = 5 m_N", :red, :right, 9))
annotate!(4000, 700000, text(L"m_L = 36 m_N", :red, :right, 9))
savefig("Exclusion_Plot_NplusL_variable_delta.png")
=#

#= #BLOCK for the simple plots
Yukawa_coupling = 1E-1 #For simple plots
Mass_delta = 12.9154966

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
