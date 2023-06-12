include("Triple_Solver.jl")

#Parameters
m_quark = 100000
Lambda_dQCD = 100

#Define parameters of implicit Euler backward solution method
x_initial = 5
x_final = 1E7
Delta_t = 1E-1
tvec = collect(log(x_initial):Delta_t:log(x_final))
xvec = exp.(tvec)
Npoints = length(xvec)

# Define which processes should be included in the analysis - THE ORDER CRUCIALLY MATTERS AND MUST NOT BE INTERCHANGED!
qqTogg = true #No.1: Quark antiquark annihilation into gluons
DDTogg = true #No.2: Diquark rearrangement annihilation into gluons 
BBTogg = true #No.3: Baryon rearrangement annihilation into gluons 
qqToDg = true #No.4: Quark capture into diquarks
qDToBg = true #No.5: Diquark capture into baryons
qBToqq = false #No.6: Rearrangement quark + baryon into 2 quarks
qBToDg = false #No.7: Rearrangement quark + baryon into diquark + gluons
DqToqg = false #No.8: Rearrangement diquark + quark into quark and gluons
DBToqg = false #No.9: Rearrangement baryon + diquark into guark and gluon 
DDToBq = false #No.10: Rearrangement 2 diquarks into baryon and quark 
BDToDq = false #No.11: Rarrangement diquark + baryon into diquark and quark 
BBToDD = false #No.12: Rearrangement 2 baryons -> 2 diquarks 
DDToqq = false #No.13: Rearrangement 2 diquarks -> 2 quarks 
BBToqq = false #No.14: Rearrangement 2 baryons -> 2 quarks 

active_processes = [qqTogg, DDTogg, BBTogg, qqToDg, qDToBg, qBToqq, qBToDg, DqToqg, DBToqg, DDToBq, BDToDq, BBToDD, DDToqq, BBToqq]

YxQuark, YxDiQuark, YxBaryon = Coupled_Freeze_Out(x_initial, x_final, Delta_t, m_quark, Lambda_dQCD, active_processes)

results_file = open("FO_data.csv", "w")
IOStream("FO_data.csv")
write(results_file, "x, Y1, Y2, Y3 \n")
for i in 1:Npoints
    write(results_file, join((xvec[i],YxQuark[i], YxDiQuark[i], YxBaryon[i]), ","),"\n")
end
close(results_file)

plot(xvec, YxQuark, xscale = :log10, yscale = :log10, xlabel = "x = m/T", ylabel = "Y(x)", label=L"Y_1(x)",title = "Freezeout", minorgrid = true, minorticks = 10)
plot!(xvec, YxDiQuark, label=L"Y_2(x)")
plot!(xvec, YxBaryon, label=L"Y_3(x)")
savefig("Test.png")
