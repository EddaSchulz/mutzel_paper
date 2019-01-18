### Simulate polyploid cells with XA repression hypothesis in stochastic Gillespie model
include("sim_functions_gil_Oct2018.jl")

println(ARGS[1])
println(ARGS[2])

pars=convert(Array{Float64},readcsv(ARGS[1]))

# Simulate all parameters from stochastic simulation for initial cond XaXa
nr_cells=100
#Simulation time
t = 100.0

ma_scan=Array{Float64}(size(pars,1)*nr_cells,33+4*(Int64(t)+1))
@time for n=1:size(pars,1)
  println(n)
  p=pars[n,:]
  println(p)
  ma_scan[(n-1)*nr_cells+1:n*nr_cells,1:33] = p'.*ones(nr_cells,33)
  ma_scan[(n-1)*nr_cells+1:n*nr_cells,34:33+4*(Int64(t)+1)] = sim_ba_gil_polypl_repXA(p,nr_cells,t)
end

writecsv(ARGS[2],ma_scan)
