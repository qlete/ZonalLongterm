using Pkg; Pkg.activate("..");
## load external modules
using CSV, DataFrames, JuMP, Gurobi, Plots, LinearAlgebra, Plots.PlotMeasures, PyPlot
using LightGraphs, Compose, Colors, GraphPlot, Gadfly, PlotlyJS
const GRB_ENV = Gurobi.Env()
const TOL = 1e-2
const plt = PyPlot

## load internal modules
push!(LOAD_PATH, "modules/")
using ReduceNetwork

## include internal files
#data
include("../data/data_3nodes.jl")
include("../data/data_cwe.jl")
#models
include("models/expansion_nodal.jl")
include("models/helpers_fbmc.jl")
include("models/expansion_fbmc_centralized.jl")
include("models/splitting_decentralized.jl")
include("models/expansion_zonal_pa.jl")
include("plot_investment.jl")

## build reduced network

write_reduced_network("../data/CWE2018_daily/", "../data/reduced_net/", 100, 20)

## CWE results
#### load data
G_Nidx, Gn, G_Zidx, Gz, isfastG, DT, MC, I, FC, MC_ex, FC_ex, MC_tilde, I_tilde, VOLL, D, 
    R, Rv, N, Z, Nz, Zn, X, X_bar, L, TC, PTDF = load_data_cwe()
# nodal 
sol_nodal = investment_nodal()
# centralized FBMC
sol_fci = investment_fbmc_centralized()
sol_fcr = redispatch_fbmc(sol_fci["x"], sol_fci["x_bar"], sol_fci["Z_tilde"])
# decentralized FBMC
sol_fdi = investment_fbmc_decentralized()
sol_fdr = redispatch_fbmc(sol_fdi["x"], sol_fdi["x_bar"], sol_fdi["Z_tilde"])
# zonal PA
sol_pai = investment_pa()
sol_par = redispatch_pa(sol_pai["x"], sol_pai["x_bar"])

#### make plots
plot_investment(sol_nodal, sol_fci, sol_fcr, sol_fdi, sol_fdr, sol_pai, sol_par)
plot_net(sol_nodal, sol_fci, sol_fdi, sol_pai)
plot_welfare(sol_nodal, sol_fci, sol_fcr, sol_fdi, sol_fdr, sol_pai, sol_par)
