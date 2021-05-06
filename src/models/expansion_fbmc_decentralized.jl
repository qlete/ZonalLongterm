using JuMP, Gurobi, LinearAlgebra, Plots, Plots.PlotMeasures, BSON

# include("../data/data_3nodes.jl")
include("../data/data_cwe.jl")
include("helpers_fbmc.jl")
include("expansion_fbmc_centralized.jl")
include("splitting_decentralized.jl")
const TOL = 1e-2
const GRB_ENV = Gurobi.Env()

G_Nidx, Gn, G_Zidx, Gz, isfastG, DT, MC, I, FC, MC_ex, FC_ex, MC_tilde, I_tilde, VOLL, D, 
    R, Rv, N, Z, Nz, Zn, X, X_bar, L, TC, PTDF, eta = load_data(true)

sol_inv = investment_fbmc_decentralized()
sol_rd = redispatch_fbmc(sol_inv["x"], sol_inv["x_bar"], sol_inv["Z_tilde"])

# save solutions to files
bson("./results/solutions/fbmc_decentralized_inv.bson", sol_inv)
bson("./results/solutions/fbmc_decentralized_rd.bson", sol_rd)

# prices = (sol["rhop"] .- sol["rhom"]) ./ DT
# 
# investment_cost = sum((I[i]+FC[i])*x_ep[i,z] for i=1:length(I), z=1:length(Z)) +
#     sum(FC_ex[g]*x_bar_ep[g] for g=1:length(X_bar)) +
#     I_tilde*sum(Z_tj[n] for n in 1:length(N))
# @show investment_cost
# total_cost = op_cost + investment_cost
# @show total_cost
# @show sum(X_bar .- x_bar_ep)
# 
# weighted_price = sum(sum(D[j,n] for n in Nz[z])*prices[j,z] for j in 1:size(D, 1), z in 1:length(Z))/sum(D[j,n] for j in 1:size(D, 1), n in 1:length(N))
# @show 1000*weighted_price
# 
# CSV.write("../results/FBMC_decent/prices.csv", Tables.table(prices))
