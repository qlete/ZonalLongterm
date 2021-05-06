function investment_fbmc_centralized()
    # Build capacity investment model
    m = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))

    # define primal variables
    @variable(m, y[1:length(I), 1:size(D, 1), 1:length(Z)] >= 0) # power produced by tech i in slice j
    @variable(m, yrv[1:length(I), 1:size(D, 1), 1:length(Z)] >= 0) # power reserved by tech i in slice j
    @variable(m, yre[1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, x[1:length(I), 1:length(Z)] >= 0) # investment in technology i
    @variable(m, s[j=1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, y_bar[1:length(X_bar), 1:size(D, 1)] >= 0)
    @variable(m, yrv_bar[1:length(X_bar), 1:size(D, 1)] >= 0)
    @variable(m, x_bar[1:length(X_bar)] >= 0)
    # define auxiliary variables
    @variable(m, y_tilde[1:size(D, 1), 1:length(N)] >= 0) # power produced by tech i in slice j
    @variable(m, x_tilde[1:length(I), 1:length(N)] >= 0) # investment in technology i
    @variable(m, -TC[l] <= f_tilde[j=1:size(D, 1), l=1:length(L)] <= TC[l]) # export from node n
    @variable(m, Z_tilde[1:length(N)] >= 0)
    @variable(m, z_tilde[1:size(D, 1), 1:length(N)] >= 0)
    @variable(m, p[1:size(D, 1), 1:length(Z)])
    @variable(m, r[1:size(D, 1), 1:length(N)])

    # define primal constraints
    @constraint(m, rho[j=1:size(D, 1), z=1:length(Z)], p[j,z] == sum(y[i,j,z]
        for i=1:length(I)) + yre[j,z] + sum(y_bar[g,j] for g=Gz[z]) -
        sum(D[j,n] for n in Nz[z]) + s[j,z])
    @constraint(m, rhoz_tilde[j=1:size(D, 1), z=1:length(Z)], p[j,z] ==
        sum(y_tilde[j,n] + z_tilde[j,n] for n in Nz[z]) - sum(D[j,n] for n in Nz[z]))
    @constraint(m, rhon_tilde[j=1:size(D, 1), n=1:length(N)], r[j,n] == z_tilde[j,n] + y_tilde[j,n] - D[j,n])
    @constraint(m, [j=1:size(D, 1)], sum(r[j,n] for n in 1:length(N)) == 0)
    @constraint(m, [j=1:size(D, 1), l=1:length(L)], f_tilde[j,l] == sum(PTDF[n,l]*r[j,n] for n=1:length(N)))
    @constraint(m, mu[i=1:length(I), j=1:size(D, 1), z=1:length(Z)], y[i,j,z] + yrv[i,j,z] <= x[i,z])
    @constraint(m, [j=1:size(D, 1), z=1:length(Z)], yre[j,z] <= sum(R[j,n] for n in Nz[z]))
    @constraint(m, [g=1:length(X_bar)], x_bar[g] <= X_bar[g])
    @constraint(m, [i=1:length(I)], sum(x[i,z] for z in 1:length(Z)) <= X[i])
    @constraint(m, [g=1:length(X_bar), j=1:size(D, 1)],
        y_bar[g,j] + yrv_bar[g,j] <= x_bar[g])
    # define auxiliary constraints
    @constraint(m, [j=1:size(D, 1), n=1:length(N)],
        y_tilde[j,n] <= sum(x_tilde[i,n] for i in 1:length(I)) +
        sum(x_bar[g] for g=Gn[n]) + R[j,n])
    @constraint(m, [g=1:length(X_bar)], x_bar[g] <= X_bar[g])
    @constraint(m, [g=1:length(X_bar), j=1:size(D, 1)],
        y_bar[g,j] <= x_bar[g])
    @constraint(m, [j=1:size(D, 1), n=1:length(N)], z_tilde[j,n] <= Z_tilde[n])
    @constraint(m, [i=1:length(I), z=1:length(Z)], sum(x_tilde[i,n] for n in Nz[z]) <= x[i,z])
    @constraint(m, rho_rv[z=1:length(Z), j=1:size(D, 1)], sum(yrv[i,j,z] for i=1:2) + 
        sum(yrv_bar[g,j] for g in Gz[z] if isfastG[g]) >= Rv[z])

    @objective(m, Min, sum(DT[j]*(VOLL*s[j,z] +
        sum(MC[i]*(y[i,j,z] + yrv[i,j,z]) for i=1:length(I))) for j=1:size(D, 1), z=1:length(Z)) +
        sum((I[i] + FC[i])*x[i,z] for i=1:length(I), z=1:length(Z)) +
        sum(DT[j]*MC_tilde*z_tilde[j,n] for j=1:size(D, 1), n=1:length(N)) +
        sum(DT[j]*MC_ex[g]*(y_bar[g,j] + yrv_bar[g,j]) for g=1:length(X_bar), j=1:size(D, 1)) + 
        sum(FC_ex[g]*x_bar[g] for g=1:length(X_bar)) +
        I_tilde*sum(Z_tilde))

    optimize!(m)
    
    sol = Dict()
    push!(sol, "p" => value.(p))
    push!(sol, "r" => value.(r))
    push!(sol, "y" => value.(y))
    push!(sol, "yrv" => value.(yrv))
    push!(sol, "yre" => value.(yre))
    push!(sol, "x" => value.(x))
    push!(sol, "s" => value.(s))
    push!(sol, "y_bar" => value.(y_bar))
    push!(sol, "yrv_bar" => value.(yrv_bar))
    push!(sol, "x_bar" => value.(x_bar))
    push!(sol, "y_tilde" => value.(y_tilde))
    push!(sol, "x_tilde" => value.(x_tilde))
    push!(sol, "f_tilde" => value.(f_tilde))
    push!(sol, "Z_tilde" => value.(Z_tilde))
    push!(sol, "z_tilde" => value.(z_tilde))
    push!(sol, "rho" => -dual.(rho))

    return sol
end

function redispatch_fbmc(x_ep::Array{Float64, 2}, x_bar_ep::Array{Float64, 1},
    Z_tj::Array{Float64, 1})
    # determine turbojet investment necessary and redispatch
    model_rd = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))

    # declare variables
    @variable(model_rd, y_rd[1:length(I), 1:size(D, 1), 1:length(N)] >= 0)
    @variable(model_rd, yre_rd[1:size(D, 1), 1:length(N)] >= 0)
    @variable(model_rd, y_bar_rd[1:length(X_bar), 1:size(D, 1)] >= 0)
    @variable(model_rd, s_rd[1:size(D, 1), 1:length(N)] >= 0)
    @variable(model_rd, x_rd[1:length(I), 1:length(N)] >= 0) # zonal investment disaggregation
    @variable(model_rd, r_rd[1:size(D, 1), 1:length(N)])
    @variable(model_rd, -TC[l] <= f_rd[j=1:size(D, 1), l=1:length(L)] <= TC[l])
    @variable(model_rd, z_tilde_rd[1:size(D, 1), 1:length(N)] >= 0)

    # declare constraints

    @constraint(model_rd, [i=1:length(I), z=1:length(Z)], sum(x_rd[i,n] for n in Nz[z]) == x_ep[i,z])
    @constraint(model_rd, [j=1:size(D, 1), n=1:length(N)], z_tilde_rd[j,n] <= Z_tj[n])
    @constraint(model_rd, [i=1:length(I), j=1:size(D, 1), n=1:length(N)], y_rd[i,j,n] <= x_rd[i,n])
    @constraint(model_rd, [j=1:size(D, 1), n=1:length(N)], yre_rd[j,n] <= R[j,n])
    @constraint(model_rd, [g=1:length(X_bar), j=1:size(D, 1)], y_bar_rd[g,j] <= x_bar_ep[g])
    @constraint(model_rd, rho[j=1:size(D, 1), n=1:length(N)], r_rd[j,n] ==  z_tilde_rd[j,n] +
        yre_rd[j,n] + sum(y_rd[i,j,n] for i=1:length(I)) + sum(y_bar_rd[g,j] for g=Gn[n]) +
        s_rd[j,n] - D[j,n])

    @constraint(model_rd, [j=1:size(D, 1)], sum(r_rd[j,n] for n=1:length(N)) == 0)
    @constraint(model_rd, [j=1:size(D, 1), l=1:length(L)], f_rd[j,l] ==
        sum(PTDF[n,l]*r_rd[j,n] for n=1:length(N)))

    # declare objective
    @objective(model_rd, Min, sum(DT[j]*(MC_tilde*z_tilde_rd[j,n] + VOLL*s_rd[j,n] +
        sum(MC[i]*y_rd[i,j,n] for i=1:length(I))) for j in 1:size(D, 1), n in 1:length(N)) +
        sum(DT[j]*MC_ex[g]*y_bar_rd[g,j] for g=1:length(X_bar), j in 1:size(D, 1)))

    optimize!(model_rd)
    
    sol = Dict()
    push!(sol, "r" => value.(r_rd))
    push!(sol, "y" => value.(y_rd))
    push!(sol, "yre" => value.(yre_rd))
    push!(sol, "x" => value.(x_rd))
    push!(sol, "s" => value.(s_rd))
    push!(sol, "y_bar" => value.(y_bar_rd))
    push!(sol, "f" => value.(f_rd))
    push!(sol, "z_tilde" => value.(z_tilde_rd))

    return sol
end