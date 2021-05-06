function investment_nodal()

    # Build capacity investment model
    m = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))

    # define primal variables
    @variable(m, y[1:length(I), 1:size(D, 1), 1:length(N)] >= 0) # power produced by tech i in slice j
    @variable(m, yrv[1:length(I), 1:size(D, 1), 1:length(N)] >= 0) # power produced by tech i in slice j
    @variable(m, yre[1:size(D, 1), 1:length(N)] >= 0) #renewable power produced in period j, node n
    @variable(m, x[1:length(I), 1:length(N)] >= 0) # investment in technology i
    @variable(m, s[j=1:size(D, 1), 1:length(N)] >= 0)
    @variable(m, y_bar[1:length(X_bar), 1:size(D, 1)] >= 0)
    @variable(m, yrv_bar[1:length(X_bar), 1:size(D, 1)] >= 0)
    @variable(m, x_bar[1:length(X_bar)] >= 0)
    # define auxiliary variables
    @variable(m, -TC[l] <= f[j=1:size(D, 1), l=1:length(L)] <= TC[l]) # export from node n
    @variable(m, r[1:size(D, 1), 1:length(N)])

    # define primal constraints
    @constraint(m, rho[j=1:size(D, 1), n=1:length(N)], r[j,n] == yre[j,n] + sum(y[i,j,n] for i=1:length(I)) +
        sum(y_bar[g,j] for g=Gn[n]) - D[j,n] + s[j,n])
    @constraint(m, [j=1:size(D, 1)], sum(r[j,n] for n in 1:length(N)) == 0)
    @constraint(m, [j=1:size(D, 1), l=1:length(L)], f[j,l] == sum(PTDF[n,l]*r[j,n] for n=1:length(N)))
    @constraint(m, [i=1:length(I), j=1:size(D, 1), n=1:length(N)], y[i,j,n] + yrv[i,j,n] <= x[i,n])
    @constraint(m, [i=1:length(I), j=1:size(D, 1), n=1:length(N)], yre[j,n] <= R[j,n])
    @constraint(m, [g=1:length(X_bar)], x_bar[g] <= X_bar[g])
    @constraint(m, [i=1:length(I)], sum(x[i,n] for n in 1:length(N)) <= X[i])
    @constraint(m, [g=1:length(X_bar), j=1:size(D, 1)], y_bar[g,j] + yrv_bar[g,j] <= x_bar[g])
    @constraint(m, reserve[z=1:length(Z), j=1:size(D, 1)], sum(yrv[i,j,n] for i=1:2, n=Nz[z]) + 
        sum(yrv_bar[g,j] for g in Gz[z] if isfastG[g]) >= Rv[z])

    @objective(m, Min, sum(DT[j]*(VOLL*s[j,n] +
        sum(MC[i]*(y[i,j,n] + yrv[i,j,n]) for i=1:length(I))) for j=1:size(D, 1), n=1:length(N)) +
        sum((I[i] + FC[i])*x[i,n] for i=1:length(I), n=1:length(N)) +
        sum(DT[j]*MC_ex[g]*(y_bar[g,j] + yrv_bar[g,j]) for g=1:length(X_bar), j=1:size(D, 1)) +
        sum(FC_ex[g]*x_bar[g] for g=1:length(X_bar)))

    optimize!(m)
    
    sol = Dict()
    push!(sol, "r" => value.(r))
    push!(sol, "y" => value.(y))
    push!(sol, "yrv" => value.(yrv))
    push!(sol, "yre" => value.(yre))
    push!(sol, "x" => value.(x))
    push!(sol, "x_bar" => value.(x_bar))
    push!(sol, "s" => value.(s))
    push!(sol, "y_bar" => value.(y_bar))
    push!(sol, "yrv_bar" => value.(yrv_bar))
    push!(sol, "x_bar" => value.(x_bar))
    push!(sol, "f" => value.(f))
    push!(sol, "rho" => -dual.(rho))

    return sol
end