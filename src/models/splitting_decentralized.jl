# this file is dedicated to an implementation of the linear splitting method
# with norm regularization

function investment_fbmc_decentralized(;sol::Union{Nothing, Dict}=nothing)
    println("#### DECENTRALIZED FBMC CAPACITY EXPANSION ####")
    # set the number of maximum iterations
    nb_it_max = 200
    
    # load the initial solution
    # if not given, start with a zero solution
    if isnothing(sol)
        status, sol = starting_point_zero()
    else
        status = MOI.OPTIMAL
    end
    
    # initialize parameters
    k = 1
    issol = false
    dist = 1.0
    oldsol = sol
    
    println("> Start the splitting algorithm")
    println("Iteration \t sol \t Distance \t Inv. \t NR \t Cap payment")
    while !isapprox(dist, 0, atol=TOL) && !issol
        try
            status, dist, sol = primal_dual_augmented_oracle(sol)
        catch e
            println(">>> Error in the primal-dual augmented oracle. Use the backup problem.")
            status, _, sol = augmented_oracle(sol)
            dist = distance(oldsol, sol)
        end
        issol = check_solution(sol)
        
        println(k, "\t", issol, "\t", round(dist, digits=2), "\t", 
            round(sum(sol["x"]), digits=2), "\t", round(sum(sol["Z_tilde"]), digits=2), 
            "\t", round(sum(sol["nu_tilde"]), digits=2))
        
        k += 1
        if k > nb_it_max
            break
        end
    end

    return sol
end

# centralized problem augmented by the fixed point term
function augmented_oracle(sol::Dict)
    # Build capacity investment model
    m = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attributes(m, "OutputFlag" => 0, "Threads" => 4)

    # define primal variables
    @variable(m, y[1:length(I), 1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, yrv[1:length(I), 1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, yre[1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, x[1:length(I), 1:length(Z)] >= 0) # investment in technology i
    @variable(m, s[j=1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, y_bar[1:length(X_bar), 1:size(D, 1)] >= 0)
    @variable(m, yrv_bar[1:length(X_bar), 1:size(D, 1)] >= 0)
    @variable(m, x_bar[1:length(X_bar)] >= 0)
    # define auxiliary variables
    @variable(m, y_tilde[1:size(D, 1), 1:length(N)] >= 0) # power produced by tech i in slice j
    @variable(m, x_tilde[1:length(I), 1:length(N)] >= 0) # investment in technology i
    @variable(m, f_tildep[j=1:size(D, 1), l=1:length(L)] >= 0)
    @variable(m, f_tildem[j=1:size(D, 1), l=1:length(L)] >= 0)
    @variable(m, pp[1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, pm[1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, rp[1:size(D, 1), 1:length(N)] >= 0)
    @variable(m, rm[1:size(D, 1), 1:length(N)] >= 0)
    @variable(m, Z_tilde[n=1:length(N)] >= 0)
    @variable(m, z_tilde[j=1:size(D, 1), n=1:length(N)] >= 0)

    # define primal constraints
    @constraint(m, mu_bar[g=1:length(X_bar), j=1:size(D, 1)], x_bar[g] - y_bar[g,j] - yrv_bar[g,j] >= 0)
    @constraint(m, delta_bar[g=1:length(X_bar)], X_bar[g] - x_bar[g] >= 0)
    @constraint(m, delta[i=1:length(I)], X[i] - sum(x[i,z] for z in 1:length(Z)) >= 0)
    @constraint(m, mure[j=1:size(D, 1), z=1:length(Z)], sum(R[j,n] for n in Nz[z]) - yre[j,z] >= 0)
    @constraint(m, rhop[j=1:size(D, 1), z=1:length(Z)], -(pp[j,z] - pm[j,z]) + yre[j,z] + sum(y[i,j,z] for i=1:length(I)) + sum(y_bar[g,j] for g=Gz[z]) - sum(D[j,n] for n in Nz[z]) + s[j,z] >= 0)
    @constraint(m, rhom[j=1:size(D, 1), z=1:length(Z)], (pp[j,z] - pm[j,z]) - yre[j,z] - sum(y[i,j,z] for i=1:length(I)) - sum(y_bar[g,j] for g=Gz[z]) + sum(D[j,n] for n in Nz[z]) - s[j,z] >= 0)
    @constraint(m, rhoz_tildep[j=1:size(D, 1), z=1:length(Z)], -(pp[j,z] - pm[j,z]) + sum(z_tilde[j,n] for n in Nz[z]) + sum(y_tilde[j,n] for n in Nz[z]) - sum(D[j,n] for n in Nz[z]) >= 0)
    @constraint(m, rhoz_tildem[j=1:size(D, 1), z=1:length(Z)], (pp[j,z] - pm[j,z]) - sum(z_tilde[j,n] for n in Nz[z]) - sum(y_tilde[j,n] for n in Nz[z]) + sum(D[j,n] for n in Nz[z]) >= 0)
    @constraint(m, rhon_tildep[j=1:size(D, 1), n=1:length(N)], (rp[j,n] - rm[j,n]) - z_tilde[j,n] - y_tilde[j,n] + D[j,n] >= 0)
    @constraint(m, rhon_tildem[j=1:size(D, 1), n=1:length(N)], -(rp[j,n] - rm[j,n]) + z_tilde[j,n] + y_tilde[j,n] - D[j,n] >= 0)
    @constraint(m, rho_rv[j=1:size(D,1), z=1:length(Z)], -Rv[z] + sum(yrv[i,j,z] for i=1:2) + sum(yrv_bar[g,j] for g in Gz[z] if isfastG[g]) >= 0)
    @constraint(m, phip[j=1:size(D, 1)], sum((rp[j,n] - rm[j,n]) for n in 1:length(N)) >= 0)
    @constraint(m, phim[j=1:size(D, 1)], -sum((rp[j,n] - rm[j,n]) for n in 1:length(N)) >= 0)
    @constraint(m, psip[j=1:size(D, 1), l=1:length(L)], (f_tildep[j,l] - f_tildem[j,l]) - sum(PTDF[n,l]*(rp[j,n] - rm[j,n]) for n=1:length(N)) >= 0)
    @constraint(m, psim[j=1:size(D, 1), l=1:length(L)], -(f_tildep[j,l] - f_tildem[j,l]) + sum(PTDF[n,l]*(rp[j,n] - rm[j,n]) for n=1:length(N)) >= 0)
    @constraint(m, mu[i=1:length(I), j=1:size(D, 1), z=1:length(Z)], x[i,z] - y[i,j,z] - yrv[i,j,z] >= 0)
    # define auxiliary constraints
    @constraint(m, mu_tilde[j=1:size(D, 1), n=1:length(N)], R[j,n] + sum(x_bar[g] for g in Gn[n]) + sum(x_tilde[i,n] for i=1:length(I)) - y_tilde[j,n] >= 0)
    @constraint(m, nu_tilde[z=1:length(Z)], sum(x[i,z] for i=1:length(I)) - sum(x_tilde[i,n] for i=1:length(I), n in Nz[z]) >= 0)
    @constraint(m, lambdap[j=1:size(D, 1), l=1:length(L)], TC[l] - (f_tildep[j,l] - f_tildem[j,l]) >= 0)
    @constraint(m, lambdam[j=1:size(D, 1), l=1:length(L)], TC[l] + (f_tildep[j,l] - f_tildem[j,l]) >= 0)
    @constraint(m, gamma_tilde[j=1:size(D, 1), n=1:length(N)], Z_tilde[n] - z_tilde[j,n] >= 0)

    @expression(m, original_obj, sum(DT[j]*(VOLL*s[j,z] + sum(MC[i]*(y[i,j,z] + yrv[i,j,z])
        for i=1:length(I))) for j=1:size(D, 1), z=1:length(Z)) +
        sum(DT[j]*MC_ex[g]*(y_bar[g,j] + yrv_bar[g,j]) for g=1:length(X_bar), j=1:size(D, 1)) +
        sum((I[i] + FC[i])*x[i,z] for i=1:length(I), z=1:length(Z)) +
        sum(FC_ex[g]*x_bar[g] for g=1:length(X_bar)) +
        sum(DT[j]*MC_tilde*z_tilde[j,n] for j=1:size(D, 1), n=1:length(N)) +
        I_tilde*sum(Z_tilde[n] for n=1:length(N)))

    @expression(m, augmented_term, sum(x[i,z]*sol["nu_tilde"][z]
        for i=1:length(I), z=1:length(Z)) +
        sum(x_bar[g]*(sum(sol["mu_tilde"][j,G_Nidx[g]] for j=1:size(D, 1)))
        for g=1:length(X_bar)))

    @objective(m, Min, original_obj + augmented_term)

    optimize!(m)

    # return solution
    sol_final = Dict()
    push!(sol_final, "pp" => value.(pp))
    push!(sol_final, "pm" => value.(pm))
    push!(sol_final, "rp" => value.(rp))
    push!(sol_final, "rm" => value.(rm))
    push!(sol_final, "y" => value.(y))
    push!(sol_final, "yrv" => value.(yrv))
    push!(sol_final, "yre" => value.(yre))
    push!(sol_final, "x" => value.(x))
    push!(sol_final, "s" => value.(s))
    push!(sol_final, "y_bar" => value.(y_bar))
    push!(sol_final, "yrv_bar" => value.(yrv_bar))
    push!(sol_final, "x_bar" => value.(x_bar))
    push!(sol_final, "delta" => abs.(dual.(delta)))
    push!(sol_final, "delta_bar" => abs.(dual.(delta_bar)))
    push!(sol_final, "mu_bar" => abs.(dual.(mu_bar)))
    push!(sol_final, "mure" => abs.(dual.(mure)))
    push!(sol_final, "y_tilde" => value.(y_tilde))
    push!(sol_final, "x_tilde" => value.(x_tilde))
    push!(sol_final, "f_tildep" => value.(f_tildep))
    push!(sol_final, "f_tildem" => value.(f_tildem))
    push!(sol_final, "mu" => abs.(dual.(mu)))
    push!(sol_final, "rhop" => abs.(dual.(rhop)))
    push!(sol_final, "rhom" => abs.(dual.(rhom)))
    push!(sol_final, "rhon_tildep" => abs.(dual.(rhon_tildep)))
    push!(sol_final, "rhon_tildem" => abs.(dual.(rhon_tildem)))
    push!(sol_final, "rhoz_tildep" => abs.(dual.(rhoz_tildep)))
    push!(sol_final, "rhoz_tildem" => abs.(dual.(rhoz_tildem)))
    push!(sol_final, "rho_rv" => abs.(dual.(rho_rv)))
    push!(sol_final, "phip" => abs.(dual.(phip)))
    push!(sol_final, "phim" => abs.(dual.(phim)))
    push!(sol_final, "psip" => abs.(dual.(psip)))
    push!(sol_final, "psim" => abs.(dual.(psim)))
    push!(sol_final, "mu_tilde" => abs.(dual.(mu_tilde)))
    push!(sol_final, "nu_tilde" => abs.(dual.(nu_tilde)))
    push!(sol_final, "lambdap" => abs.(dual.(lambdap)))
    push!(sol_final, "lambdam" => abs.(dual.(lambdam)))
    push!(sol_final, "Z_tilde" => value.(Z_tilde))
    push!(sol_final, "gamma_tilde" => abs.(dual.(gamma_tilde)))
    push!(sol_final, "z_tilde" => value.(z_tilde))

    return termination_status(m), objective_value(m), sol_final
end

# augmented problem with minimization of distance to previous solution
function primal_dual_augmented_oracle(sol::Dict)
    # Build capacity investment model
    m = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attributes(m, "OutputFlag" => 0,
        "NumericFocus" => 3, "Presolve" => 1, "Threads" => 4)
        
    # define primal variables
    @variable(m, pp[j=1:size(D, 1), z=1:length(Z)] >= 0)
    @variable(m, pm[j=1:size(D, 1), z=1:length(Z)] >= 0)
    @variable(m, rp[1:size(D, 1), 1:length(N)] >= 0)
    @variable(m, rm[1:size(D, 1), 1:length(N)] >= 0)
    @variable(m, y[1:length(I), 1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, yrv[1:length(I), 1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, yre[1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, y_bar[1:length(X_bar), 1:size(D, 1)] >= 0)
    @variable(m, yrv_bar[1:length(X_bar), 1:size(D, 1)] >= 0)
    @variable(m, x[i=1:length(I), z=1:length(Z)] >= 0) # investment in technology i
    @variable(m, x_bar[g=1:length(X_bar)] >= 0) # capacity of exisiting generator conserved
    @variable(m, s[j=1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, y_tilde[1:size(D, 1), 1:length(N)] >= 0) # power produced by tech i in slice j
    @variable(m, x_tilde[1:length(I), 1:length(N)] >= 0) # investment in technology i
    @variable(m, f_tildep[1:size(D, 1), l=1:length(L)] >= 0) # flows on line l in period j
    @variable(m, f_tildem[1:size(D, 1), l=1:length(L)] >= 0) # flows on line l in period j

    # define dual variables
    @variable(m, mu[1:length(I), 1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, mure[1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, mu_bar[1:length(X_bar), 1:size(D, 1)] >= 0)
    @variable(m, delta_bar[1:length(X_bar)] >= 0) # dual variable of the exisitng capacity constraint
    @variable(m, delta[1:length(I)] >= 0)
    @variable(m, rhop[j=1:size(D, 1), z=1:length(Z)] >= 0)
    @variable(m, rhom[j=1:size(D, 1), z=1:length(Z)] >= 0)
    @variable(m, rhon_tildep[1:size(D, 1), 1:length(N)] >= 0)
    @variable(m, rhon_tildem[1:size(D, 1), 1:length(N)] >= 0)
    @variable(m, rhoz_tildep[1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, rhoz_tildem[1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, rho_rv[1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, phip[1:size(D, 1)] >= 0)
    @variable(m, phim[1:size(D, 1)] >= 0)
    @variable(m, psip[1:size(D, 1), 1:length(L)] >= 0)
    @variable(m, psim[1:size(D, 1), 1:length(L)] >= 0)
    @variable(m, mu_tilde[1:size(D,1), 1:length(N)] >= 0)
    @variable(m, nu_tilde[1:length(Z)] >= 0)
    @variable(m, lambdap[1:size(D, 1), 1:length(L)] >= 0)
    @variable(m, lambdam[1:size(D, 1), 1:length(L)] >= 0)

    @variable(m, Z_tilde[n=1:length(N)] >= 0)
    @variable(m, gamma_tilde[1:size(D, 1), 1:length(N)] >= 0)
    @variable(m, z_tilde[1:size(D, 1), 1:length(N)] >= 0)

    # primal constraints
    @constraint(m, Fmu_bar[g=1:length(X_bar), j=1:size(D, 1)], x_bar[g] - y_bar[g,j] - yrv_bar[g,j] >= 0)
    @constraint(m, Fmure[j=1:size(D, 1), z=1:length(Z)], sum(R[j,n] for n in Nz[z]) - yre[j,z] >= 0)
    @constraint(m, Fdelta_bar[g=1:length(X_bar)], X_bar[g] - x_bar[g] >= 0)
    @constraint(m, Fdelta[i=1:length(I)], X[i] - sum(x[i,z] for z in 1:length(Z)) >= 0)
    @constraint(m, Frhop[j=1:size(D, 1), z=1:length(Z)], -(pp[j,z] - pm[j,z]) + yre[j,z] + sum(y[i,j,z] for i=1:length(I)) + sum(y_bar[g,j] for g=Gz[z]) - sum(D[j,n] for n in Nz[z]) + s[j,z] >= 0)
    @constraint(m, Frhom[j=1:size(D, 1), z=1:length(Z)], (pp[j,z] - pm[j,z]) - yre[j,z] - sum(y[i,j,z] for i=1:length(I)) - sum(y_bar[g,j] for g=Gz[z]) + sum(D[j,n] for n in Nz[z]) - s[j,z] >= 0)
    @constraint(m, Frhoz_tildep[j=1:size(D, 1), z=1:length(Z)], -(pp[j,z] - pm[j,z]) + sum(z_tilde[j,n] for n in Nz[z]) + sum(y_tilde[j,n] for n in Nz[z]) - sum(D[j,n] for n in Nz[z]) >= 0)
    @constraint(m, Frhoz_tildem[j=1:size(D, 1), z=1:length(Z)], (pp[j,z] - pm[j,z]) - sum(z_tilde[j,n] for n in Nz[z]) - sum(y_tilde[j,n] for n in Nz[z]) + sum(D[j,n] for n in Nz[z]) >= 0)
    @constraint(m, Frhon_tildep[j=1:size(D, 1), n=1:length(N)], (rp[j,n] - rm[j,n]) - z_tilde[j,n] - y_tilde[j,n] + D[j,n] >= 0)
    @constraint(m, Frhon_tildem[j=1:size(D, 1), n=1:length(N)], -(rp[j,n] - rm[j,n]) + z_tilde[j,n] + y_tilde[j,n] - D[j,n] >= 0)
    @constraint(m, Frho_rv[j=1:size(D,1), z=1:length(Z)], -Rv[z] + sum(yrv[i,j,z] for i=1:2) + sum(yrv_bar[g,j] for g in Gz[z] if isfastG[g]) >= 0)
    @constraint(m, Fphip[j=1:size(D, 1)], sum((rp[j,n] - rm[j,n]) for n=1:length(N)) >= 0)
    @constraint(m, Fphim[j=1:size(D, 1)], -sum((rp[j,n] - rm[j,n]) for n=1:length(N)) >= 0)
    @constraint(m, Fpsip[j=1:size(D, 1), l=1:length(L)], (f_tildep[j,l] - f_tildem[j,l]) - sum(PTDF[n,l]*(rp[j,n] - rm[j,n]) for n=1:length(N)) >= 0)
    @constraint(m, Fpsim[j=1:size(D, 1), l=1:length(L)], -(f_tildep[j,l] - f_tildem[j,l]) + sum(PTDF[n,l]*(rp[j,n] - rm[j,n]) for n=1:length(N)) >= 0)
    @constraint(m, Fmu[i=1:length(I), j=1:size(D, 1), z=1:length(Z)], x[i,z] - y[i,j,z] - yrv[i,j,z] >= 0)
    @constraint(m, Fmu_tilde[j=1:size(D, 1), n=1:length(N)], R[j,n] + sum(x_bar[g] for g=Gn[n]) + sum(x_tilde[i,n] for i=1:length(I)) - y_tilde[j,n] >= 0)
    @constraint(m, Fnu_tilde[z=1:length(Z)], sum(x[i,z] for i in 1:length(I)) - sum(x_tilde[i,n] for i in 1:length(I), n in Nz[z]) >= 0)
    @constraint(m, Flambdap[j=1:size(D,1), l=1:length(L)], TC[l] - (f_tildep[j,l] - f_tildem[j,l]) >= 0)
    @constraint(m, Flambdam[j=1:size(D,1), l=1:length(L)], TC[l] + (f_tildep[j,l] - f_tildem[j,l]) >= 0)
    @constraint(m, Fgamma_tilde[j=1:size(D, 1), n=1:length(N)], Z_tilde[n] - z_tilde[j,n] >= 0)

    # define dual constraints
    @constraint(m, Fy[i=1:length(I), j=1:size(D, 1), z=1:length(Z)], mu[i,j,z] + DT[j]*MC[i] - (rhop[j,z] - rhom[j,z]) >= 0)
    @constraint(m, Fyrv[i=1:length(I), j=1:size(D, 1), z=1:length(Z)], mu[i,j,z] + DT[j]*MC[i] - rho_rv[j,z]*(i in [1, 2]) >= 0)
    @constraint(m, Fyre[j=1:size(D, 1), z=1:length(Z)], mure[j,z] - (rhop[j,z] - rhom[j,z]) >= 0)
    @constraint(m, Fy_bar[g=1:length(X_bar), j=1:size(D, 1)], mu_bar[g,j] + DT[j]*MC_ex[g] - (rhop[j,G_Zidx[g]] - rhom[j,G_Zidx[g]]) >= 0)
    @constraint(m, Fyrv_bar[g=1:length(X_bar), j=1:size(D, 1)], mu_bar[g,j] + DT[j]*MC_ex[g] - rho_rv[j,G_Zidx[g]]*isfastG[g] >= 0)
    @constraint(m, Fx[i=1:length(I), z=1:length(Z)], I[i] + FC[i] - sum(mu[i,j,z] for j=1:size(D, 1)) + delta[i] - nu_tilde[z] + sol["nu_tilde"][z]  >= 0)
    @constraint(m, Fx_bar[g=1:length(X_bar)], FC_ex[g] - sum(mu_bar[g,j] for j=1:size(D, 1)) + delta_bar[g] - sum(mu_tilde[j,G_Nidx[g]] for j=1:size(D, 1)) + sum(sol["mu_tilde"][j,G_Nidx[g]] for j=1:size(D, 1)) >= 0)
    @constraint(m, Fs[j=1:size(D, 1), z=1:length(Z)], DT[j]*VOLL - (rhop[j,z] - rhom[j,z]) >= 0)
    @constraint(m, Fy_tilde[j=1:size(D, 1), n=1:length(N)], mu_tilde[j,n] - (rhoz_tildep[j,Zn[n]] - rhoz_tildem[j,Zn[n]]) + (rhon_tildep[j,n] - rhon_tildem[j,n]) >= 0)
    @constraint(m, Fx_tilde[i=1:length(I), n=1:length(N)], nu_tilde[Zn[n]] - sum(mu_tilde[j,n] for j=1:size(D, 1)) >= 0)
    @constraint(m, Fpp[j=1:size(D, 1), z=1:length(Z)], (rhop[j,z] - rhom[j,z]) + (rhoz_tildep[j,z] - rhoz_tildem[j,z]) >= 0)
    @constraint(m, Fpm[j=1:size(D, 1), z=1:length(Z)], -(rhop[j,z] - rhom[j,z]) - (rhoz_tildep[j,z] - rhoz_tildem[j,z]) >= 0)
    @constraint(m, Frp[j=1:size(D, 1), n=1:length(N)], -(rhon_tildep[j,n] - rhon_tildem[j,n]) - (phip[j] - phim[j]) + sum(PTDF[n,l]*(psip[j,l] - psim[j,l]) for l=1:length(L)) >= 0)
    @constraint(m, Frm[j=1:size(D, 1), n=1:length(N)], (rhon_tildep[j,n] - rhon_tildem[j,n]) + (phip[j] - phim[j]) - sum(PTDF[n,l]*(psip[j,l] - psim[j,l]) for l=1:length(L)) >= 0)
    @constraint(m, Ff_tildep[j=1:size(D, 1), l=1:length(L)], -(psip[j,l] - psim[j,l]) + (lambdap[j,l] - lambdam[j,l]) >= 0)
    @constraint(m, Ff_tildem[j=1:size(D, 1), l=1:length(L)], (psip[j,l] - psim[j,l]) - (lambdap[j,l] - lambdam[j,l]) >= 0)
    @constraint(m, FZ_tilde[n=1:length(N)], I_tilde - sum(gamma_tilde[j,n] for j=1:size(D, 1)) >= 0)
    @constraint(m, Fz_tilde[j=1:size(D, 1), n=1:length(N)], DT[j]*MC_tilde + gamma_tilde[j,n] + (rhon_tildep[j,n] - rhon_tildem[j,n]) - (rhoz_tildep[j,Zn[n]] - rhoz_tildem[j,Zn[n]]) >= 0)

    @expression(m, original_obj, sum(DT[j]*(VOLL*s[j,z] + sum(MC[i]*(y[i,j,z] + yrv[i,j,z])
        for i=1:length(I))) for j=1:size(D, 1), z=1:length(Z)) +
        sum(DT[j]*MC_ex[g]*(y_bar[g,j] + yrv_bar[g,j]) for g=1:length(X_bar), j=1:size(D, 1)) +
        sum((I[i] + FC[i])*x[i,z] for i=1:length(I), z=1:length(Z)) +
        sum(FC_ex[g]*x_bar[g] for g=1:length(X_bar)) +
        sum(DT[j]*MC_tilde*z_tilde[j,n] for j=1:size(D, 1), n=1:length(N)) +
        I_tilde*sum(Z_tilde[n] for n=1:length(N)))

    @expression(m, augmented_term, sum(x[i,z]*sol["nu_tilde"][z]
        for i=1:length(I), z=1:length(Z)) +
        sum(x_bar[g]*(sum(sol["mu_tilde"][j,G_Nidx[g]] for j=1:size(D, 1)))
        for g=1:length(X_bar)))

    @expression(m, dual_obj,
        -sum(mure[j,z]*sum(R[j,n] for n=Nz[z]) for j=1:size(D,1), z=1:length(Z)) -
        sum(mu_tilde[j,n]*R[j,n] for j=1:size(D, 1), n=1:length(N)) -
        sum(X_bar[g]*delta_bar[g] for g=1:length(X_bar)) -
        sum(X[i]*delta[i] for i=1:length(I)) +
        sum(sum(D[j,n] for n=Nz[z])*(rhop[j,z] - rhom[j,z]) for j=1:size(D, 1), z=1:length(Z)) +
        sum(sum(D[j,n] for n=Nz[z])*(rhoz_tildep[j,z] - rhoz_tildem[j,z]) for j=1:size(D, 1), z=1:length(Z)) -
        sum(D[j,n]*(rhon_tildep[j,n] - rhon_tildem[j,n]) for j=1:size(D, 1), n=1:length(N)) +
        sum(Rv[z]*rho_rv[j,z] for j=1:size(D,1), z=1:length(Z)) -
        sum(TC[l]*sum(lambdap[j,l] + lambdam[j,l] for j=1:size(D, 1)) for l=1:length(L)))

    @constraint(m, primal_dual, original_obj + augmented_term <= dual_obj)

    @expression(m, dist_primal,
        sum((yre[j,z]-round(sol["yre"][j,z], digits=4))^2 for j=1:size(D, 1), z=1:length(Z)) +
        sum((y[i,j,z]-round(sol["y"][i,j,z], digits=4))^2 for i=1:length(I), j=1:size(D, 1), z=1:length(Z)) +
        sum((yrv[i,j,z]-round(sol["yrv"][i,j,z], digits=4))^2 for i=1:length(I), j=1:size(D, 1), z=1:length(Z)) +
        sum((x[i,z]-round(sol["x"][i,z], digits=4))^2 for i=1:length(I), z=1:length(Z)) +
        sum((s[j,z]-round(sol["s"][j,z], digits=4))^2 for j=1:size(D, 1), z=1:length(Z)) +
        sum((y_bar[g,j]-round(sol["y_bar"][g,j], digits=4))^2 for g=1:length(X_bar), j=1:size(D, 1)) +
        sum((yrv_bar[g,j]-round(sol["yrv_bar"][g,j], digits=4))^2 for g=1:length(X_bar), j=1:size(D, 1)) +
        sum((x_bar[g]-round(sol["x_bar"][g], digits=4))^2 for g=1:length(X_bar)) +
        sum((y_tilde[j,n]-round(sol["y_tilde"][j,n], digits=4))^2 for j=1:size(D, 1), n=1:length(N)) +
        sum((x_tilde[i,n]-round(sol["x_tilde"][i,n], digits=4))^2 for i=1:length(I), n=1:length(N)) +
        sum((f_tildep[j,l]-round(sol["f_tildep"][j,l], digits=4))^2 for j=1:size(D,1), l=1:length(L)) +
        sum((f_tildem[j,l]-round(sol["f_tildem"][j,l], digits=4))^2 for j=1:size(D,1), l=1:length(L)) +
        sum((pp[j,z]-round(sol["pp"][j,z], digits=4))^2 for j=1:size(D, 1), z=1:length(Z)) +
        sum((pm[j,z]-round(sol["pm"][j,z], digits=4))^2 for j=1:size(D, 1), z=1:length(Z)) +
        sum((rp[j,n]-round(sol["rp"][j,n], digits=4))^2 for j=1:size(D, 1), n=1:length(N)) +
        sum((rm[j,n]-round(sol["rm"][j,n], digits=4))^2 for j=1:size(D, 1), n=1:length(N)) +
        sum((Z_tilde[n]-round(sol["Z_tilde"][n], digits=4))^2 for n=1:length(N)) +
        sum((z_tilde[j,n]-round(sol["z_tilde"][j,n], digits=4))^2 for j=1:size(D,1), n=1:length(N))
    )

    @expression(m, dist_dual,
        sum((delta_bar[g]-round(sol["delta_bar"][g], digits=4))^2 for g=1:length(X_bar)) +
        sum((delta[i]-round(sol["delta"][i], digits=4))^2 for i=1:length(I)) +
        sum((rhop[j,z]-round(sol["rhop"][j,z], digits=4))^2 for j=1:size(D, 1), z=1:length(Z)) +
        sum((rhom[j,z]-round(sol["rhom"][j,z], digits=4))^2 for j=1:size(D, 1), z=1:length(Z)) +
        sum((rhoz_tildep[j,z]-round(sol["rhoz_tildep"][j,z], digits=4))^2 for j=1:size(D, 1), z=1:length(Z)) +
        sum((rhoz_tildem[j,z]-round(sol["rhoz_tildem"][j,z], digits=4))^2 for j=1:size(D, 1), z=1:length(Z)) +
        sum((rhon_tildep[j,n]-round(sol["rhon_tildep"][j,n], digits=4))^2 for j=1:size(D, 1), n=1:length(N)) +
        sum((rhon_tildem[j,n]-round(sol["rhon_tildem"][j,n], digits=4))^2 for j=1:size(D, 1), n=1:length(N)) +
        sum((rho_rv[j,z]-round(sol["rho_rv"][j,z], digits=4))^2 for j=1:size(D, 1), z=1:length(Z)) +
        sum((phip[j]-round(sol["phip"][j], digits=4))^2 for j=1:size(D, 1)) +
        sum((phim[j]-round(sol["phim"][j], digits=4))^2 for j=1:size(D, 1)) +
        sum((psip[j,l]-round(sol["psip"][j,l], digits=4))^2 for j=1:size(D, 1), l=1:length(L)) +
        sum((psim[j,l]-round(sol["psim"][j,l], digits=4))^2 for j=1:size(D, 1), l=1:length(L)) +
        sum((mu[i,j,z]-round(sol["mu"][i,j,z], digits=4))^2 for i=1:length(I), j=1:size(D, 1), z=1:length(Z)) +
        sum((mure[j,z]-round(sol["mure"][j,z], digits=4))^2 for j=1:size(D, 1), z=1:length(Z)) +
        sum((mu_bar[g,j]-round(sol["mu_bar"][g,j], digits=4))^2 for g=1:length(X_bar), j=1:size(D, 1)) +
        sum((mu_tilde[j,n]-round(sol["mu_tilde"][j,n], digits=4))^2 for j=1:size(D, 1), n=1:length(N)) +
        sum((nu_tilde[z]-round(sol["nu_tilde"][z], digits=4))^2 for z=1:length(Z)) +
        sum((lambdap[j,l]-round(sol["lambdap"][j,l], digits=4))^2 for j=1:size(D, 1), l=1:length(L)) +
        sum((lambdam[j,l]-round(sol["lambdam"][j,l], digits=4))^2 for j=1:size(D, 1), l=1:length(L)) +
        sum((gamma_tilde[j,n]-round(sol["gamma_tilde"][j,n], digits=4))^2 for j=1:size(D,1), n=1:length(N))
    )

    @objective(m, Min, dist_primal + dist_dual)

    optimize!(m)

    # return solution
    sol_final = Dict()
    push!(sol_final, "pp" => value.(pp))
    push!(sol_final, "pm" => value.(pm))
    push!(sol_final, "rp" => value.(rp))
    push!(sol_final, "rm" => value.(rm))
    push!(sol_final, "y" => value.(y))
    push!(sol_final, "yrv" => value.(yrv))
    push!(sol_final, "yre" => value.(yre))
    push!(sol_final, "x" => value.(x))
    push!(sol_final, "s" => value.(s))
    push!(sol_final, "y_bar" => value.(y_bar))
    push!(sol_final, "yrv_bar" => value.(yrv_bar))
    push!(sol_final, "x_bar" => value.(x_bar))
    push!(sol_final, "delta" => value.(delta))
    push!(sol_final, "delta_bar" => value.(delta_bar))
    push!(sol_final, "mu_bar" => value.(mu_bar))
    push!(sol_final, "y_tilde" => value.(y_tilde))
    push!(sol_final, "x_tilde" => value.(x_tilde))
    push!(sol_final, "f_tildep" => value.(f_tildep))
    push!(sol_final, "f_tildem" => value.(f_tildem))
    push!(sol_final, "mu" => value.(mu))
    push!(sol_final, "mure" => value.(mure))
    push!(sol_final, "rhop" => value.(rhop))
    push!(sol_final, "rhom" => value.(rhom))
    push!(sol_final, "rhon_tildep" => value.(rhon_tildep))
    push!(sol_final, "rhon_tildem" => value.(rhon_tildem))
    push!(sol_final, "rhoz_tildep" => value.(rhoz_tildep))
    push!(sol_final, "rhoz_tildem" => value.(rhoz_tildem))
    push!(sol_final, "rho_rv" => value.(rho_rv))
    push!(sol_final, "phip" => value.(phip))
    push!(sol_final, "phim" => value.(phim))
    push!(sol_final, "psip" => value.(psip))
    push!(sol_final, "psim" => value.(psim))
    push!(sol_final, "mu_tilde" => value.(mu_tilde))
    push!(sol_final, "nu_tilde" => value.(nu_tilde))
    push!(sol_final, "lambdap" => value.(lambdap))
    push!(sol_final, "lambdam" => value.(lambdam))
    push!(sol_final, "Z_tilde" => value.(Z_tilde))
    push!(sol_final, "gamma_tilde" => value.(gamma_tilde))
    push!(sol_final, "z_tilde" => value.(z_tilde))

    return termination_status(m), objective_value(m), sol_final
end

function starting_point_zero()
    sol_final = Dict()
    push!(sol_final, "pp" => zeros(size(D, 1), length(Z)))
    push!(sol_final, "pm" => zeros(size(D, 1), length(Z)))
    push!(sol_final, "rp" => zeros(size(D, 1), length(N)))
    push!(sol_final, "rm" => zeros(size(D, 1), length(N)))
    push!(sol_final, "y" => zeros(length(I), size(D, 1), length(Z)))
    push!(sol_final, "yrv" => zeros(length(I), size(D, 1), length(Z)))
    push!(sol_final, "yre" => zeros(size(D, 1), length(Z)))
    push!(sol_final, "x" => zeros(length(I), length(Z)))
    push!(sol_final, "s" => zeros(size(D, 1), length(Z)))
    push!(sol_final, "y_bar" => zeros(length(X_bar), size(D, 1)))
    push!(sol_final, "yrv_bar" => zeros(length(X_bar), size(D, 1)))
    push!(sol_final, "x_bar" => zeros(length(X_bar)))
    push!(sol_final, "delta_bar" => zeros(length(X_bar)))
    push!(sol_final, "delta" => zeros(length(I)))
    push!(sol_final, "mu_bar" => zeros(length(X_bar), size(D, 1)))
    push!(sol_final, "y_tilde" => zeros(size(D, 1), length(N)))
    push!(sol_final, "x_tilde" => zeros(length(I), length(N)))
    push!(sol_final, "f_tildep" => zeros(size(D, 1), length(L)))
    push!(sol_final, "f_tildem" => zeros(size(D, 1), length(L)))
    push!(sol_final, "mu" => zeros(length(I), size(D, 1), length(Z)))
    push!(sol_final, "mure" => zeros(size(D, 1), length(Z)))
    push!(sol_final, "rhop" => zeros(size(D, 1), length(Z)))
    push!(sol_final, "rhom" => zeros(size(D, 1), length(Z)))
    push!(sol_final, "rhon_tildep" => zeros(size(D, 1), length(N)))
    push!(sol_final, "rhon_tildem" => zeros(size(D, 1), length(N)))
    push!(sol_final, "rhoz_tildep" => zeros(size(D, 1), length(Z)))
    push!(sol_final, "rhoz_tildem" => zeros(size(D, 1), length(Z)))
    push!(sol_final, "rho_rv" => zeros(size(D, 1), length(Z)))
    push!(sol_final, "phip" => zeros(size(D, 1)))
    push!(sol_final, "phim" => zeros(size(D, 1)))
    push!(sol_final, "psip" => zeros(size(D, 1), length(L)))
    push!(sol_final, "psim" => zeros(size(D, 1), length(L)))
    push!(sol_final, "mu_tilde" => zeros(size(D, 1), length(N)))
    push!(sol_final, "nu_tilde" => zeros(length(Z)))
    push!(sol_final, "lambdap" => zeros(size(D, 1), length(L)))
    push!(sol_final, "lambdam" => zeros(size(D, 1), length(L)))
    push!(sol_final, "Z_tilde" => zeros(length(N)))
    push!(sol_final, "gamma_tilde" => zeros(size(D,1), length(N)))
    push!(sol_final, "z_tilde" => zeros(size(D,1), length(N)))

    return MOI.OPTIMAL, sol_final
end

# computes the distance between two solutions
function distance(sol1::Dict, sol2::Dict)

    dist = 0
    for k in keys(sol1)
        dist += sum((sol1[k] .- sol2[k]).^2)
    end
    return dist
end
