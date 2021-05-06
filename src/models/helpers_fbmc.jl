# function to check whether a given vector of values is a solution to the LCP
function check_solution(sol::Dict)::Bool
    issol = Vector{Bool}()
    # primal constraints
    for k in keys(sol)
        push!(issol, all(sol[k] .>= -TOL))
    end

    # define dual mappings
    Fy(mu, rhop, rhom) = collect(mu[i,j,z] + DT[j]*MC[i] - (rhop[j,z] - rhom[j,z]) for i=1:length(I), j=1:size(D, 1), z=1:length(Z))
    Fyrv(mu, rho_rv) = collect(mu[i,j,z] + DT[j]*MC[i] - rho_rv[j,z]*(i in [1, 2]) for i=1:length(I), j=1:size(D, 1), z=1:length(Z))
    Fx(mu, delta) = collect(I[i] + FC[i] - sum(mu[i,j,z] for j=1:size(D, 1)) + delta[i] for i=1:length(I), z=1:length(Z))
    Fs(rhop, rhom) = collect(DT[j]*VOLL - (rhop[j,z] - rhom[j,z]) for j=1:size(D, 1), z=1:length(Z))
    Fmu(x, y, yrv) = collect(x[i,z] - y[i,j,z] - yrv[i,j,z] for i=1:length(I), j=1:size(D, 1), z=1:length(Z))

    Fy_bar(mu_bar, rhop, rhom) = collect(mu_bar[g,j] + DT[j]*MC_ex[g] - (rhop[j,G_Zidx[g]] - rhom[j,G_Zidx[g]]) for g=1:length(X_bar), j=1:size(D, 1))
    Fyrv_bar(mu_bar, rho_rv) = collect(mu_bar[g,j] + DT[j]*MC_ex[g] - rho_rv[j,G_Zidx[g]]*isfastG[g] for g=1:length(X_bar), j=1:size(D, 1))
    Fx_bar(mu_bar, delta_bar) = collect(FC_ex[g] - sum(mu_bar[g,j] for j=1:size(D, 1)) + delta_bar[g] for g=1:length(X_bar))
    Fdelta_bar(x_bar) = collect(X_bar[g] - x_bar[g] for g=1:length(X_bar))
    Fdelta(x) = collect(X[i] - sum(x[i,z] for z in 1:length(Z)) for i=1:length(I))
    Fmu_bar(x_bar, y_bar, yrv_bar) = collect(x_bar[g] - y_bar[g,j] - yrv_bar[g,j] for g=1:length(X_bar), j=1:size(D, 1))

    Fnu_tilde(x, x_tilde) = collect(sum(x[i,z] for i=1:length(I)) - sum(x_tilde[i,n] for i=1:length(I), n in Nz[z]) for z=1:length(Z))

    Frhop(pp, pm, y, y_bar, s) = collect(-(pp[j,z] - pm[j,z]) + sum(y[i,j,z] for i=1:length(I)) + (isempty(Gz[z]) ? 0 : sum(y_bar[g,j] for g=Gz[z])) - sum(D[j,n] for n in Nz[z]) + s[j,z] for j=1:size(D, 1), z=1:length(Z))
    Frhom(pp, pm, y, y_bar, s) = collect((pp[j,z] - pm[j,z]) - sum(y[i,j,z] for i=1:length(I)) - (isempty(Gz[z]) ? 0 : sum(y_bar[g,j] for g=Gz[z])) + sum(D[j,n] for n in Nz[z]) - s[j,z] for j=1:size(D, 1), z=1:length(Z))
    Frhoz_tildep(pp, pm, z_tilde, y_tilde) = collect(-(pp[j,z] - pm[j,z]) + sum(z_tilde[j,n] for n in Nz[z]) + sum(y_tilde[j,n] for n in Nz[z]) - sum(D[j,n] for n in Nz[z]) for j=1:size(D, 1), z=1:length(Z))
    Frhoz_tildem(pp, pm, z_tilde, y_tilde) = collect((pp[j,z] - pm[j,z]) - sum(z_tilde[j,n] for n in Nz[z]) - sum(y_tilde[j,n] for n in Nz[z]) + sum(D[j,n] for n in Nz[z]) for j=1:size(D, 1), z=1:length(Z))
    Frho_rv(yrv, yrv_bar) = collect(-Rv[z] + sum(yrv[i,j,z] for i=1:2) + ((isempty(Gz[z]) || !any(isfastG[g] for g in Gz[z])) ? 0 : sum(yrv_bar[g,j] for g in Gz[z] if isfastG[g])) for j=1:size(D,1), z=1:length(Z))
    Fy_tilde(mu_tilde, rhoz_tildep, rhoz_tildem, rhon_tildep, rhon_tildem) = collect(mu_tilde[j,n] - (rhoz_tildep[j,Zn[n]] - rhoz_tildem[j,Zn[n]]) + (rhon_tildep[j,n] - rhon_tildem[j,n]) for j=1:size(D, 1), n=1:length(N))
    Fx_tilde(nu_tilde, mu_tilde) = collect(nu_tilde[Zn[n]] - sum(mu_tilde[j,n] for j=1:size(D, 1)) for i=1:length(I), n=1:length(N))
    Fpp(rhop, rhom, rhoz_tildep, rhoz_tildem) = collect((rhop[j,z] - rhom[j,z]) + (rhoz_tildep[j,z] - rhoz_tildem[j,z]) for j=1:size(D, 1), z=1:length(Z))
    Fpm(rhop, rhom, rhoz_tildep, rhoz_tildem) = collect(-(rhop[j,z] - rhom[j,z]) - (rhoz_tildep[j,z] - rhoz_tildem[j,z]) for j=1:size(D, 1), z=1:length(Z))
    Frp(rhon_tildep, rhon_tildem, phip, phim, psip, psim) = collect(-(rhon_tildep[j,n] - rhon_tildem[j,n]) - (phip[j] - phim[j]) + sum(PTDF[n,l]*(psip[j,l] - psim[j,l]) for l=1:length(L)) for j=1:size(D, 1), n=1:length(N))
    Frm(rhon_tildep, rhon_tildem, phip, phim, psip, psim) = collect((rhon_tildep[j,n] - rhon_tildem[j,n]) + (phip[j] - phim[j]) - sum(PTDF[n,l]*(psip[j,l] - psim[j,l]) for l=1:length(L)) for j=1:size(D, 1), n=1:length(N))
    Ff_tildep(psip, psim, lambdap, lambdam) = collect(-(psip[j,l] - psim[j,l]) + (lambdap[j,l] - lambdam[j,l]) for j=1:size(D, 1), l=1:length(L))
    Ff_tildem(psip, psim, lambdap, lambdam) = collect((psip[j,l] - psim[j,l]) - (lambdap[j,l] - lambdam[j,l]) for j=1:size(D, 1), l=1:length(L))
    Fmu_tilde(x_bar, x_tilde, y_tilde) = collect(R[j,n] + (isempty(Gn[n]) ? 0 : sum(x_bar[g] for g=Gn[n])) + sum(x_tilde[i,n] for i=1:length(I)) - y_tilde[j,n] for j=1:size(D, 1), n=1:length(N))
    Frhon_tildep(rp, rm, z_tilde, y_tilde) = collect((rp[j,n] - rm[j,n]) - z_tilde[j,n] - y_tilde[j,n] + D[j,n] for j=1:size(D, 1), n=1:length(N))
    Frhon_tildem(rp, rm, z_tilde, y_tilde) = collect(-(rp[j,n] - rm[j,n]) + z_tilde[j,n] + y_tilde[j,n] - D[j,n] for j=1:size(D, 1), n=1:length(N))
    Fphip(rp, rm) = collect(sum((rp[j,n] - rm[j,n]) for n=1:length(N)) for j=1:size(D, 1))
    Fphim(rp, rm) = collect(sum(-(rp[j,n] - rm[j,n]) for n=1:length(N)) for j=1:size(D, 1))
    Fpsip(f_tildep, f_tildem, rp, rm) = collect((f_tildep[j,l] - f_tildem[j,l]) - sum(PTDF[n,l]*(rp[j,n] - rm[j,n]) for n=1:length(N)) for j=1:size(D, 1), l=1:length(L))
    Fpsim(f_tildep, f_tildem, rp, rm) = collect(-(f_tildep[j,l] - f_tildem[j,l]) + sum(PTDF[n,l]*(rp[j,n] - rm[j,n]) for n=1:length(N)) for j=1:size(D, 1), l=1:length(L))
    Flambdap(f_tildep, f_tildem) = collect(TC[l] - (f_tildep[j,l] - f_tildem[j,l]) for j=1:size(D,1), l=1:length(L))
    Flambdam(f_tildep, f_tildem) = collect(TC[l] + (f_tildep[j,l] - f_tildem[j,l]) for j=1:size(D,1), l=1:length(L))

    FZ_tilde(gamma_tilde) = collect(I_tilde - sum(gamma_tilde[j,n] for j=1:size(D, 1)) for n=1:length(N))
    Fz_tilde(gamma_tilde, rhon_tildep, rhon_tildem, rhoz_tildep, rhoz_tildem) = collect(DT[j]*MC_tilde + gamma_tilde[j,n] + (rhon_tildep[j,n] - rhon_tildem[j,n]) - (rhoz_tildep[j,Zn[n]] - rhoz_tildem[j,Zn[n]]) for j=1:size(D, 1), n=1:length(N))
    Fgamma_tilde(Z_tilde, z_tilde) = collect(Z_tilde[n] - z_tilde[j,n] for j=1:size(D, 1), n=1:length(N))
    # dual constraints
    push!(issol, all(Fy(sol["mu"], sol["rhop"], sol["rhom"]) .>= -TOL))
    push!(issol, all(Fyrv(sol["mu"], sol["rho_rv"]) .>= -TOL))
    push!(issol, all(Fx(sol["mu"], sol["delta"]) .>= -TOL))
    push!(issol, all(Fs(sol["rhop"], sol["rhom"]) .>= -TOL))
    push!(issol, all(Fmu(sol["x"], sol["y"], sol["yrv"]) .>= -TOL))

    push!(issol, all(Fy_bar(sol["mu_bar"], sol["rhop"], sol["rhom"]) .>= -TOL))
    push!(issol, all(Fyrv_bar(sol["mu_bar"], sol["rho_rv"]) .>= -TOL))
    push!(issol, all(Fx_bar(sol["mu_bar"], sol["delta_bar"]) .>= -TOL))
    push!(issol, all(Fdelta(sol["x"]) .>= -TOL))
    push!(issol, all(Fdelta_bar(sol["x_bar"]) .>= -TOL))
    push!(issol, all(Fmu_bar(sol["x_bar"], sol["y_bar"], sol["yrv_bar"]) .>= -TOL))

    push!(issol, all(Frhop(sol["pp"], sol["pm"], sol["y"], sol["y_bar"], sol["s"]) .>= -TOL))
    push!(issol, all(Frhom(sol["pp"], sol["pm"], sol["y"], sol["y_bar"], sol["s"]) .>= -TOL))
    push!(issol, all(Frhoz_tildep(sol["pp"], sol["pm"], sol["z_tilde"], sol["y_tilde"]) .>= -TOL))
    push!(issol, all(Frhoz_tildem(sol["pp"], sol["pm"], sol["z_tilde"], sol["y_tilde"]) .>= -TOL))
    push!(issol, all(Frho_rv(sol["yrv"], sol["yrv_bar"]) .>= -TOL))
    push!(issol, all(Fy_tilde(sol["mu_tilde"], sol["rhoz_tildep"], sol["rhoz_tildem"], sol["rhon_tildep"], sol["rhon_tildem"]) .>= -TOL))
    push!(issol, all(Fx_tilde(sol["nu_tilde"], sol["mu_tilde"]) .>= -TOL))
    push!(issol, all(Fpp(sol["rhop"], sol["rhom"], sol["rhoz_tildep"], sol["rhoz_tildem"]) .>= -TOL))
    push!(issol, all(Fpm(sol["rhop"], sol["rhom"], sol["rhoz_tildep"], sol["rhoz_tildem"]) .>= -TOL))
    push!(issol, all(Frp(sol["rhon_tildep"], sol["rhon_tildem"], sol["phip"], sol["phim"], sol["psip"], sol["psim"]) .>= -TOL))
    push!(issol, all(Frm(sol["rhon_tildep"], sol["rhon_tildem"], sol["phip"], sol["phim"], sol["psip"], sol["psim"]) .>= -TOL))
    push!(issol, all(Ff_tildep(sol["psip"], sol["psim"], sol["lambdap"], sol["lambdam"]) .>= -TOL))
    push!(issol, all(Ff_tildem(sol["psip"], sol["psim"], sol["lambdap"], sol["lambdam"]) .>= -TOL))

    push!(issol, all(Fmu_tilde(sol["x_bar"], sol["x_tilde"], sol["y_tilde"]) .>= -TOL))
    push!(issol, all(Frhon_tildep(sol["rp"], sol["rm"], sol["z_tilde"], sol["y_tilde"]) .>= -TOL))
    push!(issol, all(Frhon_tildem(sol["rp"], sol["rm"], sol["z_tilde"], sol["y_tilde"]) .>= -TOL))
    push!(issol, all(Fphip(sol["rp"], sol["rm"]) .>= -TOL))
    push!(issol, all(Fphim(sol["rp"], sol["rm"]) .>= -TOL))
    push!(issol, all(Fpsip(sol["f_tildep"], sol["f_tildem"], sol["rp"], sol["rm"], ) .>= -TOL))
    push!(issol, all(Fpsim(sol["f_tildep"], sol["f_tildem"], sol["rp"], sol["rm"], ) .>= -TOL))
    push!(issol, all(Flambdap(sol["f_tildep"], sol["f_tildem"]) .>= -TOL))
    push!(issol, all(Flambdam(sol["f_tildep"], sol["f_tildem"]) .>= -TOL))

    push!(issol, all(FZ_tilde(sol["gamma_tilde"]) .>= -TOL))
    push!(issol, all(Fz_tilde(sol["gamma_tilde"], sol["rhon_tildep"], sol["rhon_tildem"], sol["rhoz_tildep"], sol["rhoz_tildem"]) .>= -TOL))
    push!(issol, all(Fgamma_tilde(sol["Z_tilde"], sol["z_tilde"]) .>= -TOL))

    # complementarities
    push!(issol, all(isapprox.(Fy(sol["mu"], sol["rhop"], sol["rhom"]) .* sol["y"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fyrv(sol["mu"], sol["rho_rv"]) .* sol["yrv"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fx(sol["mu"], sol["delta"]) .* sol["x"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fs(sol["rhop"], sol["rhom"]) .* sol["s"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fmu(sol["x"], sol["y"], sol["yrv"]) .* sol["mu"], 0, atol=TOL)))

    push!(issol, all(isapprox.(Fy_bar(sol["mu_bar"], sol["rhop"], sol["rhom"]) .* sol["y_bar"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fyrv_bar(sol["mu_bar"], sol["rho_rv"]) .* sol["yrv_bar"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fx_bar(sol["mu_bar"], sol["delta_bar"]) .* sol["x_bar"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fdelta(sol["x"]) .* sol["delta"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fdelta_bar(sol["x_bar"]) .* sol["delta_bar"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fmu_bar(sol["x_bar"], sol["y_bar"], sol["yrv_bar"]) .* sol["mu_bar"], 0, atol=TOL)))

    push!(issol, all(isapprox.(Frhop(sol["pp"], sol["pm"], sol["y"], sol["y_bar"], sol["s"]) .* sol["rhop"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Frhom(sol["pp"], sol["pm"], sol["y"], sol["y_bar"], sol["s"]) .* sol["rhom"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Frhoz_tildep(sol["pp"], sol["pm"], sol["z_tilde"], sol["y_tilde"]) .* sol["rhoz_tildep"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Frhoz_tildem(sol["pp"], sol["pm"], sol["z_tilde"], sol["y_tilde"]) .* sol["rhoz_tildem"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Frho_rv(sol["yrv"], sol["yrv_bar"]) .* sol["rho_rv"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fy_tilde(sol["mu_tilde"], sol["rhoz_tildep"], sol["rhoz_tildem"], sol["rhon_tildep"], sol["rhon_tildem"]) .* sol["y_tilde"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fx_tilde(sol["nu_tilde"], sol["mu_tilde"]) .* sol["x_tilde"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fpp(sol["rhop"], sol["rhom"], sol["rhoz_tildep"], sol["rhoz_tildem"]) .* sol["pp"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fpm(sol["rhop"], sol["rhom"], sol["rhoz_tildep"], sol["rhoz_tildem"]) .* sol["pm"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Frp(sol["rhon_tildep"], sol["rhon_tildem"], sol["phip"], sol["phim"], sol["psip"], sol["psim"]) .* sol["rp"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Frm(sol["rhon_tildep"], sol["rhon_tildem"], sol["phip"], sol["phim"], sol["psip"], sol["psim"]) .* sol["rm"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Ff_tildep(sol["psip"], sol["psim"], sol["lambdap"], sol["lambdam"]) .* sol["f_tildep"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Ff_tildem(sol["psip"], sol["psim"], sol["lambdap"], sol["lambdam"]) .* sol["f_tildem"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fmu_tilde(sol["x_bar"], sol["x_tilde"], sol["y_tilde"]) .* sol["mu_tilde"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Frhon_tildep(sol["rp"], sol["rm"], sol["z_tilde"], sol["y_tilde"]) .* sol["rhon_tildep"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Frhon_tildem(sol["rp"], sol["rm"], sol["z_tilde"], sol["y_tilde"]) .* sol["rhon_tildem"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fphip(sol["rp"], sol["rm"]) .* sol["phip"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fphim(sol["rp"], sol["rm"]) .* sol["phim"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fpsip(sol["f_tildep"], sol["f_tildem"], sol["rp"], sol["rm"]) .* sol["psip"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fpsim(sol["f_tildep"], sol["f_tildem"], sol["rp"], sol["rm"]) .* sol["psim"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Flambdap(sol["f_tildep"], sol["f_tildem"]) .* sol["lambdap"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Flambdam(sol["f_tildep"], sol["f_tildem"]) .* sol["lambdam"], 0, atol=TOL)))
    push!(issol, all(isapprox.(FZ_tilde(sol["gamma_tilde"]) .* sol["Z_tilde"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fz_tilde(sol["gamma_tilde"], sol["rhon_tildep"], sol["rhon_tildem"], sol["rhoz_tildep"], sol["rhoz_tildem"]) .* sol["z_tilde"], 0, atol=TOL)))
    push!(issol, all(isapprox.(Fgamma_tilde(sol["Z_tilde"], sol["z_tilde"]) .* sol["gamma_tilde"], 0, atol=TOL)))

    # @show issol
    return all(issol)
end

# function to compute the determinant of the Jacobian at a particular equilibrium solution
function equilibirum_determinant(numbering::Dict, sol::Dict)
    tot_length = numbering["z_tilde"][size(D,1), length(N)]
    @show tot_length
    M = zeros(tot_length, tot_length)
    for j=1:size(D, 1), z=1:length(Z)
        gidx = numbering["pp"][j,z]
        if isapprox(sol["pp"][j,z], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["rhop"][j,z]] = 1
            M[gidx, numbering["rhom"][j,z]] = -1
            M[gidx, numbering["rhoz_tildep"][j,z]] = 1
            M[gidx, numbering["rhoz_tildem"][j,z]] = -1
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        gidx = numbering["pm"][j,z]
        if isapprox(sol["pm"][j,z], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["rhop"][j,z]] = -1
            M[gidx, numbering["rhom"][j,z]] = 1
            M[gidx, numbering["rhoz_tildep"][j,z]] = -1
            M[gidx, numbering["rhoz_tildem"][j,z]] = 1
        end
    end
    for j in 1:size(D, 1), n in 1:length(N)
        gidx = numbering["rp"][j,n]
        if isapprox(sol["rp"][j,n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["rhon_tildep"][j,n]] = -1
            M[gidx, numbering["rhon_tildem"][j,n]] = 1
            M[gidx, numbering["phip"][j]] = -1
            M[gidx, numbering["phim"][j]] = 1
            for l in 1:length(L)
                M[gidx, numbering["psip"][j,l]] = PTDF[n,l]
                M[gidx, numbering["psim"][j,l]] = -PTDF[n,l]
            end
        end
    end
    for j in 1:size(D, 1), n in 1:length(N)
        gidx = numbering["rm"][j,n]
        if isapprox(sol["rm"][j,n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["rhon_tildep"][j,n]] = 1
            M[gidx, numbering["rhon_tildem"][j,n]] = -1
            M[gidx, numbering["phip"][j]] = 1
            M[gidx, numbering["phim"][j]] = -1
            for l in 1:length(L)
                M[gidx, numbering["psip"][j,l]] = -PTDF[n,l]
                M[gidx, numbering["psim"][j,l]] = PTDF[n,l]
            end
        end
    end
    for i=1:length(I), j=1:size(D, 1), z=1:length(Z)
        gidx = numbering["y"][i,j,z]
        if isapprox(sol["y"][i,j,z], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["mu"][i,j,z]] = 1
            M[gidx, numbering["rhop"][j,z]] = -1
            M[gidx, numbering["rhom"][j,z]] = 1
        end
    end
    for i=1:length(I), z=1:length(Z)
        gidx = numbering["x"][i,z]
        if isapprox(sol["x"][i,z], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            for j=1:size(D,1)
                M[gidx, numbering["mu"][i,j,z]] = -1
            end
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        gidx = numbering["s"][j,z]
        if isapprox(sol["s"][j,z], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["rhop"][j,z]] = -1
            M[gidx, numbering["rhom"][j,z]] = 1
        end
    end
    for i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)
        gidx = numbering["y_bar"][i,j,n]
        if isapprox(sol["y_bar"][i,j,n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["mu_bar"][i,j,n]] = 1
            M[gidx, numbering["rhop"][j,Zn[n]]] = -1
            M[gidx, numbering["rhom"][j,Zn[n]]] = 1
        end
    end
    for i=1:length(MC_ex), n=1:length(N)
        gidx = numbering["x_bar"][i,n]
        if isapprox(sol["x_bar"][i,n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["delta"][i,n]] = 1
            for j=1:size(D,1)
                M[gidx, numbering["mu_bar"][i,j,n]] = -1
            end
        end
    end
    for i=1:length(MC_ex), n=1:length(N)
        gidx = numbering["delta"][i,n]
        if isapprox(sol["delta"][i,n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["x_bar"][i,n]] = -1
        end
    end
    for i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)
        gidx = numbering["mu_bar"][i,j,n]
        if isapprox(sol["mu_bar"][i,j,n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["x_bar"][i,n]] = 1
            M[gidx, numbering["y_bar"][i,j,n]] = -1
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        gidx = numbering["y_tilde"][j,n]
        if isapprox(sol["y_tilde"][j,n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["mu_tilde"][j,n]] = 1
            M[gidx, numbering["rhoz_tildep"][j,Zn[n]]] = -1
            M[gidx, numbering["rhoz_tildem"][j,Zn[n]]] = 1
            M[gidx, numbering["rhon_tildep"][j,n]] = 1
            M[gidx, numbering["rhon_tildem"][j,n]] = -1
        end
    end
    for i=1:length(I), n=1:length(N)
        gidx = numbering["x_tilde"][i,n]
        if isapprox(sol["x_tilde"][i,n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["nu_tilde"][i,Zn[n]]] = 1
            for j=1:size(D,1)
                M[gidx, numbering["mu_tilde"][j,n]] = -1
            end
        end
    end
    for j=1:size(D, 1), l=1:length(L)
        gidx = numbering["f_tildep"][j,l]
        if isapprox(sol["f_tildep"][j,l], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["psip"][j,l]] = -1
            M[gidx, numbering["psim"][j,l]] = 1
            M[gidx, numbering["lambdap"][j,l]] = TC[l]
            M[gidx, numbering["lambdam"][j,l]] = -TC[l]
        end
    end
    for j=1:size(D, 1), l=1:length(L)
        gidx = numbering["f_tildem"][j,l]
        if isapprox(sol["f_tildem"][j,l], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["psip"][j,l]] = 1
            M[gidx, numbering["psim"][j,l]] = -1
            M[gidx, numbering["lambdap"][j,l]] = -TC[l]
            M[gidx, numbering["lambdam"][j,l]] = TC[l]
        end
    end
    for i=1:length(I), j=1:size(D, 1), z=1:length(Z)
        gidx = numbering["mu"][i,j,z]
        if isapprox(sol["mu"][i,j,z], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["x"][i,z]] = 1
            M[gidx, numbering["y"][i,j,z]] = -1
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        gidx = numbering["rhop"][j,z]
        if isapprox(sol["rhop"][j,z], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["pp"][j,z]] = -1
            M[gidx, numbering["pm"][j,z]] = 1
            for i=1:length(I)
                M[gidx, numbering["y"][i,j,z]] = 1
            end
            M[gidx, numbering["s"][j,z]] = 1
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        gidx = numbering["rhom"][j,z]
        if isapprox(sol["rhom"][j,z], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["pp"][j,z]] = 1
            M[gidx, numbering["pm"][j,z]] = -1
            for i=1:length(I)
                M[gidx, numbering["y"][i,j,z]] = -1
            end
            M[gidx, numbering["s"][j,z]] = -1
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        gidx = numbering["rhon_tildep"][j,n]
        if isapprox(sol["rhon_tildep"][j,n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["rp"][j,n]] = 1
            M[gidx, numbering["rm"][j,n]] = -1
            M[gidx, numbering["z_tilde"][j,n]] = -1
            M[gidx, numbering["y_tilde"][j,n]] = -1
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        gidx = numbering["rhon_tildem"][j,n]
        if isapprox(sol["rhon_tildem"][j,n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["rp"][j,n]] = -1
            M[gidx, numbering["rm"][j,n]] = 1
            M[gidx, numbering["z_tilde"][j,n]] = 1
            M[gidx, numbering["y_tilde"][j,n]] = 1
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        gidx = numbering["rhoz_tildep"][j,z]
        if isapprox(sol["rhoz_tildep"][j,z], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["pp"][j,z]] = -1
            M[gidx, numbering["pm"][j,z]] = 1
            for n in Nz[z]
                M[gidx, numbering["y_tilde"][j,n]] = 1
            end
            for n in Nz[z]
                M[gidx, numbering["z_tilde"][j,n]] = 1
            end
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        gidx = numbering["rhoz_tildem"][j,z]
        if isapprox(sol["rhoz_tildem"][j,z], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["pp"][j,z]] = 1
            M[gidx, numbering["pm"][j,z]] = -1
            for n in Nz[z]
                M[gidx, numbering["y_tilde"][j,n]] = -1
            end
            for n in Nz[z]
                M[gidx, numbering["z_tilde"][j,n]] = -1
            end
        end
    end
    for j=1:size(D, 1)
        gidx = numbering["phip"][j]
        if isapprox(sol["phip"][j], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            for n=1:length(N)
                M[gidx,numbering["rp"][j,n]] = 1
                M[gidx,numbering["rm"][j,n]] = -1
            end
        end
    end
    for j=1:size(D, 1)
        gidx = numbering["phim"][j]
        if isapprox(sol["phim"][j], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            for n=1:length(N)
                M[gidx,numbering["rp"][j,n]] = -1
                M[gidx,numbering["rm"][j,n]] = 1
            end
        end
    end
    for j=1:size(D, 1), l=1:length(L)
        gidx = numbering["psip"][j,l]
        if isapprox(sol["psip"][j,l], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["f_tildep"][j,l]] = 1
            M[gidx, numbering["f_tildem"][j,l]] = -1
            for n=1:length(N)
                M[gidx, numbering["rp"][j,n]] = -PTDF[n,l]
                M[gidx, numbering["rm"][j,n]] = PTDF[n,l]
            end
        end
    end
    for j=1:size(D, 1), l=1:length(L)
        gidx = numbering["psim"][j,l]
        if isapprox(sol["psim"][j,l], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["f_tildep"][j,l]] = -1
            M[gidx, numbering["f_tildem"][j,l]] = 1
            for n=1:length(N)
                M[gidx, numbering["rp"][j,n]] = PTDF[n,l]
                M[gidx, numbering["rm"][j,n]] = -PTDF[n,l]
            end
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        gidx = numbering["mu_tilde"][j,n]
        if isapprox(sol["mu_tilde"][j,n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            for i=1:length(MC_ex)
                M[gidx, numbering["x_bar"][i,n]] = 1
            end
            for i=1:length(I)
                M[gidx, numbering["x_tilde"][n]] = 1
            end
            M[gidx, numbering["y_tilde"][j,n]] = -1
        end
    end
    for i=1:length(I), z=1:length(Z)
        gidx = numbering["nu_tilde"][i,z]
        if isapprox(sol["nu_tilde"][i,z], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["x"][i,z]] = 1
            for n in Nz[z]
                M[gidx, numbering["x_tilde"][i,n]] = -1
            end
        end
    end
    for j=1:size(D, 1), l=1:length(L)
        gidx = numbering["lambdap"][j,l]
        if isapprox(sol["lambdap"][j,l], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["f_tildep"][j,l]] = -1
            M[gidx, numbering["f_tildem"][j,l]] = 1
        end
    end
    for j=1:size(D, 1), l=1:length(L)
        gidx = numbering["lambdam"][j,l]
        if isapprox(sol["lambdam"][j,l], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["f_tildep"][j,l]] = 1
            M[gidx, numbering["f_tildem"][j,l]] = -1
        end
    end
    for n=1:length(N)
        gidx = numbering["Z_tilde"][n]
        if isapprox(sol["Z_tilde"][n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            for j=1:size(D,1)
                M[gidx, numbering["gamma_tilde"][j,n]] = -1
            end
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        gidx = numbering["gamma_tilde"][j,n]
        if isapprox(sol["gamma_tilde"][j,n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["Z_tilde"][n]] = 1
            M[gidx, numbering["z_tilde"][j,n]] = -1
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        gidx = numbering["z_tilde"][j,n]
        if isapprox(sol["z_tilde"][j,n], 0, atol=TOL)
            M[gidx, gidx] = 1
        else
            M[gidx, numbering["gamma_tilde"][j,n]] = 1
            M[gidx, numbering["rhon_tildep"][j,n]] = 1
            M[gidx, numbering["rhon_tildem"][j,n]] = -1
            M[gidx, numbering["rhoz_tildep"][j,Zn[n]]] = -1
            M[gidx, numbering["rhoz_tildem"][j,Zn[n]]] = 1
        end
    end
    @show any((sum(abs.(M), dims=2)) .== 0)

    return det(M)
end

# function to determine the 'alpha' of the solution, that is the indices of complementarity conditions for which the constraint (as opposed to the variable) is tight at 0.
function get_alpha(numbering::Dict, sol::Dict)
    alpha = Vector{Int}()
    for j=1:size(D, 1), z=1:length(Z)
        if !isapprox(sol["pp"][j,z], 0, atol=TOL)
            push!(alpha, numbering["pp"][j,z])
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        if !isapprox(sol["pm"][j,z], 0, atol=TOL)
            push!(alpha, numbering["pm"][j,z])
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        if !isapprox(sol["rp"][j,n], 0, atol=TOL)
            push!(alpha, numbering["rp"][j,n])
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        if !isapprox(sol["rm"][j,n], 0, atol=TOL)
            push!(alpha, numbering["rm"][j,n])
        end
    end
    for i=1:length(I), j=1:size(D, 1), z=1:length(Z)
        if !isapprox(sol["y"][i,j,z], 0, atol=TOL)
            push!(alpha, numbering["y"][i,j,z])
        end
    end
    for i=1:length(I), z=1:length(Z)
        if !isapprox(sol["x"][i,z], 0, atol=TOL)
            push!(alpha, numbering["x"][i,z])
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        if !isapprox(sol["s"][j,z], 0, atol=TOL)
            push!(alpha, numbering["s"][j,z])
        end
    end
    for i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)
        if !isapprox(sol["y_bar"][i,j,n], 0, atol=TOL)
            push!(alpha, numbering["y_bar"][i,j,n])
        end
    end
    for i=1:length(MC_ex), n=1:length(N)
        if !isapprox(sol["x_bar"][i,n], 0, atol=TOL)
            push!(alpha, numbering["x_bar"][i,n])
        end
    end
    for i=1:length(MC_ex), n=1:length(N)
        if !isapprox(sol["delta"][i,n], 0, atol=TOL)
            push!(alpha, numbering["delta"][i,n])
        end
    end
    for i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)
        if !isapprox(sol["mu_bar"][i,j,n], 0, atol=TOL)
            push!(alpha, numbering["mu_bar"][i,j,n])
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        if !isapprox(sol["y_tilde"][j,n], 0, atol=TOL)
            push!(alpha, numbering["y_tilde"][j,n])
        end
    end
    for i=1:length(I), n=1:length(N)
        if !isapprox(sol["x_tilde"][i,n], 0, atol=TOL)
            push!(alpha, numbering["x_tilde"][i,n])
        end
    end
    for j=1:size(D,1), l=1:length(L)
        if !isapprox(sol["f_tildep"][j,l], 0, atol=TOL)
            push!(alpha, numbering["f_tildep"][j,l])
        end
    end
    for j=1:size(D,1), l=1:length(L)
        if !isapprox(sol["f_tildem"][j,l], 0, atol=TOL)
            push!(alpha, numbering["f_tildem"][j,l])
        end
    end
    for i=1:length(I), j=1:size(D, 1), z=1:length(Z)
        if !isapprox(sol["mu"][i,j,z], 0, atol=TOL)
            push!(alpha, numbering["mu"][i,j,z])
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        if !isapprox(sol["rhop"][j,z], 0, atol=TOL)
            push!(alpha, numbering["rhop"][j,z])
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        if !isapprox(sol["rhom"][j,z], 0, atol=TOL)
            push!(alpha, numbering["rhom"][j,z])
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        if !isapprox(sol["rhon_tildep"][j,n], 0, atol=TOL)
            push!(alpha, numbering["rhon_tildep"][j,n])
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        if !isapprox(sol["rhon_tildem"][j,n], 0, atol=TOL)
            push!(alpha, numbering["rhon_tildem"][j,n])
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        if !isapprox(sol["rhoz_tildep"][j,z], 0, atol=TOL)
            push!(alpha, numbering["rhoz_tildep"][j,z])
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        if !isapprox(sol["rhoz_tildem"][j,z], 0, atol=TOL)
            push!(alpha, numbering["rhoz_tildem"][j,z])
        end
    end
    for j=1:size(D, 1)
        if !isapprox(sol["phip"][j], 0, atol=TOL)
            push!(alpha, numbering["phip"][j])
        end
    end
    for j=1:size(D, 1)
        if !isapprox(sol["phim"][j], 0, atol=TOL)
            push!(alpha, numbering["phim"][j])
        end
    end
    for j=1:size(D, 1), l=1:length(L)
        if !isapprox(sol["psip"][j,l], 0, atol=TOL)
            push!(alpha, numbering["psip"][j,l])
        end
    end
    for j=1:size(D, 1), l=1:length(L)
        if !isapprox(sol["psim"][j,l], 0, atol=TOL)
            push!(alpha, numbering["psim"][j,l])
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        if !isapprox(sol["mu_tilde"][j,n], 0, atol=TOL)
            push!(alpha, numbering["mu_tilde"][j,n])
        end
    end
    for i=1:length(I), z=1:length(Z)
        if !isapprox(sol["nu_tilde"][i,z], 0, atol=TOL)
            push!(alpha, numbering["nu_tilde"][i,z])
        end
    end
    for j=1:size(D, 1), l=1:length(L)
        if !isapprox(sol["lambdap"][j,l], 0, atol=TOL)
            push!(alpha, numbering["lambdap"][j,l])
        end
    end
    for j=1:size(D, 1), l=1:length(L)
        if !isapprox(sol["lambdam"][j,l], 0, atol=TOL)
            push!(alpha, numbering["lambdam"][j,l])
        end
    end
    for n=1:length(N)
        if !isapprox(sol["Z_tilde"][n], 0, atol=TOL)
            push!(alpha, numbering["Z_tilde"][n])
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        if !isapprox(sol["gamma_tilde"][j,n], 0, atol=TOL)
            push!(alpha, numbering["gamma_tilde"][j,n])
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        if !isapprox(sol["z_tilde"][j,n], 0, atol=TOL)
            push!(alpha, numbering["z_tilde"][j,n])
        end
    end

    return alpha
end

# function to solve the LCP on the polytope of a given alpha
# if sol0 is provided, maximize the distance to this solution
function solve_alpha(alpha::Vector{Int}; sol0::Union{Nothing, Dict}=nothing)
    # create model
    m = Model(optimizer_with_attributes(Gurobi.Optimizer, "FeasibilityTol" => 1e-2, "OutputFlag" => 0))

    # define primal variables
    @variable(m, p[j=1:size(D, 1), z=1:length(Z)])
    @variable(m, r[1:size(D, 1), 1:length(N)])
    @variable(m, y[1:length(I), 1:size(D, 1), 1:length(Z)] >= 0) # power produced by tech i in slice j
    @variable(m, x[i=1:length(I), z=1:length(Z)] >= 0) # investment in technology i
    @variable(m, s[j=1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, y_tilde[1:size(D, 1), 1:length(N)] >= 0) # power produced in slice j
    @variable(m, x_tilde[1:length(I), 1:length(N)] >= 0) # investment in node n
    @variable(m, y_bar[1:length(MC_ex), 1:size(D, 1), 1:length(N)] >= 0) # power produced by existing tech i in slice j
    @variable(m, x_bar[1:length(MC_ex), 1:length(N)] >= 0) # part of existing technology i kept
    @variable(m, f_tilde[1:size(D, 1), l=1:length(L)]) # export from node n

    # define dual variables
    @variable(m, mu[1:length(I), 1:size(D, 1), 1:length(Z)] >= 0)
    @variable(m, mu_bar[1:length(MC_ex), 1:size(D, 1), 1:length(N)] >= 0)
    @variable(m, rho[j=1:size(D, 1), 1:length(Z)])
    @variable(m, rhon_tilde[1:size(D, 1), 1:length(N)])
    @variable(m, rhoz_tilde[1:size(D, 1), 1:length(Z)])
    @variable(m, phi[1:size(D, 1)])
    @variable(m, psi[1:size(D, 1), 1:length(L)])
    @variable(m, mu_tilde[1:size(D,1), 1:length(N)] >= 0)
    @variable(m, nu_tilde[1:length(I), 1:length(Z)] >= 0)
    @variable(m, lambdap[1:size(D, 1), 1:length(L)] >= 0)
    @variable(m, lambdam[1:size(D, 1), 1:length(L)] >= 0)
    @variable(m, delta[1:length(I)] >= 0)
    @variable(m, delta_bar[1:length(MC_ex), 1:length(N)] >= 0)

    @variable(m, Z_tilde[n=1:length(N)] >= 0)
    @variable(m, gamma_tilde[1:size(D, 1), 1:length(N)] >= 0)
    @variable(m, z_tilde[1:size(D, 1), 1:length(N)] >= 0)

    # define expressions
    @expression(m, Fy[i=1:length(I), j=1:size(D, 1), z=1:length(Z)], mu[i,j,z] + DT[j]*MC[i] - rho[j,z])
    @expression(m, Fy_bar[i=1:length(I), j=1:size(D, 1), n=1:length(N)], mu_bar[i,j,n] + DT[j]*MC_ex[i] - (rhop[j,Zn[n]] - rhom[j,Zn[n]]))
    @expression(m, Fx[i=1:length(I), z=1:length(Z)], I[i] + FC[i] - sum(mu[i,j,z] for j=1:size(D, 1)) + delta[i])
    @expression(m, Fx_bar[i=1:length(I), n=1:length(N)], FC_ex[i] - sum(mu_bar[i,j,n] for j=1:size(D, 1)) + delta_bar[i,n])
    @expression(m, Fs[j=1:size(D, 1), z=1:length(Z)], DT[j]*VOLL - rho[j,z])
    @expression(m, Fmu[i=1:length(I), j=1:size(D, 1), z=1:length(Z)], x[i,z] - y[i,j,z])
    @expression(m, Fmu_bar[i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)], x_bar[i,n] - y_bar[i,j,n])

    @expression(m, Fnu_tilde[i=1:length(I), z=1:length(Z)], x[i,z] - sum(x_tilde[i,n] for n in Nz[z]))

    @expression(m, Frho[j=1:size(D, 1), z=1:length(Z)], -(pp[j,z] - pm[j,z]) + sum(y[i,j,z] for i=1:length(I)) +
        sum(y_bar[i,j,n] for i=1:length(MC_ex), n=Nz[z]) - sum(D[j,n] for n in Nz[z]) + s[j,z])
    @expression(m, Frhoz_tilde[j=1:size(D, 1), z=1:length(Z)], -(pp[j,z] - pm[j,z]) + sum(z_tilde[j,n] for n in Nz[z]) +
        sum(y_tilde[j,n] for n in Nz[z]) - sum(D[j,n] for n in Nz[z]))
    @expression(m, Fy_tilde[j=1:size(D, 1), n=1:length(N)], mu_tilde[j,n] - rhoz_tilde[j,Zn[n]] + rhon_tilde[j,n])
    @expression(m, Fx_tilde[i=1:length(I), n=1:length(N)], nu_tilde[i, Zn[n]] - sum(mu_tilde[j,n] for j=1:size(D, 1)))
    @expression(m, Fp[j=1:size(D, 1), z=1:length(Z)], rho[j,z] + rhoz_tilde[j,z])
    @expression(m, Fr[j=1:size(D, 1), n=1:length(N)], -rhon_tilde[j,n] - phi[j] + sum(PTDF[n,l]*psi[j,l]
        for l=1:length(L)))
    @expression(m, Ff_tilde[j=1:size(D, 1), l=1:length(L)], -psi[j,l] + (lambdap[j,l] - lambdam[j,l]))
    @expression(m, Fmu_tilde[j=1:size(D, 1), n=1:length(N)], sum(x_bar[i,n] for i in 1:length(MC_ex)) + sum(x_tilde[i,n] for i=1:length(I)) - y_tilde[j,n])
    @expression(m, Frhon_tilde[j=1:size(D, 1), n=1:length(N)], (rp[j,n] - rm[j,n]) - z_tilde[j,n] - y_tilde[j,n] + D[j,n])
    @expression(m, Fphi[j=1:size(D, 1)], sum((rp[j,n] - rm[j,n]) for n=1:length(N)))
    @expression(m, Fpsi[j=1:size(D, 1), l=1:length(L)], f_tilde[j,l] - sum(PTDF[n,l]*(rp[j,n] - rm[j,n]) for n=1:length(N)))
    @expression(m, Flambdap[j=1:size(D,1), l=1:length(L)], TC[l] - f_tilde[j,l])
    @expression(m, Flambdam[j=1:size(D,1), l=1:length(L)], TC[l] + f_tilde[j,l])
    @expression(m, Fdelta[i=1:length(MC_ex), n=1:length(N)], X[i,n] - x_bar[i,n])

    @expression(m, FZ_tilde[n=1:length(N)], I_tilde - sum(gamma_tilde[j,n] for j=1:size(D, 1)))
    @expression(m, Fz_tilde[j=1:size(D, 1), n=1:length(N)], DT[j]*MC_tilde + gamma_tilde[j,n] + rhon_tilde[j,n] -
        rhoz_tilde[j,Zn[n]])
    @expression(m, Fgamma_tilde[j=1:size(D, 1), n=1:length(N)], Z_tilde[n] - z_tilde[j,n])

    # add constraints
    for j=1:size(D, 1), z=1:length(Z)
        @constraint(m, Fp[j,z] == 0)
    end
    for j in 1:size(D, 1), n in 1:length(N)
        @constraint(m, Fr[j,n] == 0)
    end
    for i=1:length(I), j=1:size(D, 1), z=1:length(Z)
        gidx = numbering["y"][i,j,z]
        if gidx ∈ alpha
            @constraint(m, Fy[i,j,z] == 0)
        else
            @constraint(m, Fy[i,j,z] >= 0)
            @constraint(m, y[i,j,z] == 0)
        end
    end
    for i=1:length(I), z=1:length(Z)
        gidx = numbering["x"][i,z]
        if gidx ∈ alpha
            @constraint(m, Fx[i,z] == 0)
        else
            @constraint(m, Fx[i,z] >= 0)
            @constraint(m, x[i,z] == 0)
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        gidx = numbering["s"][j,z]
        if gidx ∈ alpha
            @constraint(m, Fs[j,z] == 0)
        else
            @constraint(m, Fs[j,z] >= 0)
            @constraint(m, s[j,z] == 0)
        end
    end
    for i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)
        gidx = numbering["y_bar"][i,j,n]
        if gidx ∈ alpha
            @constraint(m, Fy_bar[i,j,n] == 0)
        else
            @constraint(m, Fy_bar[i,j,n] >= 0)
            @constraint(m, y_bar[i,j,n] == 0)
        end
    end
    for i=1:length(MC_ex), n=1:length(N)
        gidx = numbering["x_bar"][i,n]
        if gidx ∈ alpha
            @constraint(m, Fx_bar[i,n] == 0)
        else
            @constraint(m, Fx_bar[i,n] >= 0)
            @constraint(m, x_bar[i,n] == 0)
        end
    end
    for i=1:length(MC_ex), n=1:length(N)
        gidx = numbering["delta"][i,n]
        if gidx ∈ alpha
            @constraint(m, Fdelta[i,n] == 0)
        else
            @constraint(m, Fdelta[i,n] >= 0)
            @constraint(m, delta[i,n] == 0)
        end
    end
    for i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)
        gidx = numbering["mu_bar"][i,j,n]
        if gidx ∈ alpha
            @constraint(m, Fmu_bar[i,j,n] == 0)
        else
            @constraint(m, Fmu_bar[i,j,n] >= 0)
            @constraint(m, mu_bar[i,j,n] == 0)
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        gidx = numbering["y_tilde"][j,n]
        if gidx ∈ alpha
            @constraint(m, Fy_tilde[j,n] == 0)
        else
            @constraint(m, Fy_tilde[j,n] >= 0)
            @constraint(m, y_tilde[j,n] == 0)
        end
    end
    for i=1:length(I), n=1:length(N)
        gidx = numbering["x_tilde"][i,n]
        if gidx ∈ alpha
            @constraint(m, Fx_tilde[i,n] == 0)
        else
            @constraint(m, Fx_tilde[i,n] >= 0)
            @constraint(m, x_tilde[i,n] == 0)
        end
    end
    for j=1:size(D, 1), l=1:length(L)
        @constraint(m, Ff_tilde[j,l] == 0)
    end
    for i=1:length(I), j=1:size(D, 1), z=1:length(Z)
        gidx = numbering["mu"][i,j,z]
        if gidx ∈ alpha
            @constraint(m, Fmu[i,j,z] == 0)
        else
            @constraint(m, Fmu[i,j,z] >= 0)
            @constraint(m, mu[i,j,z] == 0)
        end
    end
    for j=1:size(D, 1), z=1:length(Z)
        @constraint(m, Frho[j,z] == 0)
    end
    for j=1:size(D, 1), n=1:length(N)
        @constraint(m, Frhon_tilde[j,n] == 0)
    end
    for j=1:size(D, 1), z=1:length(Z)
        @constraint(m, Frhoz_tilde[j,z] == 0)
    end
    for j=1:size(D, 1)
        @constraint(m, Fphi[j] == 0)
    end
    for j=1:size(D, 1), l=1:length(L)
        @constraint(m, Fpsi[j,l] == 0)
    end
    for j=1:size(D, 1), n=1:length(N)
        gidx = numbering["mu_tilde"][j,n]
        if gidx ∈ alpha
            @constraint(m, Fmu_tilde[j,n] == 0)
        else
            @constraint(m, Fmu_tilde[j,n] >= 0)
            @constraint(m, mu_tilde[j,n] == 0)
        end
    end
    for i=1:length(I), z=1:length(Z)
        gidx = numbering["nu_tilde"][i,z]
        if gidx ∈ alpha
            @constraint(m, Fnu_tilde[i,z] == 0)
        else
            @constraint(m, Fnu_tilde[i,z] >= 0)
            @constraint(m, nu_tilde[i,z] == 0)
        end
    end
    for j=1:size(D, 1), l=1:length(L)
        gidx = numbering["lambdap"][j,l]
        if gidx ∈ alpha
            @constraint(m, Flambdap[j,l] == 0)
        else
            @constraint(m, Flambdap[j,l] >= 0)
            @constraint(m, lambdap[j,l] == 0)
        end
    end
    for j=1:size(D, 1), l=1:length(L)
        gidx = numbering["lambdam"][j,l]
        if gidx ∈ alpha
            @constraint(m, Flambdam[j,l] == 0)
        else
            @constraint(m, Flambdam[j,l] >= 0)
            @constraint(m, lambdam[j,l] == 0)
        end
    end
    for n=1:length(N)
        gidx = numbering["Z_tilde"][n]
        if gidx ∈ alpha
            @constraint(m, FZ_tilde[n] == 0)
        else
            @constraint(m, FZ_tilde[n] >= 0)
            @constraint(m, Z_tilde[n] == 0)
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        gidx = numbering["gamma_tilde"][j,n]
        if gidx ∈ alpha
            @constraint(m, Fgamma_tilde[j,n] == 0)
        else
            @constraint(m, Fgamma_tilde[j,n] >= 0)
            @constraint(m, gamma_tilde[j,n] == 0)
        end
    end
    for j=1:size(D, 1), n=1:length(N)
        gidx = numbering["z_tilde"][j,n]
        if gidx ∈ alpha
            @constraint(m, Fz_tilde[j,n] == 0)
        else
            @constraint(m, Fz_tilde[j,n] >= 0)
            @constraint(m, z_tilde[j,n] == 0)
        end
    end

    # if sol0 is provided, add new variables for the L1 distance to this solution
    if !isnothing(sol0)

        bigM = 1e6
        # define the absolute values variables
        @variable(m, tp[j=1:size(D, 1), z=1:length(Z)])
        @variable(m, tr[1:size(D, 1), 1:length(N)])
        @variable(m, ty[1:length(I), 1:size(D, 1), 1:length(Z)])
        @variable(m, tx[i=1:length(I), z=1:length(Z)])
        @variable(m, ts[j=1:size(D, 1), 1:length(Z)])
        @variable(m, ty_bar[1:length(MC_ex), 1:size(D, 1), 1:length(N)])
        @variable(m, tx_bar[i=1:length(MC_ex), n=1:length(N)])
        @variable(m, tdelta[i=1:length(MC_ex), n=1:length(N)])
        @variable(m, tmu_bar[1:length(MC_ex), 1:size(D, 1), 1:length(N)])
        @variable(m, ty_tilde[1:size(D, 1), 1:length(N)])
        @variable(m, tx_tilde[1:length(I), 1:length(N)])
        @variable(m, tf_tilde[1:size(D, 1), l=1:length(L)])
        @variable(m, tmu[1:length(I), 1:size(D, 1), 1:length(Z)])
        @variable(m, trho[j=1:size(D, 1), 1:length(Z)])
        @variable(m, trhon_tilde[1:size(D, 1), 1:length(N)])
        @variable(m, trhoz_tilde[1:size(D, 1), 1:length(Z)])
        @variable(m, tphi[1:size(D, 1)])
        @variable(m, tpsi[1:size(D, 1), 1:length(L)])
        @variable(m, tmu_tilde[1:size(D,1), 1:length(N)])
        @variable(m, tnu_tilde[1:length(I), 1:length(Z)])
        @variable(m, tlambdap[1:size(D, 1), 1:length(L)])
        @variable(m, tlambdam[1:size(D, 1), 1:length(L)])
        @variable(m, tZ_tilde[n=1:length(N)])
        @variable(m, tgamma_tilde[1:size(D, 1), 1:length(N)])
        @variable(m, tz_tilde[1:size(D, 1), 1:length(N)])

        @variable(m, spp[j=1:size(D, 1), z=1:length(Z)] >= 0)
        @variable(m, spr[1:size(D, 1), 1:length(N)] >= 0)
        @variable(m, spy[1:length(I), 1:size(D, 1), 1:length(Z)] >= 0)
        @variable(m, spx[i=1:length(I), z=1:length(Z)] >= 0)
        @variable(m, sps[j=1:size(D, 1), 1:length(Z)] >= 0)
        @variable(m, spy_bar[1:length(MC_ex), 1:size(D, 1), 1:length(N)] >= 0)
        @variable(m, spx_bar[i=1:length(MC_ex), n=1:length(N)] >= 0)
        @variable(m, spdelta[i=1:length(MC_ex), n=1:length(N)] >= 0)
        @variable(m, spmu_bar[1:length(MC_ex), 1:size(D, 1), 1:length(N)] >= 0)
        @variable(m, spy_tilde[1:size(D, 1), 1:length(N)] >= 0)
        @variable(m, spx_tilde[1:length(I), 1:length(N)] >= 0)
        @variable(m, spf_tilde[1:size(D, 1), l=1:length(L)] >= 0)
        @variable(m, spmu[1:length(I), 1:size(D, 1), 1:length(Z)] >= 0)
        @variable(m, sprho[j=1:size(D, 1), 1:length(Z)] >= 0)
        @variable(m, sprhon_tilde[1:size(D, 1), 1:length(N)] >= 0)
        @variable(m, sprhoz_tilde[1:size(D, 1), 1:length(Z)] >= 0)
        @variable(m, spphi[1:size(D, 1)] >= 0)
        @variable(m, sppsi[1:size(D, 1), 1:length(L)] >= 0)
        @variable(m, spmu_tilde[1:size(D,1), 1:length(N)] >= 0)
        @variable(m, spnu_tilde[1:length(I), 1:length(Z)] >= 0)
        @variable(m, splambdap[1:size(D, 1), 1:length(L)] >= 0)
        @variable(m, splambdam[1:size(D, 1), 1:length(L)] >= 0)
        @variable(m, spZ_tilde[n=1:length(N)] >= 0)
        @variable(m, spgamma_tilde[1:size(D, 1), 1:length(N)] >= 0)
        @variable(m, spz_tilde[1:size(D, 1), 1:length(N)] >= 0)

        @variable(m, smp[j=1:size(D, 1), z=1:length(Z)] >= 0)
        @variable(m, smr[1:size(D, 1), 1:length(N)] >= 0)
        @variable(m, smy[1:length(I), 1:size(D, 1), 1:length(Z)] >= 0)
        @variable(m, smx[i=1:length(I), z=1:length(Z)] >= 0)
        @variable(m, sms[j=1:size(D, 1), 1:length(Z)] >= 0)
        @variable(m, smy_bar[1:length(MC_ex), 1:size(D, 1), 1:length(N)] >= 0)
        @variable(m, smx_bar[i=1:length(MC_ex), n=1:length(N)] >= 0)
        @variable(m, smdelta[i=1:length(MC_ex), n=1:length(N)] >= 0)
        @variable(m, smmu_bar[1:length(MC_ex), 1:size(D, 1), 1:length(N)] >= 0)
        @variable(m, smy_tilde[1:size(D, 1), 1:length(N)] >= 0)
        @variable(m, smx_tilde[1:length(I), 1:length(N)] >= 0)
        @variable(m, smf_tilde[1:size(D, 1), l=1:length(L)] >= 0)
        @variable(m, smmu[1:length(I), 1:size(D, 1), 1:length(Z)] >= 0)
        @variable(m, smrho[j=1:size(D, 1), 1:length(Z)] >= 0)
        @variable(m, smrhon_tilde[1:size(D, 1), 1:length(N)] >= 0)
        @variable(m, smrhoz_tilde[1:size(D, 1), 1:length(Z)] >= 0)
        @variable(m, smphi[1:size(D, 1)] >= 0)
        @variable(m, smpsi[1:size(D, 1), 1:length(L)] >= 0)
        @variable(m, smmu_tilde[1:size(D,1), 1:length(N)] >= 0)
        @variable(m, smnu_tilde[1:length(I), 1:length(Z)] >= 0)
        @variable(m, smlambdap[1:size(D, 1), 1:length(L)] >= 0)
        @variable(m, smlambdam[1:size(D, 1), 1:length(L)] >= 0)
        @variable(m, smZ_tilde[n=1:length(N)] >= 0)
        @variable(m, smgamma_tilde[1:size(D, 1), 1:length(N)] >= 0)
        @variable(m, smz_tilde[1:size(D, 1), 1:length(N)] >= 0)

        # binary variables
        @variable(m, dp[j=1:size(D, 1), z=1:length(Z)], Bin)
        @variable(m, dr[1:size(D, 1), 1:length(N)], Bin)
        @variable(m, dy[1:length(I), 1:size(D, 1), 1:length(Z)], Bin)
        @variable(m, dx[i=1:length(I), z=1:length(Z)], Bin)
        @variable(m, ds[j=1:size(D, 1), 1:length(Z)], Bin)
        @variable(m, dy_bar[1:length(MC_ex), 1:size(D, 1), 1:length(N)], Bin)
        @variable(m, dx_bar[i=1:length(MC_ex), n=1:length(N)], Bin)
        @variable(m, ddelta[i=1:length(MC_ex), n=1:length(N)], Bin)
        @variable(m, dmu_bar[1:length(MC_ex), 1:size(D, 1), 1:length(N)], Bin)
        @variable(m, dy_tilde[1:size(D, 1), 1:length(N)], Bin)
        @variable(m, dx_tilde[1:length(I), 1:length(N)], Bin)
        @variable(m, df_tilde[1:size(D, 1), l=1:length(L)], Bin)
        @variable(m, dmu[1:length(I), 1:size(D, 1), 1:length(Z)], Bin)
        @variable(m, drho[j=1:size(D, 1), 1:length(Z)], Bin)
        @variable(m, drhon_tilde[1:size(D, 1), 1:length(N)], Bin)
        @variable(m, drhoz_tilde[1:size(D, 1), 1:length(Z)], Bin)
        @variable(m, dphi[1:size(D, 1)], Bin)
        @variable(m, dpsi[1:size(D, 1), 1:length(L)], Bin)
        @variable(m, dmu_tilde[1:size(D,1), 1:length(N)], Bin)
        @variable(m, dnu_tilde[1:length(I), 1:length(Z)], Bin)
        @variable(m, dlambdap[1:size(D, 1), 1:length(L)], Bin)
        @variable(m, dlambdam[1:size(D, 1), 1:length(L)], Bin)
        @variable(m, dZ_tilde[n=1:length(N)], Bin)
        @variable(m, dgamma_tilde[1:size(D, 1), 1:length(N)], Bin)
        @variable(m, dz_tilde[1:size(D, 1), 1:length(N)], Bin)

        # define the abs constraints
        # define variable t as the difference between the new and old solutions
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], tp[j,z] == p[j,z] - sol0["p"][j,z])
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], tr[j,n] == r[j,n] - sol0["r"][j,n])
        @constraint(m, [i=1:length(I), j=1:size(D, 1), z=1:length(Z)], ty[i,j,z] == y[i,j,z] - sol0["y"][i,j,z])
        @constraint(m, [i=1:length(I), z=1:length(Z)], tx[i,z] == x[i,z] - sol0["x"][i,z])
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], ts[j,z] == s[j,z] - sol0["s"][j,z])
        @constraint(m, [i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)], ty_bar[i,j,n] == y_bar[i,j,n] - sol0["y_bar"][i,j,n])
        @constraint(m, [i=1:length(MC_ex), n=1:length(N)], tx_bar[i,n] == x_bar[i,n] - sol0["x_bar"][i,n])
        @constraint(m, [i=1:length(MC_ex), n=1:length(N)], tdelta[i,n] == delta[i,n] - sol0["delta"][i,n])
        @constraint(m, [i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)], tmu_bar[i,j,n] == mu_bar[i,j,n] - sol0["mu_bar"][i,j,n])
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], ty_tilde[j,n] == y_tilde[j,n] - sol0["y_tilde"][j,n])
        @constraint(m, [i=1:length(I), n=1:length(N)], tx_tilde[i,n] == x_tilde[i,n] - sol0["x_tilde"][i,n])
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], tf_tilde[j,l] == f_tilde[j,l] - sol0["f_tilde"][j,l])
        @constraint(m, [i=1:length(I), j=1:size(D, 1), z=1:length(Z)], tmu[i,j,z] == mu[i,j,z] - sol0["mu"][i,j,z])
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], trho[j,z] == rho[j,z] - sol0["rho"][j,z])
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], trhon_tilde[j,n] == rhon_tilde[j,n] - sol0["rhon_tilde"][j,n])
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], trhoz_tilde[j,z] == rhoz_tilde[j,z] - sol0["rhoz_tilde"][j,z])
        @constraint(m, [j=1:size(D, 1)], tphi[j] == phi[j] - sol0["phi"][j])
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], tpsi[j,l] == psi[j,l] - sol0["psi"][j,l])
        @constraint(m, [j=1:size(D,1), n=1:length(N)], tmu_tilde[j,n] == mu_tilde[j,n] - sol0["mu_tilde"][j,n])
        @constraint(m, [i=1:length(I), z=1:length(Z)], tnu_tilde[i,z] == nu_tilde[i,z] - sol0["nu_tilde"][i,z])
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], tlambdap[j,l] == lambdap[j,l] - sol0["lambdap"][j,l])
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], tlambdam[j,l] == lambdam[j,l] - sol0["lambdam"][j,l])
        @constraint(m, [n=1:length(N)], tZ_tilde[n] == Z_tilde[n] - sol0["Z_tilde"][n])
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], tgamma_tilde[j,n] == gamma_tilde[j,n] - sol0["gamma_tilde"][j,n])
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], tz_tilde[j,n] == z_tilde[j,n] - sol0["z_tilde"][j,n])

        # define t = sp - sm
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], tp[j,z] == spp[j,z] - smp[j,z])
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], tr[j,n] == spr[j,n] - smr[j,n])
        @constraint(m, [i=1:length(I), j=1:size(D, 1), z=1:length(Z)], ty[i,j,z] == spy[i,j,z] - smy[i,j,z])
        @constraint(m, [i=1:length(I), z=1:length(Z)], tx[i,z] == spx[i,z] - smx[i,z])
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], ts[j,z] == sps[j,z] - sms[j,z])
        @constraint(m, [i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)], ty_bar[i,j,n] == spy_bar[i,j,n] - smy_bar[i,j,n])
        @constraint(m, [i=1:length(MC_ex), n=1:length(N)], tx_bar[i,n] == spx_bar[i,n] - smx_bar[i,n])
        @constraint(m, [i=1:length(MC_ex), n=1:length(N)], tdelta[i,n] == spdelta[i,n] - smdelta[i,n])
        @constraint(m, [i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)], tmu_bar[i,j,n] == spmu_bar[i,j,n] - smmu_bar[i,j,n])
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], ty_tilde[j,n] == spy_tilde[j,n] - smy_tilde[j,n])
        @constraint(m, [i=1:length(I), n=1:length(N)], tx_tilde[i,n] == spx_tilde[i,n] - smx_tilde[i,n])
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], tf_tilde[j,l] == spf_tilde[j,l] - smf_tilde[j,l])
        @constraint(m, [i=1:length(I), j=1:size(D, 1), z=1:length(Z)], tmu[i,j,z] == spmu[i,j,z] - smmu[i,j,z])
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], trho[j,z] == sprho[j,z] - smrho[j,z])
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], trhon_tilde[j,n] == sprhon_tilde[j,n] - smrhon_tilde[j,n])
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], trhoz_tilde[j,z] == sprhoz_tilde[j,z] - smrhoz_tilde[j,z])
        @constraint(m, [j=1:size(D, 1)], tphi[j] == spphi[j] - smphi[j])
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], tpsi[j,l] == sppsi[j,l] - smpsi[j,l])
        @constraint(m, [j=1:size(D,1), n=1:length(N)], tmu_tilde[j,n] == spmu_tilde[j,n] - smmu_tilde[j,n])
        @constraint(m, [i=1:length(I), z=1:length(Z)], tnu_tilde[i,z] == spnu_tilde[i,z] - smnu_tilde[i,z])
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], tlambdap[j,l] == splambdap[j,l] - smlambdap[j,l])
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], tlambdam[j,l] == splambdam[j,l] - smlambdam[j,l])
        @constraint(m, [n=1:length(N)], tZ_tilde[n] == spZ_tilde[n] - smZ_tilde[n])
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], tgamma_tilde[j,n] == spgamma_tilde[j,n] - smgamma_tilde[j,n])
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], tz_tilde[j,n] == spz_tilde[j,n] - smz_tilde[j,n])

        # constraints on spp and smp with binary variables
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], spp[j,z] <= dp[j,z]*bigM)
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], spr[j,n] <= dr[j,n]*bigM)
        @constraint(m, [i=1:length(I), j=1:size(D, 1), z=1:length(Z)], spy[i,j,z] <= dy[i,j,z]*bigM)
        @constraint(m, [i=1:length(I), z=1:length(Z)], spx[i,z] <= dx[i,z]*bigM)
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], sps[j,z] <= ds[j,z]*bigM)
        @constraint(m, [i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)], spy_bar[i,j,n] <= dy_bar[i,j,n]*bigM)
        @constraint(m, [i=1:length(MC_ex), n=1:length(N)], spx_bar[i,n] <= dx_bar[i,n]*bigM)
        @constraint(m, [i=1:length(MC_ex), n=1:length(N)], spdelta[i,n] <= ddelta[i,n]*bigM)
        @constraint(m, [i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)], spmu_bar[i,j,n] <= dmu_bar[i,j,n]*bigM)
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], spy_tilde[j,n] <= dy_tilde[j,n]*bigM)
        @constraint(m, [i=1:length(I), n=1:length(N)], spx_tilde[i,n] <= dx_tilde[i,n]*bigM)
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], spf_tilde[j,l] <= df_tilde[j,l]*bigM)
        @constraint(m, [i=1:length(I), j=1:size(D, 1), z=1:length(Z)], spmu[i,j,z] <= dmu[i,j,z]*bigM)
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], sprho[j,z] <= drho[j,z]*bigM)
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], sprhon_tilde[j,n] <= drhon_tilde[j,n]*bigM)
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], sprhoz_tilde[j,z] <= drhoz_tilde[j,z]*bigM)
        @constraint(m, [j=1:size(D, 1)], spphi[j] <= dphi[j]*bigM)
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], sppsi[j,l] <= dpsi[j,l]*bigM)
        @constraint(m, [j=1:size(D,1), n=1:length(N)], spmu_tilde[j,n] <= dmu_tilde[j,n]*bigM)
        @constraint(m, [i=1:length(I), z=1:length(Z)], spnu_tilde[i,z] <= dnu_tilde[i,z]*bigM)
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], splambdap[j,l] <= dlambdap[j,l]*bigM)
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], splambdam[j,l] <= dlambdam[j,l]*bigM)
        @constraint(m, [n=1:length(N)], spZ_tilde[n] <= dZ_tilde[n]*bigM)
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], spgamma_tilde[j,n] <= dgamma_tilde[j,n]*bigM)
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], spz_tilde[j,n] <= dz_tilde[j,n]*bigM)

        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], smp[j,z] <= (1-dp[j,z])*bigM)
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], smr[j,n] <= (1-dr[j,n])*bigM)
        @constraint(m, [i=1:length(I), j=1:size(D, 1), z=1:length(Z)], smy[i,j,z] <= (1-dy[i,j,z])*bigM)
        @constraint(m, [i=1:length(I), z=1:length(Z)], smx[i,z] <= (1-dx[i,z])*bigM)
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], sms[j,z] <= (1-ds[j,z])*bigM)
        @constraint(m, [i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)], smy_bar[i,j,n] <= (1-dy_bar[i,j,n])*bigM)
        @constraint(m, [i=1:length(MC_ex), n=1:length(N)], smx_bar[i,n] <= (1-dx_bar[i,n])*bigM)
        @constraint(m, [i=1:length(MC_ex), n=1:length(N)], smdelta[i,n] <= (1-ddelta[i,n])*bigM)
        @constraint(m, [i=1:length(MC_ex), j=1:size(D, 1), n=1:length(N)], smmu_bar[i,j,n] <= (1-dmu_bar[i,j,n])*bigM)
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], smy_tilde[j,n] <= (1-dy_tilde[j,n])*bigM)
        @constraint(m, [i=1:length(I), n=1:length(N)], smx_tilde[i,n] <= (1-dx_tilde[i,n])*bigM)
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], smf_tilde[j,l] <= (1-df_tilde[j,l])*bigM)
        @constraint(m, [i=1:length(I), j=1:size(D, 1), z=1:length(Z)], smmu[i,j,z] <= (1-dmu[i,j,z])*bigM)
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], smrho[j,z] <= (1-drho[j,z])*bigM)
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], smrhon_tilde[j,n] <= (1-drhon_tilde[j,n])*bigM)
        @constraint(m, [j=1:size(D, 1), z=1:length(Z)], smrhoz_tilde[j,z] <= (1-drhoz_tilde[j,z])*bigM)
        @constraint(m, [j=1:size(D, 1)], smphi[j] <= (1-dphi[j])*bigM)
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], smpsi[j,l] <= (1-dpsi[j,l])*bigM)
        @constraint(m, [j=1:size(D,1), n=1:length(N)], smmu_tilde[j,n] <= (1-dmu_tilde[j,n])*bigM)
        @constraint(m, [i=1:length(I), z=1:length(Z)], smnu_tilde[i,z] <= (1-dnu_tilde[i,z])*bigM)
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], smlambdap[j,l] <= (1-dlambdap[j,l])*bigM)
        @constraint(m, [j=1:size(D, 1), l=1:length(L)], smlambdam[j,l] <= (1-dlambdam[j,l])*bigM)
        @constraint(m, [n=1:length(N)], smZ_tilde[n] <= (1-dZ_tilde[n])*bigM)
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], smgamma_tilde[j,n] <= (1-dgamma_tilde[j,n])*bigM)
        @constraint(m, [j=1:size(D, 1), n=1:length(N)], smz_tilde[j,n] <= (1-dz_tilde[j,n])*bigM)

        @expression(m, pos_sum, sum(sum(spp) + sum(spr) + sum(spx) + sum(spy) + sum(sps) + sum(spx_bar) + sum(spy_bar) + sum(spdelta) + sum(spmu_bar) +
            sum(spmu) + sum(spy_tilde) + sum(spx_tilde) + sum(spf_tilde) + sum(sprho) + sum(sprhon_tilde) +
            sum(sprhoz_tilde) + sum(spphi) + sum(sppsi) + sum(spmu_tilde) + sum(spnu_tilde) + sum(splambdap) +
            sum(splambdam) + sum(spZ_tilde) + sum(spgamma_tilde) + sum(spz_tilde)))

        @expression(m, neg_sum, sum(sum(smp) + sum(smr) + sum(smx) + sum(smy) + sum(sms) + sum(smx_bar) + sum(smy_bar) + sum(smdelta) + sum(smmu_bar) +
            sum(smmu) + sum(smy_tilde) + sum(smx_tilde) + sum(smf_tilde) + sum(smrho) + sum(smrhon_tilde) +
            sum(smrhoz_tilde) + sum(smphi) + sum(smpsi) + sum(smmu_tilde) + sum(smnu_tilde) + sum(smlambdap) +
            sum(smlambdam) + sum(smZ_tilde) + sum(smgamma_tilde) + sum(smz_tilde)))

        @constraint(m, pos_sum + neg_sum <= 1e6)

        @objective(m, Max, pos_sum + neg_sum)

    else
        @objective(m, Max, sum(DT[j]*(VOLL*s[j,z] + sum(MC[i]*y[i,j,z] for i=1:length(I))) for j=1:size(D, 1), z=1:length(Z)) +
            sum((I[i]+FC[i])*x[i,z] for i=1:length(I), z=1:length(Z)))
    end

    for var in [spp, spr, spy, spx, sps, spy_bar, spx_bar, spdelta, spmu_bar, spy_tilde, spx_tilde, spf_tilde, spmu, sprho,
        sprhon_tilde, sprhoz_tilde, spphi, sppsi, spmu_tilde, spnu_tilde, splambdap, splambdam, spZ_tilde, spgamma_tilde, spz_tilde]
        setupperbound.(var, 0)
    end
    for var in [smp, smr, smy, smx, sms, smy_bar, smx_bar, smdelta, smmu_bar, smy_tilde, smx_tilde, smf_tilde, smmu, smrho,
        smrhon_tilde, smrhoz_tilde, smphi, smpsi, smmu_tilde, smnu_tilde, smlambdap, smlambdam, smZ_tilde, smgamma_tilde, smz_tilde]
        setupperbound.(var, 0)
    end

    for var in [spp, spr, spy, spx, sps, spy_bar, spx_bar, spdelta, spmu_bar, spy_tilde, spx_tilde, spf_tilde, spmu, sprho,
        sprhon_tilde, sprhoz_tilde, spphi, sppsi, spmu_tilde, spnu_tilde, splambdap, splambdam, spZ_tilde, spgamma_tilde, spz_tilde]
        for v in var
            setupperbound(v, Inf)
            optimize!(m)
            if isapprox(objective_value(m), 1e6, atol=1e-2)
                println(v)
            end
            setupperbound(v, 0)
        end
    end
    for var in [smp, smr, smy, smx, sms, smy_bar, smx_bar, smdelta, smmu_bar, smy_tilde, smx_tilde, smf_tilde, smmu, smrho,
        smrhon_tilde, smrhoz_tilde, smphi, smpsi, smmu_tilde, smnu_tilde, smlambdap, smlambdam, smZ_tilde, smgamma_tilde, smz_tilde]
        for v in var
            setupperbound(v, Inf)
            optimize!(m)
            if isapprox(objective_value(m), 1e6, atol=1e-2)
                println(v)
            end
            setupperbound(v, 0)
        end
    end

    optimize!(m)

    @show objective_value(m)

    sol_final = Dict()
    push!(sol_final, "p" => value.(p))
    push!(sol_final, "r" => value.(r))
    push!(sol_final, "y" => value.(y))
    push!(sol_final, "x" => value.(x))
    push!(sol_final, "s" => value.(s))
    push!(sol_final, "y_bar" => value.(y_bar))
    push!(sol_final, "x_bar" => value.(x_bar))
    push!(sol_final, "delta" => value.(delta))
    push!(sol_final, "mu_bar" => value.(mu_bar))
    push!(sol_final, "y_tilde" => value.(y_tilde))
    push!(sol_final, "x_tilde" => value.(x_tilde))
    push!(sol_final, "f_tilde" => value.(f_tilde))
    push!(sol_final, "mu" => value.(mu))
    push!(sol_final, "rho" => value.(rho))
    push!(sol_final, "rhon_tilde" => value.(rhon_tilde))
    push!(sol_final, "rhoz_tilde" => value.(rhoz_tilde))
    push!(sol_final, "phi" => value.(phi))
    push!(sol_final, "psi" => value.(psi))
    push!(sol_final, "mu_tilde" => value.(mu_tilde))
    push!(sol_final, "nu_tilde" => value.(nu_tilde))
    push!(sol_final, "lambdap" => value.(lambdap))
    push!(sol_final, "lambdam" => value.(lambdam))
    push!(sol_final, "Z_tilde" => value.(Z_tilde))
    push!(sol_final, "gamma_tilde" => value.(gamma_tilde))
    push!(sol_final, "z_tilde" => value.(z_tilde))

    return termination_status(m), sol_final
end
