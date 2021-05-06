function plot_investment(sol_nodal::Dict, sol_fci::Dict, sol_fcr::Dict, 
    sol_fdi::Dict, sol_fdr::Dict, sol_pi::Dict, sol_pr::Dict)

    G = DataFrame(CSV.File("../data/CWE2018_daily/Generators.csv"))

    investment_costs = Vector{Float64}()
    op_costs = Vector{Float64}()
    total_costs = Vector{Float64}()
    investment_per_type = Matrix{Int}(undef, length(I), 4)
    decom_per_zone = zeros(length(Z), 4)
    investment_per_zone = zeros(length(Z), 4)
    investment_nr = Vector{Int}()
    decommissioning = Vector{Int}()
    decom_per_type = Matrix{Int}(undef, length(unique(G[!, :Technology])), 4)
    decom_per_zone_type = zeros(length(Z), length(unique(G[!, :Technology])), 4)
    investment_per_zone_type = zeros(length(Z), length(I)+1, 4)

    # NODAL RESULTS
    decom_long = X_bar .- sol_nodal["x_bar"]
    investment_cost = sum((I[i]+FC[i])*sol_nodal["x"][i,n] for i=1:length(I), n=1:length(N)) +
        sum(FC_ex[g]*sol_nodal["x_bar"][g] for g=1:length(X_bar))
    op_cost = sum(DT[j]*(VOLL*sol_nodal["s"][j,n] +
        sum(MC[i]*(sol_nodal["y"][i,j,n]) for i=1:length(I))) for j=1:size(D, 1), n=1:length(N)) +
        sum(DT[j]*MC_ex[g]*(sol_nodal["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1))
    total_cost = investment_cost + op_cost
    push!(investment_costs, investment_cost)
    push!(op_costs, op_cost)
    push!(total_costs, total_cost)
    investment_per_type[:, 1] = round.(Int, sum(sol_nodal["x"], dims=2) .* 1000)
    x_nodes = sum(sol_nodal["x"], dims=1)
    for n in 1:length(N)
        investment_per_zone[Zn[n], 1] += 1000*x_nodes[n]
    end
    push!(investment_nr, 0)
    push!(decommissioning, round.(Int, sum(X_bar.- sol_nodal["x_bar"]) .* 1000))
    decom_res = DataFrame(name=G[!, :Generator], tech=G[!, :Technology], decom=1000*decom_long)
    gd = groupby(decom_res, :tech)
    decom_per_type[:, 1] = round.(Int, sum.(map(x -> x[!, :decom], gd)))
    techs = first.(map(x -> x[!, :tech], gd))
    decom_per_node = collect(isempty(Gn[n]) ? 0 : sum(decom_long[g] for g in Gn[n]) for n in 1:length(N))
    for n in 1:length(N)
        decom_per_zone[Zn[n], 1] += 1000*decom_per_node[n]
    end
    for z in 1:length(Z), t in 1:length(techs)
        decom_per_zone_type[z, t, 1] = any(decom_res[g, :tech] == techs[t] for g in Gz[z]) ? sum(decom_long[g] for g in Gz[z] if decom_res[g, :tech] == techs[t]) : 0
    end
    for z in 1:length(Z)
        for i in 1:length(I)
            investment_per_zone_type[z, i, 1] = sum(sol_nodal["x"][i,n] for n in Nz[z])
        end
        investment_per_zone_type[z, length(I)+1, 1] = 0
    end

    # FBMC CENTRALIZED RESULTS

    decom_long = X_bar .- sol_fci["x_bar"]
    investment_cost = sum((I[i]+FC[i])*sol_fci["x"][i,z] for i=1:length(I), z=1:length(Z)) +
        sum(FC_ex[g]*sol_fci["x_bar"][g] for g=1:length(X_bar)) +
        I_tilde*sum(sol_fci["Z_tilde"][n] for n in 1:length(N))
    op_cost = sum(DT[j]*(VOLL*sol_fcr["s"][j,n] +
        sum(MC[i]*(sol_fcr["y"][i,j,n]) for i=1:length(I))) for j=1:size(D, 1), n=1:length(N)) +
        sum(DT[j]*MC_ex[g]*(sol_fcr["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1)) +
        sum(DT[j]*MC_tilde*sol_fcr["z_tilde"][j,n] for j=1:size(D, 1), n=1:length(N))
    total_cost = investment_cost + op_cost
    push!(investment_costs, investment_cost)
    push!(op_costs, op_cost)
    push!(total_costs, total_cost)
    investment_per_type[:, 2] = round.(Int, sum(sol_fci["x"], dims=2) .* 1000)
    x_nodes = sum(sol_fcr["x"], dims=1)
    for n in 1:length(N)
        investment_per_zone[Zn[n], 2] += 1000*x_nodes[n]
    end
    push!(investment_nr, round.(Int, sum(sol_fci["Z_tilde"]) .* 1000))
    push!(decommissioning, round.(Int, sum(X_bar.- sol_fci["x_bar"]) .* 1000))
    decom_res = DataFrame(name=G[!, :Generator], tech=G[!, :Technology], decom=1000*decom_long)
    gd = groupby(decom_res, :tech)
    decom_per_type[:, 2] = round.(Int, sum.(map(x -> x[!, :decom], gd)))
    decom_per_node = collect(isempty(Gn[n]) ? 0 : sum(X_bar[g] - sol_fci["x_bar"][g] for g in Gn[n]) for n in 1:length(N))
    for n in 1:length(N)
        decom_per_zone[Zn[n], 2] += 1000*decom_per_node[n]
    end
    for z in 1:length(Z), t in 1:length(techs)
        decom_per_zone_type[z, t, 2] = any(decom_res[g, :tech] == techs[t] for g in Gz[z]) ? sum(decom_long[g] for g in Gz[z] if decom_res[g, :tech] == techs[t]) : 0
    end
    for z in 1:length(Z)
        for i in 1:length(I)
            investment_per_zone_type[z, i, 2] = sum(sol_fcr["x"][i,n] for n in Nz[z])
        end
        investment_per_zone_type[z, length(I)+1, 2] = sum(sol_fci["Z_tilde"][n] for n in Nz[z])
    end
    z_ger = findfirst(Z .== "DE/LX")
    decom_DE_CCGT = collect((isempty(Gn[n]) || !any(G_Zidx[g] == z_ger && decom_res[g, :tech] == "CCGT" for g in Gn[n])) ? 0 : sum(decom_long[g] for g in Gn[n] if (G_Zidx[g] == z_ger && decom_res[g, :tech] == "CCGT")) for n in 1:length(N))
    investment_DE_CCGT = sol_fci["x"][1, :]


    # FBMC DECENTRALIZED RESULTS

    decom_long = X_bar .- sol_fdi["x_bar"]
    investment_cost = sum((I[i]+FC[i])*sol_fdi["x"][i,z] for i=1:length(I), z=1:length(Z)) +
        sum(FC_ex[g]*sol_fdi["x_bar"][g] for g=1:length(X_bar)) +
        I_tilde*sum(sol_fdi["Z_tilde"][n] for n in 1:length(N))
    op_cost = sum(DT[j]*(VOLL*sol_fdr["s"][j,n] +
        sum(MC[i]*(sol_fdr["y"][i,j,n]) for i=1:length(I))) for j=1:size(D, 1), n=1:length(N)) +
        sum(DT[j]*MC_ex[g]*(sol_fdr["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1)) +
        sum(DT[j]*MC_tilde*sol_fdr["z_tilde"][j,n] for j=1:size(D, 1), n=1:length(N))
    total_cost = investment_cost + op_cost
    push!(investment_costs, investment_cost)
    push!(op_costs, op_cost)
    push!(total_costs, total_cost)
    investment_per_type[:, 3] = round.(Int, sum(sol_fdi["x"], dims=2) .* 1000)
    x_nodes = sum(sol_fdr["x"], dims=1)
    for n in 1:length(N)
        investment_per_zone[Zn[n], 3] += 1000*x_nodes[n]
    end
    push!(investment_nr, round.(Int, sum(sol_fdi["Z_tilde"]) .* 1000))
    push!(decommissioning, round.(Int, sum(X_bar.- sol_fdi["x_bar"]) .* 1000))
    decom_res = DataFrame(name=G[!, :Generator], tech=G[!, :Technology], decom=1000*decom_long)
    gd = groupby(decom_res, :tech)
    decom_per_type[:, 3] = round.(Int, sum.(map(x -> x[!, :decom], gd)))
    decom_per_node = collect(isempty(Gn[n]) ? 0 : sum(X_bar[g] - sol_fdi["x_bar"][g] for g in Gn[n]) for n in 1:length(N))
    for n in 1:length(N)
        decom_per_zone[Zn[n], 3] += 1000*decom_per_node[n]
    end
    for z in 1:length(Z), t in 1:length(techs)
        decom_per_zone_type[z, t, 3] = any(decom_res[g, :tech] == techs[t] for g in Gz[z]) ? sum(decom_long[g] for g in Gz[z] if decom_res[g, :tech] == techs[t]) : 0
    end
    for z in 1:length(Z)
        for i in 1:length(I)
            investment_per_zone_type[z, i, 3] = sum(sol_fdr["x"][i,n] for n in Nz[z])
        end
        investment_per_zone_type[z, length(I)+1, 3] = sum(sol_fdi["Z_tilde"][n] for n in Nz[z])
    end

    ## ZONAL PA RESULTS

    decom_long = X_bar .- sol_pi["x_bar"]
    investment_cost = sum((I[i]+FC[i])*sol_pi["x"][i,z] for i=1:length(I), z=1:length(Z)) +
        sum(FC_ex[g]*sol_pi["x_bar"][g] for g=1:length(X_bar)) +
        I_tilde*sum(sol_pr["Z_tilde"][n] for n in 1:length(N))
    op_cost = sum(DT[j]*(VOLL*sol_pr["s"][j,n] +
        sum(MC[i]*(sol_pr["y"][i,j,n]) for i=1:length(I))) for j=1:size(D, 1), n=1:length(N)) +
        sum(DT[j]*MC_ex[g]*(sol_pr["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1)) +
        sum(DT[j]*MC_tilde*sol_pr["z_tilde"][j,n] for j=1:size(D, 1), n=1:length(N))
    total_cost = investment_cost + op_cost
    push!(investment_costs, investment_cost)
    push!(op_costs, op_cost)
    push!(total_costs, total_cost)
    investment_per_type[:, 4] = round.(Int, sum(sol_pi["x"], dims=2) .* 1000)
    x_nodes = sum(sol_pr["x"], dims=1)
    for n in 1:length(N)
        investment_per_zone[Zn[n], 4] += 1000*x_nodes[n]
    end
    push!(investment_nr, round.(Int, sum(sol_pr["Z_tilde"]) .* 1000))
    push!(decommissioning, round.(Int, sum(X_bar.- sol_pi["x_bar"]) .* 1000))
    decom_res = DataFrame(name=G[!, :Generator], tech=G[!, :Technology], decom=1000*decom_long)
    gd = groupby(decom_res, :tech)
    decom_per_type[:, 4] = round.(Int, sum.(map(x -> x[!, :decom], gd)))
    decom_per_node = collect(isempty(Gn[n]) ? 0 : sum(X_bar[g] - sol_pi["x_bar"][g] for g in Gn[n]) for n in 1:length(N))
    for n in 1:length(N)
        decom_per_zone[Zn[n], 4] += 1000*decom_per_node[n]
    end
    for z in 1:length(Z), t in 1:length(techs)
        decom_per_zone_type[z, t, 4] = any(decom_res[g, :tech] == techs[t] for g in Gz[z]) ? sum(decom_long[g] for g in Gz[z] if decom_res[g, :tech] == techs[t]) : 0
    end
    for z in 1:length(Z)
        for i in 1:length(I)
            investment_per_zone_type[z, i, 4] = sum(sol_pr["x"][i,n] for n in Nz[z])
        end
        investment_per_zone_type[z, length(I)+1, 4] = sum(sol_pr["Z_tilde"][n] for n in Nz[z])
    end

    # Display results
    df_res = DataFrame()
    df_res.Policy = ["Nodal", "FBMC-C", "FBMC-D", "PA-NR"]
    df_res.Decom = decommissioning
    df_res.NR = investment_nr
    df_res.OpCosts = op_costs
    df_res.InvCosts = investment_costs
    df_res.TotalCosts = total_costs
    println(df_res)

    # decomission per zone and type (3D plots)
    fig = plt.figure(figsize=[9.25, 7.64])
    agen = [221, 222, 223, 224]
    policies = ["Nodal", "FBMC-C", "FBMC-D", "PA-NR"]
    for i = 1:4
        ax = fig[:add_subplot](agen[i], projection="3d", title=policies[i])
        ax.set_xticks(collect(1:length(techs)))
        ax.set_xticklabels(techs, fontsize="xx-small", rotation="vertical")
        ax.set_yticks(10*collect(0:length(Z)-1))
        ax.set_yticklabels(Z)
        ax.set_xlabel("Technology", labelpad=10)
        ax.set_ylabel("Zone")
        ax.set_zlabel("Decommission [GW]")

        for z in 1:length(Z)
            xs = collect(1:length(techs))
            ys = decom_per_zone_type[z, :, i]
            ax.bar(xs, ys, zs=10*(z-1), zdir="y", alpha=0.8)
        end
    end
    PyPlot.savefig("../results/decom_per_zone_type.png")
    # investment per zone and type (3D plots)
    fig2 = plt.figure(figsize=[9.77, 8.3])
    agen = [221, 222, 223, 224]
    policies = ["Nodal", "FBMC-C", "FBMC-D", "PA-NR"]
    fuels = ["CCGT", "OCGT", "CCGT&CHP", "NR"]
    for i = 1:4
        ax = fig2[:add_subplot](agen[i], projection="3d", title=policies[i])
        ax.set_xticks(collect(1:length(I)+1))
        ax.set_xticklabels(fuels, fontsize="xx-small")
        ax.set_yticks(10*collect(0:length(Z)-1))
        ax.set_yticklabels(Z)
        ax.set_xlabel("Technology")
        ax.set_ylabel("Zone")
        ax.set_zlabel("Investment [GW]")

        for z in 1:length(Z)
            xs = collect(1:length(I)+1)
            ys = investment_per_zone_type[z, :, i]
            ax.bar(xs, ys, zs=10*(z-1), zdir="y", alpha=0.8)
        end
    end
    PyPlot.savefig("../results/investment_per_zone_type.png")
end

function plot_net(sol_nodal::Dict, sol_fc::Dict, sol_fd::Dict, sol_pa::Dict)

    # scales
    scale_euro = 1e-6 # put everything in m€
    scale_mw = 1e-3 # put everything in GW
    scale_cost = scale_euro/scale_mw

    # read N and L from files
    Generators, _, _, _, _, Buses, Lines, _ = ReadSystemTables("../data/CWE2018_daily/")

    L_FromBusidx = indexin(Lines[!, :FromBus], Buses[!, :Bus])
    L_ToBusidx = indexin(Lines[!, :ToBus], Buses[!, :Bus])

    # read clusters
    clusters = convert(Matrix, DataFrame(
        CSV.File("../data/reduced_net/clusters.csv")))
    clusters = reshape(clusters, size(clusters, 2), 1)

    # build nodal graph
    g = LightGraphs.SimpleGraph(length(clusters))
    for l = 1:size(Lines,1)
        if clusters[L_FromBusidx[l]] != clusters[L_ToBusidx[l]]
            LightGraphs.add_edge!(g, clusters[L_FromBusidx[l]], clusters[L_ToBusidx[l]])
        end
    end
    # make the node size proportional to the amount of COAL or Lignite on the node
    COAL_amount = zeros(size(Buses, 1))
    for b in 1:size(Buses, 1)
        if any((Generators.BusGenerator .== Buses[b, :Bus]) .& 
                            ((Generators.Technology .== "Lignite") .| (Generators.Technology .== "Coal")))
            sub_df = Generators[(Generators.BusGenerator .== Buses[b, :Bus]) .& 
                                ((Generators.Technology .== "Lignite") .| (Generators.Technology .== "Coal")), :]
            COAL_amount[b] = sum(sub_df[!, :Capacity])
        end
    end
    COAL_amount_red = collect(mapreduce(x -> clusters[x] == c ? COAL_amount[x] : 0, 
        +, 1:size(Buses, 1)) for c in 1:maximum(clusters))
    nodesize = ones(maximum(clusters))
    mcar = maximum(COAL_amount_red)
    for c in 1:maximum(clusters)
        nodesize[c] = 1 + 5*COAL_amount_red[c]/mcar
    end

    xcoord = Vector{Float64}(undef, maximum(clusters))
    ycoord = Vector{Float64}(undef, maximum(clusters))

    for c=1:maximum(clusters)
        xcoord[c] = 0.5*(maximum(Buses[i, :Xcoord] for i in 1:length(clusters) if clusters[i] == c) + 
            minimum(Buses[i, :Xcoord] for i in 1:length(clusters) if clusters[i] == c))
        ycoord[c] = 0.5*(maximum(Buses[i, :Ycoord] for i in 1:length(clusters) if clusters[i] == c) + 
            minimum(Buses[i, :Ycoord] for i in 1:length(clusters) if clusters[i] == c))
    end

    rho_nodal = sum(sol_nodal["rho"][j,:] for j in 1:size(D, 1))/8760 ./ scale_cost
    rho_fc = sum(sol_fc["rho"][j,:] for j in 1:size(D, 1))/8760 ./ scale_cost
    rho_fc = collect(rho_fc[Zn[n]] for n in 1:length(N))
    rho_fd = sum((sol_fd["rhop"][j,:] - sol_fd["rhom"][j,:]) for j in 1:size(D, 1))/8760 ./ scale_cost
    rho_fd = collect(rho_fd[Zn[n]] for n in 1:length(N))
    rho_pa = sum(sol_pa["rho"][j,:] for j in 1:size(D, 1))/8760 ./ scale_cost
    rho_pa = collect(rho_pa[Zn[n]] for n in 1:length(N))
    rho = [rho_nodal, rho_fc, rho_fd, rho_pa]

    # colors = colormap("RdBu", 100)
    colors = range(HSL(colorant"green"), stop=HSL(colorant"red"), length=100)

    minPrice = minimum(vcat(rho_nodal, rho_fc, rho_fd, rho_pa))
    maxPrice = maximum(vcat(rho_nodal, rho_fc, rho_fd, rho_pa))

    nodefillc = Vector{Vector}(undef, 4)
    for i in 1:4
        colorvec = Vector()
        for n=1:size(N,1)
            if minPrice == maxPrice
                ind = round(Int, length(colors)/2)
            else
                ind = round(Int, 100*(rho[i][n]-minPrice)/(maxPrice-minPrice))
            end
            if ind <= 0
                ind = 1
            elseif ind > 100
                ind = 100
            end
            push!(colorvec, colors[ind])
        end
        nodefillc[i] = colorvec
    end

    yticks = [1, 50, 100]
    labels = Dict(zip(yticks, string.(round.(range(maxPrice, minPrice, length=3), digits=2))))

    v = reshape(range(maxPrice, minPrice, length=100), 100, 1)
    hm = Gadfly.spy(reshape(repeat(v, 10), 100, 10), Guide.xlabel(""), Guide.ylabel(""), 
        Guide.xticks(label=false), Guide.yticks(ticks=yticks), Scale.y_continuous(labels = y -> labels[y]),
        Theme(key_position = :none, grid_line_width=0pt), 
        Scale.color_continuous(colormap=Scale.lab_gradient(colors...)))

    composition = compose(context(),
            (context(0.0, 0.1, 0.4, 0.4), gplot(g, xcoord, -ycoord, nodefillc=nodefillc[1], edgestrokec=colorant"darkgray", NODESIZE=3/sqrt(nv(g)), nodesize=nodesize)),
            (context(0.4, 0.1, 0.4, 0.4), gplot(g, xcoord, -ycoord, nodefillc=nodefillc[2], edgestrokec=colorant"darkgray", NODESIZE=3/sqrt(nv(g)), nodesize=nodesize)),
            (context(0.4, 0.6, 0.4, 0.4), gplot(g, xcoord, -ycoord, nodefillc=nodefillc[4], edgestrokec=colorant"darkgray", NODESIZE=3/sqrt(nv(g)), nodesize=nodesize)), 
            (context(0.0, 0.6, 0.4, 0.4), gplot(g, xcoord, -ycoord, nodefillc=nodefillc[3], edgestrokec=colorant"darkgray", NODESIZE=3/sqrt(nv(g)), nodesize=nodesize)),
            (context(), Compose.text(0.2, 0.1, "Nodal"), fill(colorant"dimgrey")),
            (context(), Compose.text(0.6, 0.1, "FBMC-C"), fill(colorant"dimgrey")),
            (context(), Compose.text(0.2, 0.6, "FBMC-D"), fill(colorant"dimgrey")),
            (context(), Compose.text(0.6, 0.6, "PA-NR"), fill(colorant"dimgrey")),
            (context(), Compose.text(0.85, 0.3, "[€/MWh]"), fill(colorant"dimgrey")),
            (context(0.7, 0.3, 0.4, 0.5), render(hm)))

    Compose.draw(SVGJS("../results/prices_allpolicies_COALandlignite.svg"), composition)
end

function plot_welfare(sol_nodal::Dict, sol_fci::Dict, sol_fcr::Dict, 
    sol_fdi::Dict, sol_fdr::Dict, sol_pi::Dict, sol_pr::Dict)

    _, _, _, _, _, Buses, _, _ = ReadSystemTables("../data/CWE2018_daily/")

    countries = unique(Buses[!, :Country])
    # read clusters
    clusters = convert(Matrix, DataFrame(CSV.File("../data/reduced_net/clusters.csv")))
    clusters = reshape(clusters, size(clusters, 2), 1)
    # build a map from N_red to countries
    Nc = Dict()
    for c in countries
        push!(Nc, c => unique(clusters[findall(Buses[!, :Country] .== c)]))
    end
    Zc = Dict()
    for c in countries
        b = findfirst(Buses[!, :Country] .== c)
        z = findfirst(Buses[b, :ZoneBus] .== Z)
        push!(Zc, c => z)
    end

    # compute the welfare allocation for each type of agent
    sol_nodal["rho"] .= sol_nodal["rho"] ./ DT
    sol_fci["rho"] .= sol_fci["rho"] ./ DT
    sol_fdi["rhop"] .= sol_fdi["rhop"] ./ DT
    sol_fdi["rhom"] .= sol_fdi["rhom"] ./ DT
    sol_pi["rho"] .= sol_pi["rho"] ./ DT

    # Nodal 
    welfare_nodal = Vector{Float64}()
    push!(welfare_nodal, 
        sum(DT[j]*(sol_nodal["rho"][j,n] - MC[i])*(sol_nodal["y"][i,j,n]) for i=1:length(I), j=1:size(D, 1), n=1:length(N)) +
        sum(DT[j]*(sol_nodal["rho"][j,G_Nidx[g]] - MC_ex[g])*(sol_nodal["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1)))
    push!(welfare_nodal, sum(DT[j]*((VOLL - sol_nodal["rho"][j,n])*(D[j,n] - sol_nodal["s"][j,n])) for j in 1:size(D, 1), n in 1:length(N)))
    push!(welfare_nodal, -sum(DT[j]*sol_nodal["rho"][j,n]*sol_nodal["r"][j,n] for j=1:size(D, 1), n=1:length(N)))
    push!(welfare_nodal, 0)

    # FBMC-C 
    welfare_fc = Vector{Float64}()
    push!(welfare_fc, 
        sum(DT[j]*(sol_fci["rho"][j,z] - MC[i])*(sol_fci["y"][i,j,z]) for i=1:length(I), j=1:size(D, 1), z=1:length(Z)) +
        sum(DT[j]*(sol_fci["rho"][j,G_Zidx[g]] - MC_ex[g])*(sol_fci["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1)))
    push!(welfare_fc, sum(DT[j]*((VOLL - sol_fci["rho"][j,z])*(sum(D[j,n] for n in Nz[z]) - sol_fci["s"][j,z])) for j in 1:size(D, 1), z in 1:length(Z)))
    push!(welfare_fc, -sum(DT[j]*sol_fci["rho"][j,z]*sol_fci["p"][j,z] for j=1:size(D, 1), z=1:length(Z)))
    push!(welfare_fc, sum(DT[j]*MC[i]*sol_fci["y"][i,j,z] for i=1:length(I), j=1:size(D, 1), z=1:length(Z)) +
        sum(DT[j]*MC_ex[g]*(sol_fci["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1)) - 
        (sum(DT[j]*MC[i]*sol_fcr["y"][i,j,n] for i=1:length(I), j=1:size(D, 1), n=1:length(N)) +
            sum(DT[j]*MC_ex[g]*(sol_fcr["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1)))
    )

    # FBMC-D
    welfare_fd = Vector{Float64}()
    push!(welfare_fd, 
        sum(DT[j]*(sol_fdi["rhop"][j,z] - sol_fdi["rhom"][j,z] - MC[i])*(sol_fdi["y"][i,j,z]) for i=1:length(I), j=1:size(D, 1), z=1:length(Z)) +
        sum(DT[j]*(sol_fdi["rhop"][j,G_Zidx[g]] - sol_fdi["rhom"][j,G_Zidx[g]] - MC_ex[g])*(sol_fdi["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1)))
    push!(welfare_fd, sum(DT[j]*((VOLL - sol_fdi["rhop"][j,z] + sol_fdi["rhom"][j,z])*(sum(D[j,n] for n in Nz[z]) - sol_fdi["s"][j,z])) for j in 1:size(D, 1), z in 1:length(Z)))
    push!(welfare_fd, -sum(DT[j]*(sol_fdi["rhop"][j,z] - sol_fdi["rhom"][j,z])*(sol_fdi["pp"][j,z] - sol_fdi["pm"][j,z]) for j=1:size(D, 1), z=1:length(Z)))
    push!(welfare_fd, sum(DT[j]*MC[i]*sol_fdi["y"][i,j,z] for i=1:length(I), j=1:size(D, 1), z=1:length(Z)) +
        sum(DT[j]*MC_ex[g]*(sol_fdi["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1)) - 
        (sum(DT[j]*MC[i]*sol_fdr["y"][i,j,n] for i=1:length(I), j=1:size(D, 1), n=1:length(N)) +
            sum(DT[j]*MC_ex[g]*(sol_fdr["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1)))
    )

    # PA-SR
    welfare_pa = Vector{Float64}()
    push!(welfare_pa, 
        sum(DT[j]*(sol_pi["rho"][j,z] - MC[i])*(sol_pi["y"][i,j,z]) for i=1:length(I), j=1:size(D, 1), z=1:length(Z)) +
        sum(DT[j]*(sol_pi["rho"][j,G_Zidx[g]] - MC_ex[g])*(sol_pi["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1)))
    push!(welfare_pa, sum(DT[j]*((VOLL - sol_pi["rho"][j,z])*(sum(D[j,n] for n in Nz[z]) - sol_pi["s"][j,z])) for j in 1:size(D, 1), z in 1:length(Z)))
    push!(welfare_pa, -sum(DT[j]*sol_pi["rho"][j,z]*sol_pi["p"][j,z] for j=1:size(D, 1), z=1:length(Z)))
    push!(welfare_pa, sum(DT[j]*MC[i]*sol_pi["y"][i,j,z] for i=1:length(I), j=1:size(D, 1), z=1:length(Z)) +
        sum(DT[j]*MC_ex[g]*(sol_pi["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1)) - 
        (sum(DT[j]*MC[i]*sol_pr["y"][i,j,n] for i=1:length(I), j=1:size(D, 1), n=1:length(N)) +
            sum(DT[j]*MC_ex[g]*(sol_pr["y_bar"][g,j]) for g=1:length(X_bar), j=1:size(D, 1)))
    )

    # produce the load-average price for each country
    average_price_nodal = Vector{Float64}(undef, length(countries))
    for (i, c) in enumerate(countries)
        average_price_nodal[i] = sum(sol_nodal["rho"][j,n]*D[j,n] for j in 1:size(D, 1), n=Nc[c])/
            sum(D[j,n] for j in 1:size(D, 1), n=Nc[c])
    end
    average_price_fc = Vector{Float64}(undef, length(countries))
    for (i, c) in enumerate(countries)
        average_price_fc[i] = sum(sol_fci["rho"][j,Zc[c]]*sum(D[j,n] for n in Nc[c]) for j in 1:size(D, 1))/
            sum(sum(D[j,n] for n in Nc[c]) for j in 1:size(D, 1))
    end
    average_price_fd = Vector{Float64}(undef, length(countries))
    for (i, c) in enumerate(countries)
        average_price_fd[i] = sum((sol_fdi["rhop"][j,Zc[c]] - sol_fdi["rhom"][j,Zc[c]])*sum(D[j,n] for n in Nc[c]) for j in 1:size(D, 1))/
            sum(sum(D[j,n] for n in Nc[c]) for j in 1:size(D, 1))
    end
    average_price_pa = Vector{Float64}(undef, length(countries))
    for (i, c) in enumerate(countries)
        average_price_pa[i] = sum(sol_pi["rho"][j,Zc[c]]*sum(D[j,n] for n in Nc[c]) for j in 1:size(D, 1))/
            sum(sum(D[j,n] for n in Nc[c]) for j in 1:size(D, 1))
    end

    colors = ["rgba(93, 164, 214, 0.85)", "rgba(255, 144, 14, 0.85)", "rgba(44, 160, 101, 0.85)", "rgba(214, 39, 40, 0.85)"]

    # Load average price
    names_agents = ["Producers", "Consumers", "TSO-congestion", "TSO-redispatch", "Total"]

    trace2 = PlotlyJS.bar(;x=names_agents, y=welfare_fc .- welfare_nodal, 
        name="FBMC-C", marker_color=colors[1])
    trace3 = PlotlyJS.bar(;x=names_agents, y=welfare_fd .- welfare_nodal, 
        name="FBMC-D", marker_color=colors[2])
    trace4 = PlotlyJS.bar(;x=names_agents, y=welfare_pa .- welfare_nodal, 
        name="PA", marker_color=colors[3])
    layout = Layout(;barmode="group",
        yaxis_title="Surplus [M€/yr]", yaxis_automargin=true, font_size=15)
    data = [trace2, trace3, trace4]
    p = PlotlyJS.plot(data, layout)
    PlotlyJS.savefig(p, "../results/welfare_reallocation.eps", scale=5)
end
