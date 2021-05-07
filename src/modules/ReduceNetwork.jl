## Reduce the network using nodes clustering

__precompile__()

module ReduceNetwork

## external modules

using JuMP, Gurobi, DataFrames, Clustering, Statistics, Ipopt
using CSV, DataFrames, LightGraphs, Tables
using Distances, LinearAlgebra

## expoted functions

export write_reduced_network, read_reduced_network, ReadSystemTables, indexsets

## define constants
const MVAbase = 100.0

## functions

## function to read system data from CSV tables

function ReadSystemTables(systemdir::AbstractString)

	# read static technical data using standard CSV.read function
	Generators = DataFrame(CSV.File(joinpath(systemdir, "Generators.csv")))
	GeneratorsRE = DataFrame(CSV.File(string(systemdir, "/GeneratorsRE.csv"),
		types=[String, String, String, Float64, Float64, String]))
	RenewableProfiles = DataFrame(CSV.File(string(systemdir, "/RenewableProfiles.csv"),
		types = vcat([String, String], repeat([Float64], 84))))
	Loads = DataFrame(CSV.File(string(systemdir, "/Loads.csv"),
		types=[String, String, String, Float64, Float64, Float64]))
	LoadProfiles = DataFrame(CSV.File(string(systemdir, "/LoadProfiles2018.csv"),
		types=[String, String, Float64, Float64, Float64, Float64, Float64, Float64]))
	Buses = DataFrame(CSV.File(string(systemdir, "/Buses.csv"),
		types=[String, String, Float64, Float64, String, Float64, Float64, Float64, Float64]))
	Lines = DataFrame(CSV.File(string(systemdir, "/Lines.csv"),
		types=[String, String, String, String, Float64, Float64, Float64, Float64]))
	Reserves = DataFrame(CSV.File(string(systemdir, "/Reserves.csv"),
		types=[String, Float64]))

	# check that static data makes sense
	if any(indexin(unique(Lines[!, :LineType]), ["AC", "DC"]) .== 0)
		error("unrecognized line type. Supported types are AC and DC.")
	end

	# return all read data frames
	return Generators, GeneratorsRE, RenewableProfiles, Loads, LoadProfiles,
		Buses, Lines, Reserves
end

## create index and sets

function indexsets(Generators::DataFrame, Loads::DataFrame, LoadProfiles::DataFrame, 
    Buses::DataFrame, Lines::DataFrame)
    
    Z = unique(Buses[!, :ZoneBus])
    G_Nidx = indexin(Generators[!, :BusGenerator], Buses[!, :Bus])
    Load_Busidx = indexin(Loads[!, :BusLoad], Buses[!, :Bus])
    Zn = indexin(Buses[!, :ZoneBus], Z)
    Line_FromBusidx = indexin(Lines[!, :FromBus], Buses[!, :Bus])
	Line_ToBusidx = indexin(Lines[!, :ToBus], Buses[!, :Bus])

    # generator subsets
	isfastG = Generators[!, :GeneratorType] .== "FAST"

    # elements grouped by nodes
    Gn = Dict{Int, Vector{Int}}()
    Line_Frombus = Dict{Int, Vector{Int}}()
    Line_Tobus = Dict{Int, Vector{Int}}()
    for n = 1:size(Buses, 1)
		Gn[n] = Int[]
	end
    for g = 1:size(Generators, 1)
        push!(Gn[G_Nidx[g]], g)
    end
    for b in 1:size(Buses, 1)
		Line_Frombus[b] = Int[]
		Line_Tobus[b] = Int[]
    end
    for l in 1:size(Lines, 1)
		push!(Line_Frombus[Line_FromBusidx[l]], l)
		push!(Line_Tobus[Line_ToBusidx[l]], l)
	end

    # elements grouped by zones
    Nz = Vector{Vector{Int}}(undef, length(Z))
    for z in 1:length(Z)
		Nz[z] = Int[]
    end
    for b in 1:size(Buses, 1)
		push!(Nz[Zn[b]], b)
    end
	
    return Z, Load_Busidx, Zn, Nz, G_Nidx, Gn, isfastG, Line_FromBusidx, Line_ToBusidx, Line_Frombus, 
        Line_Tobus
end

## function to build the PDF matrix of a network

function ptdf(L::DataFrame, N::DataFrame, L_FromNidx::Vector, L_ToNidx::Vector,
	L_Fromn::Dict, L_Ton::Dict)
	# we assume that node size(N, 1) is the hub node
	# compute matrix T
	T = Matrix{Float64}(undef, size(N, 1)-1, size(N, 1)-1)
	for m = 1:size(N, 1)-1
		for n = 1:size(N, 1)-1
			if m == n
				T[m, n] = 0
				if !isempty(L_Fromn[n])
					T[m,n] += sum(L[!, :Susceptance][l] for l=L_Fromn[n])
				end
				if !isempty(L_Ton[n])
					T[m,n] += sum(L[!, :Susceptance][l] for l=L_Ton[n])
				end
			elseif ((!isempty(L_Fromn[m])) && (n ∈ [L_ToNidx[l] for l=L_Fromn[m]])) || ((!isempty(L_Ton[m])) && (n ∈ [L_FromNidx[l] for l=L_Ton[m]]))
				ind = Vector{Int}()
				for l = L_Fromn[m]
					if n == L_ToNidx[l]
						push!(ind, l)
					end
				end
				for l = L_Ton[m]
					if n == L_FromNidx[l]
						push!(ind, l)
					end
				end
				T[m, n] = sum(-L[!, :Susceptance][i] for i=ind)
			else
				T[m, n] = 0
			end
		end
	end
	# compute M matrix
	M = zeros(size(L, 1), size(N, 1)-1)
	for l = 1:size(L, 1)
		if L_FromNidx[l] != size(N,1)
			M[l, L_FromNidx[l]] = L[!, :Susceptance][l]
		end
		if L_ToNidx[l] != size(N,1)
			M[l, L_ToNidx[l]] = -L[!, :Susceptance][l]
		end
	end
	Tinv = inv(T)
	return hcat(M*inv(T), zeros(size(L, 1)))
end


function write_reduced_network(sys_dir::AbstractString, out_dir::AbstractString,
    nb_clusters::Int, nb_pieces::Int)

    VOLL = 3000.0
    # read data from files
    Generators, GeneratorsRE, RenewableProfiles, Loads, LoadProfiles, Buses, 
        Lines, Reserves = ReadSystemTables(sys_dir)
    Z, Load_Busidx, Zn, Nz, G_Nidx, Gn, isfastG, Line_FromBusidx, Line_ToBusidx, Line_Frombus,
        Line_Tobus = indexsets(Generators, Loads, LoadProfiles, Buses, Lines)

    # gather all lines in a tuple Vector. To be used in the topology preserving
    # nodes_clustering
    lines_tuple = Vector{Tuple{Int, Int}}()
    for l=1:size(Lines, 1)
        push!(lines_tuple, (Line_FromBusidx[l], Line_ToBusidx[l]))
    end

    # nodal net demand
    hourly_D = zeros(size(LoadProfiles, 1), size(Buses, 1))
    for h in 1:size(LoadProfiles, 1)
        for l in 1:size(Loads, 1)
            hourly_D[h, Load_Busidx[l]] += Loads[l, :FractionLoadProfile]*
                LoadProfiles[h, Loads[l, :LoadProfile]]
        end
    end
    # nodal renewable production
    # hourly_R_init is used for nodes clustering 
    # hourly_R is the renewable scenario for 2030
    hourly_R_init = zeros(size(LoadProfiles, 1), size(Buses, 1))
    for h in 1:size(LoadProfiles, 1)
        for gre in 1:size(GeneratorsRE, 1)
            n = findfirst(Buses[!, :Bus] .== GeneratorsRE[gre, :BusGeneratorRE])
            hourly_R_init[h,n] += RenewableProfiles[h, GeneratorsRE[gre, :DynamicProfileRE]]*GeneratorsRE[gre, :FractionREProfile]
        end
    end
    # move negative renewable production to demand
    hourly_D[hourly_R_init .< 0] .= hourly_D[hourly_R_init .< 0] .- hourly_R_init[hourly_R_init .< 0]
    hourly_R_init[hourly_R_init .< 0] .= 0    
    hourly_R = zeros(size(LoadProfiles, 1), size(Buses, 1))
    # Renewable increase per country
    # sources: 
    # Austria
    # APG MASTERPLAN 2030
    # Belgium 
    # https://economie.fgov.be/sites/default/files/Files/Energy/Adequacy-and-flexibility-study-for-Belgium-2020-2030-Elia.pdf
    # Germany and Luxembourg
    # https://doi.org/10.1016/j.eneco.2020.104879
    # France 
    # https://www.ecologie.gouv.fr/sites/default/files/PPE-Executive%20summary.pdf
    # Netherlands 
    # https://www.tennet.eu/news/detail/roadmap-for-development-of-offshore-grid-in-2024-2030-period/
    # https://www.pv-magazine.com/2019/11/04/netherlands-to-reach-27-gw-of-solar-by-2030/
    solar_increase = Dict(
        "A" => 2.91, 
        "B" => 2.15, 
        "D" => 1.97, 
        "F" => 2.16, 
        "Lx" => 1.97, 
        "NL" => 3)
    wind_increase = Dict(
        "A" => 1.28, 
        "B" => 1.67, 
        "D" => 1.97, 
        "F" => 1.27, 
        "Lx" => 1.97, 
        "NL" => 1.55)
    for h in 1:size(LoadProfiles, 1)
        for gre in 1:size(GeneratorsRE, 1)
            n = findfirst(Buses[!, :Bus] .== GeneratorsRE[gre, :BusGeneratorRE])
            Buses[n, :Country]
            coef = 1
            if GeneratorsRE[gre, :FuelGeneratorRE] == "SOLAR"
                coef = solar_increase[Buses[n, :Country]]
            elseif GeneratorsRE[gre, :FuelGeneratorRE] == "WIND"
                coef = wind_increase[Buses[n, :Country]]
            end
            hourly_R[h,n] += coef*RenewableProfiles[h, GeneratorsRE[gre, :DynamicProfileRE]]*GeneratorsRE[gre, :FractionREProfile]
        end
    end
    hourly_R[hourly_R .< 0] .= 0

    # demand
    total_load = sum(hourly_D, dims=2)[:] - sum(hourly_R_init, dims=2)[:]
    periods_permutation = sortperm(total_load, rev=true)
    # build the first piece manually to have a high peak (100 highest hours)
    interval_lengths, approx_values = pwconstant(total_load, Int(8760/10), nb_pieces)

    breaks = collect(sum(interval_lengths[1:i]) for i in 1:nb_pieces)
    # DT = interval_lengths ./ size(LoadProfiles, 1)
    DT = interval_lengths ./ 1.0

    D = zeros(nb_pieces, size(Buses, 1)) # nodal demand
    for n in 1:size(Buses, 1)
        sorted_load = hourly_D[:,n][periods_permutation]
        d_aggregated = Vector{Float64}(undef, nb_pieces)
        d_aggregated[1] = mean(sorted_load[j] for j in 1:breaks[1])
        for i in 2:nb_pieces-1
            d_aggregated[i] = mean(sorted_load[j] for j in breaks[i-1]+1:breaks[i])
        end
        d_aggregated[nb_pieces] = mean(sorted_load[j] for
            j in breaks[nb_pieces]:size(LoadProfiles, 1))
        D[:,n] = d_aggregated
    end

    R_init = zeros(nb_pieces, size(Buses, 1)) # nodal res
    for n in 1:size(Buses, 1)
        r_sorted = hourly_R_init[:,n][periods_permutation]
        r_aggregated = Vector{Float64}(undef, nb_pieces)
        r_aggregated[1] = mean(r_sorted[j] for j in 1:breaks[1])
        for i in 2:nb_pieces-1
            r_aggregated[i] = mean(r_sorted[j] for j in breaks[i-1]+1:breaks[i])
        end
        r_aggregated[nb_pieces] = mean(r_sorted[j] for
            j in breaks[nb_pieces]:size(LoadProfiles, 1))
        R_init[:,n] = r_aggregated
    end
    
    R = zeros(nb_pieces, size(Buses, 1)) # nodal res
    for n in 1:size(Buses, 1)
        r_sorted = hourly_R[:,n][periods_permutation]
        r_aggregated = Vector{Float64}(undef, nb_pieces)
        r_aggregated[1] = mean(r_sorted[j] for j in 1:breaks[1])
        for i in 2:nb_pieces-1
            r_aggregated[i] = mean(r_sorted[j] for j in breaks[i-1]+1:breaks[i])
        end
        r_aggregated[nb_pieces] = mean(r_sorted[j] for
            j in breaks[nb_pieces]:size(LoadProfiles, 1))
        R[:,n] = r_aggregated
    end

    # network data
    N = Buses[!, :Bus]
    L = Lines[!, :Line]
    TC = Lines[!, :FlowLimitForw]

    # add susceptance to lines
    Lines[!, :Susceptance] = MVAbase./Lines[!, :Reactance]
    PTDF = ptdf(Lines, Buses, Line_FromBusidx, Line_ToBusidx, Line_Frombus,
        Line_Tobus)
    # regularize ptdfs: ignore coefficients smaller than 1e-4
    PTDF[abs.(PTDF) .<= 1e-4] .= 0

    # different types of costs
    fuels, _, _, MC = costs()

    # existing gens cap and costs
    X_bar = Generators[!, :Capacity]
    FC_ex = Generators[!, :FixedCost]
    MC_ex = Generators[!, :VariableCost]
    
    # reduce network
    # also return the initial network's LMP to be used in finding the reduced
    # capacities
    clusters, initial_LMPs = nodes_clustering(N, L, Zn, TC, DT, D, R_init, PTDF, G_Nidx, 
        Gn, X_bar, MC_ex, VOLL, lines_tuple, Line_FromBusidx, Line_ToBusidx, 
        Line_Frombus, Line_Tobus, nb_clusters)

    LMPs_reduced = zeros(size(D,1), maximum(clusters))
    for j in 1:size(D,1)
        for c in 1:maximum(clusters)
            LMPs_reduced[j,c] = mean(initial_LMPs[j,clusters .== c])
        end
    end

    L_red, N_red, TC_init, PTDF_red, D_red, R_red, G_Nidx_red, Gn_red = reduce_net(clusters,
        N, L, PTDF, Line_FromBusidx, Line_ToBusidx, D, R, G_Nidx, Gn, TC)
    PTDF_red[abs.(PTDF_red) .<= 1e-4] .= 0

    # Generate TC by minimizing LMP distance from initial net to reduce net
    TC_red = get_reduced_capacities(LMPs_reduced ./ DT, N_red, L_red, TC_init, 
        DT, D_red, R_red, PTDF_red, G_Nidx_red, Gn_red, X_bar, MC_ex, VOLL)

    # write tables to csv files
    CSV.write(joinpath(out_dir, "DT.csv"), Tables.table(Matrix(DT')))
    CSV.write(joinpath(out_dir, "D.csv"), Tables.table(D_red, header=N_red))
    CSV.write(joinpath(out_dir, "R.csv"), Tables.table(R_red, header=N_red))
    CSV.write(joinpath(out_dir, "G_Nidx_red.csv"), Tables.table(Matrix(G_Nidx_red')))
    CSV.write(joinpath(out_dir, "clusters.csv"), Tables.table(Matrix(clusters'), header=N))
    CSV.write(joinpath(out_dir, "TC.csv"), Tables.table(Matrix(TC_red'), header=L_red))
    CSV.write(joinpath(out_dir, "PTDF.csv"), Tables.table(Matrix(PTDF_red'), header=L_red))

    return
end

function read_reduced_network(in_dir::AbstractString, Z::Vector, Zn::Vector, Nz::Vector)
    N = names(DataFrame(CSV.File(joinpath(in_dir, "D.csv"))))
    L = names(DataFrame(CSV.File(joinpath(in_dir, "TC.csv"))))
    DT_df = DataFrame(CSV.File(joinpath(in_dir, "DT.csv")))
    DT = collect(DT_df[1, i] for i in 1:size(DT_df, 2))
    D_df = DataFrame(CSV.File(joinpath(in_dir, "D.csv")))
    D = collect(D_df[i,j] for i in 1:size(D_df, 1), j in 1:size(D_df, 2))
    R_df = DataFrame(CSV.File(joinpath(in_dir, "R.csv")))
    R = collect(R_df[i,j] for i in 1:size(R_df, 1), j in 1:size(R_df, 2))
    clusters_df = DataFrame(CSV.File(joinpath(in_dir, "clusters.csv")))
    clusters = collect(clusters_df[1, i] for i in 1:size(clusters_df, 2))
    TC_df = DataFrame(CSV.File(joinpath(in_dir, "TC.csv")))
    TC = collect(TC_df[1, i] for i in 1:size(TC_df, 2))
    G_Nidx_df = DataFrame(CSV.File(joinpath(in_dir, "G_Nidx_red.csv")))
    G_Nidx = collect(G_Nidx_df[1, i] for i in 1:size(G_Nidx_df, 2))
    PTDF_df = DataFrame(CSV.File(joinpath(in_dir, "PTDF.csv")))
    PTDF = collect(PTDF_df[i,j] for i in 1:size(PTDF_df, 1), j in 1:size(PTDF_df, 2))
    Gn = Dict{Int, Vector{Int}}()
    for c in 1:maximum(clusters)
        Gn[c] = Int[]
    end
    for g = 1:size(G_Nidx, 1)
        push!(Gn[G_Nidx[g]], g)
    end
    # update Zn and Nz
    Zn_red = Vector(undef, length(N))
    for i in 1:length(Zn)
        # check that all nodes of a cluster are in the same zone
        if isassigned(Zn_red, clusters[i])
            @assert Zn_red[clusters[i]] == Zn[i]
        end
        Zn_red[clusters[i]] = Zn[i]
    end
    Nz_red = Vector{Vector{Int}}(undef, length(Z))
    for z in 1:length(Z)
        Nz_red[z] = Int[]
    end
    for n in 1:length(N)
        push!(Nz_red[Zn_red[n]], n)
    end
    G_Zidx = collect(Zn_red[n] for n in G_Nidx)
    Gz = Dict{Int, Vector{Int}}()
    for c in 1:maximum(clusters)
        Gz[c] = Int[]
    end
    for g = 1:size(G_Zidx, 1)
        push!(Gz[G_Zidx[g]], g)
    end

    return N, L, DT, D, R, G_Nidx, Gn, G_Zidx, Gz, clusters, TC, PTDF, Zn_red, Nz_red
end

# define exogeneous costs parameters
function costs()
    # types of fuel
    fuels = ["Lignite", "Coal", "CCGT", "OCGT", "CCGT&CHP"]
    # costs source: https://doi.org/10.1016/j.eneco.2020.104879
    investment_cost = [285230, 202330, 80100, 56330, 94392] # €/MW
    fixed_cost = [101500, 46286, 16500, 9333, 16500]
    variable_cost_up = [42.12, 58, 70.62, 121.37, 50.95] # €/MWh
    variable_cost_down = [36.7, 44.5, 58.18, 93.42, 38.18] # €/MWh
    variable_cost = variable_cost_down + 0.5 .* (variable_cost_up .-
        variable_cost_down)

    return fuels, investment_cost, fixed_cost, variable_cost
end

# the goal of this function is to obtain a clustering of the nodes
# based on the LMPs obtained for each hour with the existing capacity
function nodes_clustering(N::Vector{String}, L::Vector{String}, Zn::Vector,
    TC::Vector{Float64}, DT::Vector{Float64}, D::T, R::T, PTDF::T, 
    G_Nidx::Vector, Gn::Dict, X_bar::Vector, MC_ex::Vector,
    VOLL::Float64, lines_tuple::Vector{Tuple{Int, Int}},
    L_FromBusidx::Vector, L_ToBusidx::Vector, 
    L_Frombus::Dict, L_Tobus::Dict,
    nb_clusters::Int) where T <: Matrix{Float64}

    LMPs = zeros(size(D, 1), length(N))
    lambdas = zeros(length(L))

    # build the nodal OPF model
    m = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0,
        "Threads" => 1))
    # define primal variables
    @variable(m, 0 <= y[g=1:length(X_bar)] <= X_bar[g])
    # produced by tech i in slice j
    @variable(m, yre[1:length(N)] >= 0)
    @variable(m, s[1:length(N)] >= 0)
    # define auxiliary variables
    @variable(m, f[l=1:length(L)]) # export
    # from node n
    @constraint(m, lambdap[l=1:length(L)], f[l] <= TC[l])
    @constraint(m, lambdam[l=1:length(L)], f[l] >= -TC[l])
    @variable(m, r[1:length(N)])

    # define primal constraints
    @constraint(m, rho[n=1:length(N)], -r[n] + sum(y[g] for
        g in Gn[n]) + yre[n] + s[n] == D[1,n])
    @constraint(m, sum(r[n] for n in 1:length(N)) == 0)
    @constraint(m, [l=1:length(L)], f[l] == sum(PTDF[l,n]*
        r[n] for n=1:length(N)))
    @constraint(m, [n=1:length(N)], yre[n] <= R[n])

    for j in 1:size(D, 1)
        @objective(m, Min, sum(DT[j]*(VOLL*s[n] + sum(MC_ex[g]*y[g] for g=Gn[n])) for n=1:length(N)))
        for n in 1:length(N)
            set_normalized_rhs(rho[n], D[j,n])
        end
        optimize!(m)
        lambdas = lambdas .+ dual.(lambdap) .+ dual.(lambdam)
        LMPs[j,:] = dual.(rho)
    end

    # build shortest path matrix 
    g = LightGraphs.SimpleGraph(length(N))
    for l = 1:length(L)
        if Zn[L_FromBusidx[l]] == Zn[L_ToBusidx[l]]
            LightGraphs.add_edge!(g, L_FromBusidx[l], L_ToBusidx[l])
        end
    end
    sp_mat = zeros(length(N), length(N))
    for n = 1:length(N)
        sp_mat[n,:] .= dijkstra_shortest_paths(g, n).dists
    end
    
    # build the weighted Laplacian matrix of the graph 
    Laplacian = zeros(length(N), length(N))
    max_val = maximum(euclidean(LMPs[:,L_FromBusidx[l]], LMPs[:,L_ToBusidx[l]]) for l in 1:length(L))
    for l in 1:length(L)
        Laplacian[L_FromBusidx[l], L_ToBusidx[l]] = (max_val - euclidean(LMPs[:,L_FromBusidx[l]], LMPs[:,L_ToBusidx[l]]))/max_val
        Laplacian[L_ToBusidx[l], L_FromBusidx[l]] = (max_val - euclidean(LMPs[:,L_FromBusidx[l]], LMPs[:,L_ToBusidx[l]]))/max_val
    end
    for n=1:length(N)
        Laplacian[n,n] = sum((max_val - euclidean(LMPs[:,L_FromBusidx[l]], LMPs[:,L_ToBusidx[l]]))/max_val
            for l in vcat(L_Frombus[n], L_Tobus[n]))
    end
    Laplacian = log.(1 .+ Laplacian)
    graph_vol = sum(diag(Laplacian)) # volume of the graph
    scaled_id_mat = sqrt(graph_vol)*Matrix(I, length(N), length(N))
    ectd_mat = pairwise(Mahalanobis(pinv(Laplacian)), scaled_id_mat, dims=2)
    dist_matrix = ectd_mat .+ 10*sp_mat
    result = hclust(dist_matrix, linkage=:average)
    clusters = cutree(result, k=nb_clusters)

    # rearrange clusters so that last node is in last cluster
    rearranged_clusters = copy(clusters)
    rearranged_clusters[clusters .== maximum(clusters)] .= clusters[end]
    rearranged_clusters[clusters .== clusters[end]] .= maximum(clusters)

    return rearranged_clusters, LMPs
end

function reduce_net(clusters::Vector{Int}, N::Vector{String}, L::Vector{String},
    PTDF::T, Line_FromBusidx::Vector, Line_ToBusidx::Vector, D::T, R::T, G_Nidx::Vector, 
    Gn::Dict, TC::Vector) where T <: Matrix{Float64}
    nz = maximum(clusters)

    Line_FromBusidx_red = Vector{Int}()
    Line_ToBusidx_red = Vector{Int}()
    L_red = Vector{String}()
    ind = 1
    for l in 1:length(L)
        if clusters[Line_FromBusidx[l]] != clusters[Line_ToBusidx[l]]
            push!(Line_FromBusidx_red, clusters[Line_FromBusidx[l]])
            push!(Line_ToBusidx_red, clusters[Line_ToBusidx[l]])
            push!(L_red, string("Line", ind))
            ind += 1
        end
    end
    N_red = collect(string("ReducedNode", i) for i in 1:nz)

    TC_red = zeros(length(L_red))
    Tf = zeros(length(L_red), length(L))
    for i in 1:length(L_red), j in 1:length(L)
        if (clusters[Line_FromBusidx[j]] == Line_FromBusidx_red[i]) &&
            (clusters[Line_ToBusidx[j]] == Line_ToBusidx_red[i])
            Tf[i,j] = 1
            TC_red[i] += TC[j]
        elseif (clusters[Line_ToBusidx[j]] == Line_FromBusidx_red[i]) &&
            (clusters[Line_FromBusidx[j]] == Line_ToBusidx_red[i])
            Tf[i,j] = -1
            TC_red[i] += TC[j]
        end
    end

    Tbz = zeros(nz, length(N))
    for i in 1:nz, j in 1:length(N)
        if clusters[j] == i
            Tbz[i,j] = 1
        end
    end
    @assert clusters[end] == maximum(clusters) # this assumption is reducing Tbz
    Tbz = Tbz[1:end-1, 1:end-1]

    # Hf is the slack-line reduced PTDF
    Hf = PTDF[:, 1:end-1]
    PTDF_red = Tf*Hf*Tbz'*inv(Tbz*Tbz')
    PTDF_red = hcat(PTDF_red, zeros(length(L_red)))

    # aggregate nodal demand and capacity according to new network
    D_red = zeros(size(D, 1), nz)
    for i in 1:size(D, 1), j in 1:size(D, 2)
        D_red[i, clusters[j]] += D[i,j]
    end
    R_red = zeros(size(R, 1), nz)
    for i in 1:size(R, 1), j in 1:size(R, 2)
        R_red[i, clusters[j]] += R[i,j]
    end
    
    # update indexsets quantities 
    G_Nidx_red = collect(clusters[G_Nidx[g]] for g in 1:length(G_Nidx))
    Gn_red = Dict{Int, Vector{Int}}()
    for c in 1:maximum(clusters)
        Gn_red[c] = Int[]
    end
    for g = 1:size(G_Nidx, 1)
        push!(Gn_red[G_Nidx_red[g]], g)
    end

    return L_red, N_red, TC_red, PTDF_red, D_red, R_red, G_Nidx_red, Gn_red
end

# return the thermal capacities of the reeduced network that minimize the distance
# between the initial and reduced LMPs
function get_reduced_capacities(LMPs::T, N::Vector{String}, L::Vector{String}, TC_init::Vector{Float64},
    DT::Vector{Float64}, D::T, R::T, PTDF::T, G_Nidx::Vector, Gn::Dict, 
    X_bar::Vector{Float64}, MC_ex::Vector{Float64}, VOLL::Float64) where T <: Matrix{Float64}

    TC_alljs = zeros(length(L), size(D, 1))
    allobjs = zeros(size(D,1))
        
    j = 10 # period on which to compute the reduction
    m = Model(optimizer_with_attributes(Gurobi.Optimizer, "NumericFocus" => 3,
        "OutputFlag" => 0))
    # define primal variables
    @variable(m, y[1:length(X_bar)] >= 0)
    @variable(m, yre[1:length(N)] >= 0)
    @variable(m, s[1:length(N)] >= 0)
    @variable(m, f[l=1:length(L)])
    @variable(m, r[1:length(N)])
    # define dual variables
    @variable(m, rho[n=1:length(N)])
    @variable(m, phi)
    @variable(m, psi[1:length(L)])
    @variable(m, mu[1:length(X_bar)] >= 0)
    @variable(m, mure[1:length(N)] >= 0)
    @variable(m, lambdap[1:length(L)] >= 0)
    @variable(m, lambdam[1:length(L)] >= 0)
    # define capacities variables
    @variable(m, 0 <= TC[l=1:length(L)] <= 3e4)

    # define primal constraints
    @constraint(m, rho_con[n=1:length(N)], -r[n] + sum(y[g] for
        g=Gn[n]) + yre[n] + s[n] - D[j,n] == 0)
    @constraint(m, phi_con, sum(r[n] for n in 1:length(N)) == 0)
    @constraint(m, psi_con[l=1:length(L)], f[l] - sum(PTDF[l,n]*
        r[n] for n=1:length(N)) == 0)
    @constraint(m, mu_con[g=1:length(X_bar)], X_bar[g] - y[g] >= 0)
    @constraint(m, mure_con[n=1:length(N)], R[j,n] - yre[n] >= 0)
    @constraint(m, lambdap_con[l=1:length(L)], TC[l] - f[l] >= 0)
    @constraint(m, lambdam_con[l=1:length(L)], f[l] + TC[l] >= 0)
    # define dual constraints
    @constraint(m, r_con[n=1:length(N)], -rho[n] + phi -
        sum(PTDF[l,n]*psi[l] for l=1:length(L)) == 0)
    @constraint(m, y_con[g=1:length(X_bar)], -DT[j]*MC_ex[g] +
        rho[G_Nidx[g]] - mu[g] <= 0)
    @constraint(m, yre_con[n=1:length(N)], rho[n] - mure[n] <= 0)
    @constraint(m, s_con[n=1:length(N)], -DT[j]*VOLL + rho[n] <= 0)
    @constraint(m, f_con[l=1:length(L)], psi[l] - lambdap[l] +
        lambdam[l] == 0)

    # strong duality
    @constraint(m, sum(DT[j]*VOLL*s[n] for n in 1:length(N)) + sum(DT[j]*MC_ex[g]*y[g] for
        g=1:length(X_bar)) <=
        sum(rho[n]*D[j,n] for n=1:length(N)) -
        sum(X_bar[g]*mu[g] for g=1:length(X_bar)) -
        sum(R[j,n]*mure[n] for n=1:length(N)) -
        sum(TC[l]*(lambdap[l] + lambdam[l]) for l=1:length(L)))

    # objective: minimize LMP distance
    @objective(m, Min, sum((LMPs[j,n]-rho[n]/DT[j])^2 for n=1:length(N)))

    current_sol = TC_init
    for l=1:length(L)
        JuMP.fix(TC[l], TC_init[l], force=true)
    end
    optimize!(m)
    status = termination_status(m)
    new_obj = objective_value(m)
    old_obj = 2*new_obj
    while (old_obj-new_obj)/old_obj >= 1e-6
        old_obj = new_obj
        # fix lambdas
        for l=1:length(L)
            JuMP.fix(lambdap[l], value(lambdap[l]), force=true)
            JuMP.fix(lambdam[l], value(lambdam[l]), force=true)
            JuMP.unfix(TC[l])
            set_lower_bound(TC[l], 0)
            set_upper_bound(TC[l], 3e4)
        end
        optimize!(m)
        status = termination_status(m)
        if primal_status(m) == MOI.FEASIBLE_POINT
            # fix TC
            for l=1:length(L)
                JuMP.unfix(lambdap[l])
                JuMP.unfix(lambdam[l])
                set_lower_bound(lambdap[l], 0)
                set_lower_bound(lambdam[l], 0)
                JuMP.fix(TC[l], value(TC[l]), force=true)
            end
            optimize!(m)
            status = termination_status(m)
            if primal_status(m) == MOI.FEASIBLE_POINT
                new_obj = objective_value(m)
                current_sol = value.(TC)
            else
                new_obj = old_obj
            end
        else
            new_obj = old_obj
        end
    end
    allobjs[j] = new_obj
    TC_alljs[:,j] = current_sol

    bestj = argmin(allobjs)
    return round.(current_sol)
end

# helper function used in the dynamic programming algorithm
function h(alpha::Int, beta::Int, f::Vector{Float64})
    return sum(f[t] for t=alpha+1:beta)
end
# recursive function to compute the best aggregation of demand in time
function F(s::Int, r::Int, f::Vector{Float64}, Ftab::Matrix, indtab::Matrix)
    if s == r
        return sum(f[j]^2 for j=1:s)
    elseif s == 1
        return h(0, r, f)^2/r
    else
        maxresult = findmax(collect((isassigned(Ftab, s-1, q) ? Ftab[s-1, q] :
            F(s-1, q, f, Ftab, indtab)) + h(q, r, f)^2/(r - q) for q in s-1:r-1))
        Ftab[s, r] = maxresult[1]
        # the term s - 1 in the next line comes from the fact that Konno1998
        # starts their counting at 0 and from the s-1:r-1 in the prev line
        indtab[s,r] = round(Int, s - 1 + maxresult[2])
        return Ftab[s, r]
    end
end

# determine the best piecewise contant functions of nb_pieces pieces to
# approximate the total_load curve
# This is an implementation of the algorithm proposed in https://doi.org/10.1016/0167-6377(88)90030-2
function pwconstant(load::Vector{Float64}, nbred::Int, nb_pieces::Int)
    # sort load
    loadperm = sortperm(load, rev=true)
    load_sorted = load[loadperm]

    # we first reduce the load by simple contant interval of nbred pieces to
    # reduce computational time (see: DOI 10.1109/EEM.2011.5953016)
    interval_length = floor(Int, length(load)/nbred)
    last_interval_length = length(load) - (nbred-1)*interval_length
    load_red = collect(s < nbred ? mean(load_sorted[s] for s in
        (s-1)*interval_length+1:s*interval_length) :
        mean(load_sorted[(s-1)*interval_length+1:end]) for s in 1:nbred)

    Ftab = Matrix(undef, nb_pieces, nbred-1)
    indtab = zeros(Int, nb_pieces, nbred-1)
    opt_val = sum(load_red[j]^2 for j=1:nbred-1) -
        F(nb_pieces, nbred-1, load_red, Ftab, indtab)
    breaks = Vector{Int}(undef, nb_pieces-1)
    q = indtab[nb_pieces, nbred-1]
    for s in nb_pieces-1:-1:1
        breaks[s] = q
        q = indtab[s, q]
    end
    approxval = Vector{Float64}(undef, nbred)
    value = mean(load_red[j] for j in 1:breaks[1])
    for j in 1:breaks[1]
        approxval[j] = value
    end
    for i in 2:nb_pieces-1
        value = mean(load_red[j] for j in breaks[i-1]+1:breaks[i])
        for j in breaks[i-1]+1:breaks[i]
            approxval[j] = value
        end
    end
    value = mean(load_red[j] for j in breaks[nb_pieces-1]+1:nbred)
    for j in breaks[nb_pieces-1]+1:nbred
        approxval[j] = value
    end

    # plot approximation
    load_red_long = Vector{Float64}()
    approxval_long = Vector{Float64}()
    for s in 1:length(load_red)
        if s < nbred
            for i in 1:interval_length
                push!(load_red_long, load_red[s])
                push!(approxval_long, approxval[s])
            end
        else
            for i in 1:last_interval_length
                push!(load_red_long, load_red[s])
                push!(approxval_long, approxval[s])
            end
        end
    end

    @assert length(breaks) == nb_pieces - 1
    interval_lengths = Vector{Int}()
    push!(interval_lengths, interval_length*breaks[1])
    for l in 2:nb_pieces-1
        push!(interval_lengths, interval_length*(breaks[l]-breaks[l-1]))
    end
    push!(interval_lengths, last_interval_length*(length(load_red)-breaks[nb_pieces-1]))
    approx_values = Vector{Float64}()
    for l in 1:nb_pieces-1
        push!(approx_values, approxval[breaks[l]])
    end
    push!(approx_values, approxval[length(load_red)])

    return interval_lengths, round.(approx_values)
end

end
