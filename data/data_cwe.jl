function costs()
    # types of fuel
    fuels = ["CCGT", "OCGT", "CCGT&CHP"]
    # costs source: https://doi.org/10.1016/j.eneco.2020.104879
    investment_cost = [80100, 56330, 94392] # €/MW
    fixed_cost = [16500, 9333, 16500]
    variable_cost_up = [70.62, 121.37, 50.95] # €/MWh
    variable_cost_down = [58.18, 93.42, 38.18] # €/MWh
    variable_cost = variable_cost_down + 0.25 .* (variable_cost_up .-
        variable_cost_down)

    return fuels, investment_cost, fixed_cost, variable_cost
end

function load_data_cwe()

    # exogeneous parameters
    VOLL = 3000.0
    # define turbojet as a slightly more expansive oc gaz power plant
    MC_tilde = 125
    I_tilde = 100000

    # read data from files
    Generators, GeneratorsRE, RenewableProfiles, Loads, LoadProfiles, Buses, 
        Lines, Reserves = ReadSystemTables("../data/CWE2018_daily/")
    Z, Load_Busidx, Zn, Nz, G_Nidx, Gn, isfastG, Line_FromBusidx, Line_ToBusidx, Line_Frombus,
        Line_Tobus = indexsets(Generators, Loads, LoadProfiles, Buses, Lines)

    # network data
    N = Buses[!, :Bus]
    L = Lines[!, :Line]

    # different types of costs
    fuels, I, FC, MC = costs()

    # existing gens cap and costs
    X_bar = Generators[!, :Capacity]
    FC_ex = Generators[!, :FixedCost]
    MC_ex = Generators[!, :VariableCost]

    # get reduced network data
    N, L, DT, D, R, G_Nidx, Gn, G_Zidx, Gz, clusters, TC, PTDF, Zn, Nz =
        read_reduced_network("../data/reduced_net/", Z, Zn, Nz)
    
    # get zonal reserve vector
    Rv = Vector{Float64}()
    for z in Z 
        ind = findfirst(Reserves[!, :Zone] .== z)
        push!(Rv, Reserves[ind, :ReserveValue])
    end

    # SCENARIO
    # source: https://eepublicdownloads.entsoe.eu/clean-documents/tyndp-documents/TYNDP%202014/141017_SOAF%202014-2030.pdf
    # decommissioning of Nuclear in Belgium and Germany and some in France
    # France source: https://cnpp.iaea.org/countryprofiles/France/France.htm 
    # EDF released the site where pairs of reactors will be shut down. The list counts FESSENHEIM + 7 sites, which 
    # yields 16 potential reactor closure whereas only 14 will be closed. We took the two closest sites among the 
    # 7 (Bugey and Tricastin, on the Rhone) and assume each will loose one reactor.  
    decom_gens = ["CHINON 1", "CHINON 2", "FESSENHEIM 1", "FESSENHEIM 2", "BLAYAIS 1", "BLAYAIS 2", "BUGEY 2", "CRUAS 1", "CRUAS 2", 
        "DAMPIERRE 1", "DAMPIERRE 2", "GRAVELINES 1", "GRAVELINES 2", "TRICASTIN 1"]
    for g in 1:size(Generators, 1)
        if Generators[g, :Generator] in decom_gens
            X_bar[g] = 0
        end
        if Generators[g, :Technology] == "Nuclear" && Generators[g, :Country] in ["BE", "DE"]
            X_bar[g] = 0
        end
    end
    # max investment in some technologies
    X = [150e3, 150e3, 15e3]

    scale_euro = 1e-6 # put everything in m€
    scale_mw = 1e-3 # put everything in GW
    scale_cost = scale_euro/scale_mw
    # costs are in m€/GW(h) = k€/MW(h)

    return G_Nidx, Gn, G_Zidx, Gz, isfastG, DT, MC*scale_cost, I*scale_cost, FC*scale_cost, MC_ex*scale_cost,
        FC_ex*scale_cost, MC_tilde*scale_cost, I_tilde*scale_cost,
        VOLL*scale_cost, D*scale_mw, R*scale_mw, Rv*scale_mw, N, Z, Nz, Zn,
        X*scale_mw, X_bar*scale_mw, L, TC*scale_mw, PTDF
end
