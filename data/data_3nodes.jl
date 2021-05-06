function load_data_3nodes()
    DT = [1760, 5500, 1500]./8760 # fraction of the time of vertical slice j
    MC = [25, 80, 6.5, 160] # marginal cost of each technology i
    I = [16, 5, 32, 2] # investment cost of each technology i
    FC = [0, 0, 0, 0]
    MC_ex = [25, 80, 6.5, 160] # marginal cost of each technology i
    FC_ex = [0, 0, 0, 0]
    MC_tilde = 0
    I_tilde = 200
    VOLL = 3000 # value of lost load
    N = [1, 2, 3] # nodes
    Z = [1, 2] # zones
    Rv = [0.0, 0.0]
    Nz = [[1,2], [3]]
    Zn = [1, 1, 2]
    D = [0 0 7086; 0 0 9004; 0 300 10869] # result_value of the deamnd in each slice j 
    R = zeros(size(D, 1), length(N))
    X = 10000*ones(length(I), length(N)) # limit on invested capacity
    X_bar = [600, 100] #existing capacity
    isfastG = [true, true]
    MC_ex = [80, 160]
    FC_ex = [0, 0]
    G_Nidx = [1, 2]
    Gn = Dict(1 => [1], 2 => [2], 3 => Int[])
    G_Zidx = [1, 1]
    Gz = Dict(1 => [1, 2], 2 => Int[])
    L = [1, 2, 3]
    TC = [200, 100, 50]
    PTDF = [1/3 1/3 2/3; 0 0 0; 2/3 -1/3 1/3] # PTDF_nl
        
    return G_Nidx, Gn, G_Zidx, Gz, isfastG, DT, MC, I, FC, MC_ex, FC_ex, MC_tilde, I_tilde, VOLL, D, 
        R, Rv, N, Z, Nz, Zn, X, X_bar, L, TC, PTDF
end