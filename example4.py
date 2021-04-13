import chemgaloo as cg

chem1 = cg.chemical(2.)
chem2 = cg.chemical(2.)
chem3 = cg.chemical()
chem4 = cg.chemical()

cg.reactor.__cstr_c_initializer__(chemicals = [chem1, chem2])

# rxn1: chem1 + chem2 -> chem3
stois1 = [1, 1, 1]
k1 = 1.0
rxn1 = cg.reaction(
    rxtns = [chem1, chem2], 
    prdts = [chem3],
    stois = stois1,
    k = k1
    )

# rxn2: chem1 + chem3 -> chem4
stois2 = [1, 1, 1]
k2 = 1.0
rxn2 = cg.reaction(
    rxtns = [chem1, chem3],
    prdts = [chem4],
    stois = stois2,
    k = k2
)

cg.reactor.cstr(
    chemicals = [chem1, chem2, chem3, chem4],
    reactions = [rxn1, rxn2],
    n_thread = 1,
    res_time_dist_f = 'Gaussian',
    rtf_param1 = 3.,
    rtf_param2 = 1e-5,
    dt = 0.01,
    conv_thr = 1e-6,
    max_iter = 20
)