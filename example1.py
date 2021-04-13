import chemgaloo as cg
import matplotlib.pyplot as plt

chem1 = cg.chemical(1.)
chem2 = cg.chemical(1.)
chem3 = cg.chemical()
chem4 = cg.chemical()

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

detector1 = cg.detector(
    attr = 'c', 
    mode = 'ratio',
    detect = [chem1, chem2],
    value = [1.0],
    ndigits = 2,
    motion = 'print'
    )

time, c_log, _, _ = cg.reactor.batch_reactor(
    chemicals = [chem1, chem2, chem3, chem4],
    reactions = [rxn1, rxn2],
    dt = 0.01,
    nstep = 1000,
    detectors = [detector1]
    )

plt.plot(time, c_log[0][:])
plt.plot(time, c_log[1][:])
plt.plot(time, c_log[2][:])
plt.plot(time, c_log[3][:])
plt.legend(labels = ['chem1', 'chem2', 'chem3', 'chem4'])
plt.show()