import chemgaloo as cg
import matplotlib.pyplot as plt

chem1 = cg.chemical()
chem2 = cg.chemical()
chem3 = cg.chemical()
chem4 = cg.chemical()

detector1 = cg.detector(
    attr = 'c', 
    mode = 'ratio',
    detect = [chem4, chem3],
    value = [0.1],
    ndigits = 1,
    motion = 'record_and_silent'
    )

detector2 = cg.detector(
    attr = 'c', 
    mode = 'ratio',
    detect = [chem4, chem3],
    value = [0.7],
    ndigits = 1,
    motion = 'record_and_quench_and_silent'
    )

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
from numpy import linspace

stois2 = [1, 1, 1]
k2s = linspace(start = 1.0, stop = 5.0, num = 20)

k_reach_record = []
delta_ts = []

for k2 in k2s:

    rxn2 = cg.reaction(
        rxtns = [chem1, chem3],
        prdts = [chem4],
        stois = stois2,
        k = k2
    )

    chem1.c = 5.0
    chem2.c = 1.0
    chem3.c = 0.0
    chem4.c = 0.0

    time, c_log, if_expire, record = cg.reactor.batch_reactor(
        chemicals = [chem1, chem2, chem3, chem4],
        reactions = [rxn1, rxn2],
        dt = 0.001,
        nstep = 5000,
        detectors = [detector1, detector2]
        )

    print('k2 = {}'.format(k2))
    print(record)
    if len(record) >= 2:
        delta_ts.append(record[-1] - record[-2])
        k_reach_record.append(k2)

print(delta_ts)
plt.plot(k_reach_record, delta_ts, 'r-')
plt.plot(k_reach_record, delta_ts, 'ro')
plt.xlabel('k2 rate constant')
plt.ylabel('time inteval from raio {} to ratio {}'.format(detector1.value[0], detector2.value[0]))
plt.show()
