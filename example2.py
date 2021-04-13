import chemgaloo as cg
import matplotlib.pyplot as plt

chem1 = cg.chemical()
chem2 = cg.chemical()
chem3 = cg.chemical()
chem4 = cg.chemical()

detector1 = cg.detector(
    attr = 'c', 
    mode = 'ratio',
    detect = [chem3, chem2],
    value = [0.2],
    ndigits = 4,
    motion = 'quench_and_silent'
    )

detector2 = cg.detector(
    attr = 'c',
    mode = 'isolated',
    detect = [chem1],
    value = [0.01],
    ndigits = 4,
    motion = 'quench'
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
k2s = linspace(start = 1.0, stop = 5.0, num = 10)

time_reach_record = []
k_reach_record = []
show_each_k2 = False
for k2 in k2s:

    rxn2 = cg.reaction(
        rxtns = [chem1, chem3],
        prdts = [chem4],
        stois = stois2,
        k = k2
    )

    chem1.c = 1.0
    chem2.c = 1.0
    chem3.c = 0.0
    chem4.c = 0.0

    time, c_log, if_expire, _ = cg.reactor.batch_reactor(
        chemicals = [chem1, chem2, chem3, chem4],
        reactions = [rxn1, rxn2],
        dt = 0.001,
        nstep = 5000,
        detectors = [detector2]
        )

    if not if_expire:
        time_reach_record.append(time[-1])
        k_reach_record.append(k2)
    
    if show_each_k2:
        plt.plot(time, c_log[0][:])
        plt.plot(time, c_log[1][:])
        plt.plot(time, c_log[2][:])
        plt.plot(time, c_log[3][:])
        plt.legend(labels = ['chem1', 'chem2', 'chem3', 'chem4'])
        plt.show()

plt.plot(time_reach_record, k_reach_record, 'r-')
plt.plot(time_reach_record, k_reach_record, 'ro')
plt.xlabel('Time to reach ratio {}'.format(detector1.value[0]))
plt.ylabel('k2 rate constant')
plt.show()