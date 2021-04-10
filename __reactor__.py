"""
Reactor library of Chemgaloo
"""

def batch_reactor(
    chemicals = [],
    reactions = [],
    c_unit = 'mol/L',
    dt = 0.01,
    nstep = 1000,
    t_unit = 'second',
    detectors = []
):

    """
    Batch reactor: the most common case\n
    Parameters:
    ---
    chemicals: 1-d array, class: chemgaloo.chemical
        all chemicals in reactor
    reactions: 1-d array, class: chemgaloo.reaction
        all possible reactions in reactor
    c_unit: str
        unit of concentration
    dt: float
        minimal time inteval used in discrete integral in numerical simulation
    nstep: int
        total number of steps of present simulation in batch reactor, total time = nstep * dt
    t_unit: str
        unit of time
    detectors: 1-d array, class: chemgaloo.detector
        detectors defined that want to use in present reaction system
    Return:
    ---
    time: 1-d array
        discrete time points in list, useful when want to plot
    chem_c_log: 2-d array
        all chemicals' concentrations recorded at every simulation step: chem_c_log[index_of_chemical][index_of_time]
    if_expire: bool
        state variable, if reactions are quenched, will turn to True
    """

    time = [0.0]
    chem_c_log = []
    if_expire = True
    if_record = False
    for detector in detectors:
        if detector.motion.startswith('record'):
            if_record = True
            detector_record = []
            break

    for ichemi in chemicals:

        chem_c_log.append([ichemi.c])
    for istep in range(nstep):

        for irxn in reactions:

            irxn.go(dt = dt)
        # update concentrations
        for idx_chemi in range(len(chemicals)):

            chem_c_log[idx_chemi].append(chemicals[idx_chemi].c)

        time.append((istep + 1)*dt)

        # Detector motions, may be integrated into another place in future version
        detector_if_print = False
        detector_if_quench = False
        
        for idx_detector in range(len(detectors)):

            if detectors[idx_detector].turn_on():
                print('CHEMGALOO| Detector-{} get expected value at time = {}'.format(idx_detector+1, time[-1]))

                if detectors[idx_detector].motion == 'print':
                    detector_if_print = True
                elif detectors[idx_detector].motion == 'quench':
                    detector_if_quench = True
                elif detectors[idx_detector].motion == 'record':
                    detector_record.append(time[-1])

                detector_if_print = True # print, anyway

                if detectors[idx_detector].motion == 'quench_and_silent':
                    detector_if_quench = True
                    detector_if_print = False
                elif detectors[idx_detector].motion == 'silent':
                    detector_if_print = False
                elif detectors[idx_detector].motion == 'record_and_silent':
                    detector_record.append(time[-1])
                    detector_if_print = False
                elif detectors[idx_detector].motion == 'record_and_quench_and_silent':
                    detector_record.append(time[-1])
                    detector_if_quench = True
                    detector_if_print = False

        if detector_if_print:
            print('CHEMGALOO| ...')
            print('='*40+'\n'+'DETECTOR(S) ACTIVATED, STATE REPORT')
            print('Simulation step = {}, time = {} {}'.format(istep, time[-1], t_unit))
            print('-'*40+'\nConcentration(s):')
            for idx_chemi in range(len(chemicals)):
                print('Chemical-{}, c = {} {}'.format(idx_chemi+1, chem_c_log[idx_chemi][-1], c_unit))
            print('Reaction rate constant(s):')
            for idx_rxn in range(len(reactions)):
                print('Reaction-{}, k = {}'.format(idx_rxn+1, reactions[idx_rxn].k))
            print('='*40)
        if detector_if_quench:
            print('CHEMGALOO| Reactions are quenched according to detector setting...')
            if_expire = False
            break
    
    if if_record:
        return time, chem_c_log, if_expire, detector_record
    else:
        return time, chem_c_log, if_expire

def cstr():
    """
    WARNING: not implemented yet!
    """
    pass

def pfr():
    """
    WARNING: not implemented yet!
    """
    pass