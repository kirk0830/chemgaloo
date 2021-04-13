"""
Reactor library of Chemgaloo
"""

def bstr(
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
#    if_record = False

# detector_record definition appears here is a little unreasonable. For i do not want to take up many
# memory, this array should be defined elsewhere...
    detector_record = []
#    for detector in detectors:
#        if detector.motion.startswith('record'):
#            if_record = True
#            detector_record = []
#            break

    for ichemi in chemicals:

        chem_c_log.append([ichemi.c])
    for istep in range(nstep):

        # accumulate change of concentration of each chemical by iteratively calling reactions
        for irxn in reactions:

            irxn.go(dt = dt, mode = 'BSTR')
        # save final concentrations to each reaction's self.__c_bak___
        for irxn in reactions:

            irxn.__save_state__()
        # update record of concentrations of reactor
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
    
    return time, chem_c_log, if_expire, detector_record

# old name of bstr used in this file, preserve it for compatibility for old version users
batch_reactor = bstr

import random as rnd
from copy import deepcopy
import math
def __cstr_get_gau1d_val__(x, mu, sigma):

    """
    in-build Gaussian probability distribution function for CSTR
    """
    return 1/(math.sqrt(2*math.pi))/sigma * math.exp(
        -(x-mu)**2/2/sigma**2
    )

def __cstr_c_initializer__(chemicals = [], iter_list = []):

    if len(iter_list) == 0:

        if_iter = False
    else:

        if_iter = True
    
    if if_iter:

        for idx_chemi in range(len(chemicals)):

            chemicals[idx_chemi].ccstr_out = iter_list[idx_chemi]
    else:
        for ichemi in chemicals:

            ichemi.ccstr_out = 1.0
    
def __vector_similar__(v1 = [], v2 = []):

    norm_v1 = 0
    norm_v2 = 0
    corr = 0
    for idx_compo in range(len(v1)):

        norm_v1 += v1[idx_compo]**2
        norm_v2 += v2[idx_compo]**2
        corr += v1[idx_compo] * v2[idx_compo]

    norm_v1 = math.sqrt(norm_v1)
    norm_v2 = math.sqrt(norm_v2)
    return corr/norm_v1/norm_v2

def cstr(
    chemicals = [],
    reactions = [],
    n_thread = 1,
    res_time_dist_f = 'Gaussian',
    rtf_param1 = 0.,
    rtf_param2 = 0.,
    dt = 0.01,
    conv_thr = 1e-3,
    max_iter = 5,
    verbosity = 'high'
):
    """
    WARNING-1: options n_thread > 1 are not implemented! DO NOT USE!\n
    WARNING-2: options res_time_dist_f != 'Gaussian' are not implemented! DO NOT USE!\n
    Continuum-stirring-tank-reactor (CSTR)\n
    Parameters:
    ---
    chemicals: 1-d array, class: chemgaloo.chemical
        chemicals that participate reactions in present reactor
    reactions: 1-d array, class:  chemgaloo.chemical
        reactions that may occur in present reactor
    n_thred: int
        number of perturbation of residual time considered, perturbation obeys residual time function that user defined
    res_time_dist_f: str
        function type of residual time distribution function
    rtf_param(i): float
        parameters of specific residual time distribution function
    dt: float
        time step of discrete integral
    conv_thr: float
        convergence threshold for iteratively solving stable state of CSTR
    max_iter: int
        maximum number of iterations
    verbosity: str
        on which level to print information, availble options: low, medium, high and debug
    """
    if n_thread == 1:

        if res_time_dist_f == 'Gaussian':

            res_time = rtf_param1
            nstep = int(res_time/dt)
        else:
            # other function type is not supported yet
            nstep = 0
        
        chemi_out_prev = []
        for ichemi in chemicals:

            chemi_out_prev.append(ichemi.ccstr_out)
        idx_conv = 0
        iter_err = 1

        if verbosity == 'debug':
            print('Enter CSTR iterative solver...')
        while (iter_err >= conv_thr) and (idx_conv <= max_iter):
            
            for istep in range(nstep):
                    
                # accumulate change of concentration of each chemical by iteratively calling reactions
                for irxn in reactions:

                    irxn.go(dt = dt, mode = 'CSTR')
                # save final concentrations to each reaction's self.__c_bak___
                for irxn in reactions:

                    irxn.__save_state__()
            chemi_out = []
            for ichemi in chemicals:

                chemi_out.append(ichemi.c)
                ichemi.c = ichemi.ccstr_in
            
            __cstr_c_initializer__(chemicals = chemicals, iter_list = chemi_out)

            iter_err = abs(1 - __vector_similar__(
                v1 = chemi_out_prev,
                v2 = chemi_out
                ))
            if verbosity != 'low':
                print('Chemgaloo| (CSTR) iteration {}:\n'.format(idx_conv)
                    +'            conv     = {}\n'.format(iter_err)
                    +'            conv_thr = {}'.format(conv_thr)
                    )
            if (verbosity != 'low') and (verbosity != 'medium'):
                print('Chemgaloo| (CSTR) outflux information:')
                for idx_chemi in range(len(chemi_out)):

                    print('            Chemical-{}, c = {}'.format(idx_chemi+1, round(chemi_out[idx_chemi], 4)))
            if verbosity != 'low':
                print('='*40)
            chemi_out_prev = deepcopy(chemi_out)
            
            idx_conv += 1

    elif n_thread > 1:

        if res_time_dist_f == 'Gaussian':

            pass
        else:

            pass
    else:

        raise ValueError
    pass

def pfr(
    chemicals = [],
    reactions = [],
    c_unit = 'mol/L',
    length = 5,
    dl = 0.05,
    l_unit = 'm',
    velocity = 10,
    v_unit = ''
):
    """
    WARNING: not implemented yet!
    """
    pass