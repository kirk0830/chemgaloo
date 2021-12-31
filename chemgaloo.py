"""
Chemgaloo: a numerical chemical kinetic simulation program\n

"""
h = 6.626E-34
NA = 6.02E23
R = 8.314
kb = R/NA

class chemical:
    """
    chemicals in chemgaloo
    """
    def __init__(self, concentration = 0.0):
        """
        concentration: initial concentration of chemical defined, if not specified explicitly, 0.0 will be given as default
        """
        self.c = concentration
        self.ccstr_in = self.c
        self.ccstr_out = 0.
        # basic stoichemistry
        self.Mr = 0.
        # electronic properties
        self.dipole = 0.
        self.quadrapole = [] # matrix
        self.octapole = [] # matrix
        # structral properties
        self.interia = 0.
        self.point_group = ''
        # thermo-chemistry
        self.U = 0.
        self.H = 0.
        self.S = 0.
        self.F = 0.
        self.G = 0.
        self.E = 0.
        self.vib_freqs = []
        self.ZPE = 0.

    def add(self, amount):

        self.c += amount
    
    def zero_point_ener(self):

        zpe = 0.
        for ivib_freq in self.vib_freqs:

            zpe += 0.5*h*ivib_freq
        self.ZPE = zpe
        return self.ZPE
    
    def internal_ener(self, temperature = 0, if_ideal_gas = True):

        if self.ZPE == 0:

            self.zero_point_ener()
        if if_ideal_gas:

            self.U = self.E + self.ZPE + 1.5*kb*temperature
        
        return self.U
    
    def enthalpy(self, temperature = 0, if_ideal_gas = True):

        if self.U == 0:

            self.internal_ener()
        if if_ideal_gas:

            self.H = self.U + kb*temperature
        return self.H

    def entropy(self, term = 'all', if_ideal_gas = True):

        if term == 't':

            pass
        elif term == 'r':

            pass
        elif term == 'v':

            pass

    def helmholtz_free_ener(self, temperature = 0):

        self.F = self.U - temperature * self.S
        return self.F
    
    def gibbs_free_ener(self, temperature = 0):

        self.G = self.H - temperature * self.S
        return self.G

class reaction:
    """
    reaction definition\n
    Parameters:
    ---
    rxtns: 1-d array, class: chemgaloo.chemical
        reactants that participate in the present reaction\n
    prdts: 1-d array, class: chemgaloo.chemical
        products that produced in the present reaction\n
    stois: 1-d array, int
        stoichemistry, coefficients of chemicals in the present reaction, note that all values
         should be input as POSITIVE values\n
    k: float
        rate constant of present reaction
    OPTIONAL thermo-chemistry parameters:
    ---
    deltaU: float
        internal energy change
    deltaH: float
        enthalpy change
    deltaS: float
        entropy change
    deltaF: float
        Helmholtz free energy change
    deltaG: float
        Gibbs free energy change
    """
    def __init__(self, rxtns = [], prdts = [], stois = [], k = 0.0,
    deltaU = 0., deltaH = 0., deltaS = 0., deltaF = 0., deltaG = 0.):

        self.rxtns = rxtns
        self.prdts = prdts
        self.n_rxtns = len(self.rxtns)
        self.n_prdts = len(self.prdts)
        self.stois = stois
        self.k = k
        self.deltaU = deltaU
        self.deltaH = deltaH
        self.deltaS = deltaS
        self.deltaF = deltaF
        self.deltaG = deltaG

        self.__c_bak__ = []
        self.__save_state__()

    def rate(self, c_from = 'present', mode = 'BSTR'):
        """
        reaction rate calculation\n
        Parameters:
        ---
        c_from: str
            options: 'present' or 'cache'\n
            define where concentrations that used to calculate rate come from, option
            'present' means current concentration will be used, useful when there is only one reaction
            in system, do not use this option when there are many reactions occur. Instead, use option
            'cache', program will calculate rate from concentrations of the last state ("the last state
            " is the same as the state produced by the last iteration for single-reaction cases but
            different for many-reaction cases, where in one iteration, not only one reaction will
            occur but all reactions should occur at the same time, so all reactions should refer
            to the same concentrations set)
        mode: str
            options: 'BSTR' or 'CSTR'\n
            BSTR: Batch-stirring-tank-reactor\n
            CSTR: Continuum-stirring-tank-reactor
        Returns:
        ---
        for case c_from == 'present': a single rate of reaction\n
        for case c_from == 'cache' and mode == 'CSTR': a 1-d array of rate of each chemical
        that participates in present reaction
        """
        if c_from == 'present':
            # ----------------------------------------------------------
            # WARNING: ONLY USE THIS OPTION FOR SINGLE-REACTION SYSTEMS!
            #----------------------------------------------------------
            r = self.k
            for idx_rxtn in range(self.n_rxtns):

                r *= (self.rxtns[idx_rxtn].c)**self.stois[idx_rxtn]
            return r
        elif c_from == 'cache':
            if mode == 'CSTR':

                rs = []
                r0 = self.k
                for idx_rxtn in range(self.n_rxtns):

                    ir = self.k
                    
                    for jdx_rxtn in range(self.n_rxtns):

                        if idx_rxtn != jdx_rxtn:

                            ir *= (self.rxtns[jdx_rxtn].ccstr_out)**self.stois[jdx_rxtn]
                        else:

                            ir *= (self.__c_bak__[jdx_rxtn])**self.stois[jdx_rxtn]
                    rs.append(ir)
                    r0 *= self.__c_bak__[idx_rxtn]**self.stois[idx_rxtn]
                    #r0 *= self.__c_bak__[idx_rxtn]
                rs.append(r0)
                return rs
            elif mode == 'BSTR':

                r = self.k
                for idx_rxtn in range(self.n_rxtns):

                    r *= (self.__c_bak__[idx_rxtn])**self.stois[idx_rxtn]
                return r
        
    def go(self, dt, c_from = 'cache', mode = 'BSTR'):

        """
        evolution for dt of present reaction\n
        Parameter:
        ---
        dt: float
            time inteval, an important parameter used in discreate integral, should not be specified as a too large number!
        """

        if mode == 'BSTR':

            r = self.rate(c_from = c_from, mode = 'BSTR')
            for i in range(self.n_rxtns):

                self.rxtns[i].c -= self.stois[i]*r*dt
            for i in range(self.n_prdts):

                self.prdts[i].c += self.stois[self.n_rxtns + i]*r*dt
        elif mode == 'CSTR':

            rs = self.rate(c_from = c_from, mode = 'CSTR')
            for i in range(self.n_rxtns):

                self.rxtns[i].c -= self.stois[i]*rs[i]*dt
            for i in range(self.n_prdts):

                self.prdts[i].c += self.stois[self.n_rxtns + i]*rs[-1]*dt
    
    def __save_state__(self):

        self.__c_bak__=[] # reset
        for irxtn in self.rxtns:

            self.__c_bak__.append(irxtn.c)
        for iprdt in self.prdts:

            self.__c_bak__.append(iprdt.c)

        return self.__c_bak__
    
    def energy_change(self, psnt = 0., update = False):
        """
        Energy change\n
        Parameters:
        ---
        psnt: float
            extent of chemical reaction
        update: bool
            if update self.deltaU
        Return:
        ---
        deltaU: float
            Energy change at present extent of reaction
        """
        deltaU = 0.
        for idx_rxtn in self.n_rxtns:

            deltaU -= self.stois[idx_rxtn] * self.rxtns[idx_rxtn].U
        for idx_prdt in self.n_prdts:

            deltaU += self.stois[idx_prdt + self.n_rxtns] * self.prdts[idx_prdt].U
        
        if update:

            self.deltaU = deltaU

        return deltaU*psnt

    def enthalpy_change(self, psnt = 0., update = False):
        """
        Enthalpy change\n
        Parameters:
        ---
        psnt: float
            extent of chemical reaction
        update: bool
            if update self.deltaH
        Return:
        ---
        deltaH: float
            Enthalpy change at present extent of reaction
        """
        deltaH = 0.
        for idx_rxtn in self.n_rxtns:

            deltaH -= self.stois[idx_rxtn] * self.rxtns[idx_rxtn].H
        for idx_prdt in self.n_prdts:

            deltaH += self.stois[idx_prdt + self.n_rxtns] * self.prdts[idx_prdt].H

        if update:

            self.deltaH = deltaH
            
        return deltaH*psnt

    def entropy_change(self, psnt = 0., update = False):
        """
        Entropy change\n
        Parameters:
        ---
        psnt: float
            extent of chemical reaction
        update: bool
            if update self.deltaS
        Return:
        ---
        deltaS: float
            Entropy change at present extent of reaction
        """
        deltaS = 0.
        for idx_rxtn in self.n_rxtns:

            deltaS -= self.stois[idx_rxtn] * self.rxtns[idx_rxtn].S
        for idx_prdt in self.n_prdts:

            deltaS += self.stois[idx_prdt + self.n_rxtns] * self.prdts[idx_prdt].S
        
        if update:

            self.deltaS = deltaS
            
        return deltaS*psnt

    def helmholtz_free_ener_change(self, psnt = 0., mode = 'direct', temperature = 0, update = False):
        """
        Helmholtz Free energy change\n
        Parameters:
        ---
        psnt: float
            extent of chemical reaction
        mode: str
            mode of how present property is calculated, available options:
        >direct: directly calculate by F = U - TS\n
        >from_chemicals: calculate from all Fs of chemicals\n
        temperature: float
            only used when `mode` == 'direct', temperature of reactor where reaction occurs
        update: bool
            if update self.deltaF
        Return:
        ---
        deltaF: float
            Helmholtz free energy change at present extent of reaction
        """
        deltaF = 0.
        if mode == 'from_chemicals':

            for idx_rxtn in self.n_rxtns:

                deltaF -= self.stois[idx_rxtn] * self.rxtns[idx_rxtn].F
            for idx_prdt in self.n_prdts:

                deltaF += self.stois[idx_prdt + self.n_rxtns] * self.prdts[idx_prdt].F
        elif mode == 'direct':

            deltaF = self.deltaU - temperature*self.deltaS
        else:

            pass

        if update:

            self.deltaF = deltaF
        
        return deltaF*psnt

    def gibbs_free_ener_change(self, psnt = 0., mode = 'direct', temperature = 0, update = False):
        """
        Gibbs Free energy change\n
        Parameters:
        ---
        psnt: float
            extent of chemical reaction
        mode: str
            mode of how present property is calculated, available options:
        >direct: directly calculate by G = H - TS\n
        >from_chemicals: calculate from all Gs of chemicals\n
        temperature: float
            only used when `mode` == 'direct', temperature of reactor where reaction occurs
        update: bool
            if update self.deltaG
        Return:
        ---
        deltaG: float
            Gibbs free energy change at present extent of reaction
        """
        deltaG = 0.
        if mode == 'from_chemicals':

            for idx_rxtn in self.n_rxtns:

                deltaG -= self.stois[idx_rxtn] * self.rxtns[idx_rxtn].G
            for idx_prdt in self.n_prdts:

                deltaG += self.stois[idx_prdt + self.n_rxtns] * self.prdts[idx_prdt].G
        elif mode == 'direct':

            deltaG = self.deltaH - temperature*self.deltaS
        else:

            pass

        if update:

            self.deltaG = deltaG
            
        return deltaG*psnt


import __reactor__ as reactor

class detector:
    """
    Detector used to monitor reactions in reactor operando\n
    Parameters:
    ---
    attr: str
        short of attribute, define the type of property detected, options:
    >'c': concentration\n
    >'time': reaction time (not implemented yet)\n
    >'temp': temperature of reaction system, need enthalpy data (not implemented yet)\n
    mode: str
        what type of signal, current available options:
    >'isolated': concentration of each chemical defined in keyword detect\n
    >'ratio': ratio of concentrations of chemicals defined in keyword detect\n
    value: 1-d array, float
        expected value(s) that want signal reach
    >`mode` == 'isolated', should give the same number of values as that of elements in keyword detect\n
    >`mode` == 'ratio', only one number is supported at present\n
    ndigits: int
        convergence threshold
    motion: str
        what will detector do if signal is detected the reaches the keyword value, current available options:
    >Batch reactor:\n
    >'print': print information of present state and run continuely\n
    >'quench': print information, and stop\n
    >'quench_and_silent': stop\n
    >'record': save present time point by appending\n
    >'record_and_quench': save present time point by appending and quench reactions\n
    >'record_and_silent': ...\n
    >'record_and_quench_and_silent': ...
    """
    def __init__(
        self, 
        attr,
        mode = '',
        detect = [],
        value = [],
        ndigits = 4,
        motion = 'print'
        ):

        self.attr = attr
        self.mode = mode
        self.detect = detect
        self.value = value
        self.motion = motion
        self.ndigits = ndigits

    def turn_on(self):

        if self.attr == 'c':
            # detect on concentration(s) of chemicals input in array stored in detect
            if self.mode == 'isolated':
                # detect on concentration(s)
                for idx_chemi in range(len(self.detect)):
                    
                    idetect = abs(self.detect[idx_chemi].c - self.value[idx_chemi])
                    if round(idetect, self.ndigits) == 0:
                        return True
                return False
            elif self.mode == 'ratio':
                # detect on ratio of two concentrations
                iratio = self.detect[0].c/self.detect[1].c
                idetect = abs(iratio - self.value[0])
                if round(idetect, self.ndigits) == 0:
                    return True
                else:
                    return False
        elif self.attr == 'time':
            # detect if reach wanted reaction time
            if self.detect[-1] > self.value[0]:
                return True
            else:
                return False
        elif self.attr == 'temp':
            # detect if temperature in reactor reaches the value given
            # enthalpy data is needed...
            pass
        else:
            pass

        return False

