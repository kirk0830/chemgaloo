"""
Chemgaloo: a numerical chemical kinetic simulation program\n

"""

class chemical:
    """
    chemicals in chemgaloo
    """
    def __init__(self, concentration = 0.0):
        """
        concentration: initial concentration of chemical defined, if not specified explicitly, 0.0 will be given as default
        """
        self.c = concentration

    def add(self, amount):

        self.c += amount
    
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
    """
    def __init__(self, rxtns = [], prdts = [], stois = [], k = 0.0):

        self.rxtns = rxtns
        self.prdts = prdts
        self.n_rxtns = len(self.rxtns)
        self.n_prdts = len(self.prdts)
        self.stois = stois
        self.k = k

    def rate(self):

        r = self.k
        for irxtn in range(self.n_rxtns):

            r *= (self.rxtns[irxtn].c)**self.stois[irxtn]
        return r

    def go(self, dt):

        """
        evolution for dt of present reaction\n
        dt: float
            time inteval, an important parameter used in discreate integral, should not be specified as a too large number!
        """
        r = self.rate()
        for i in range(self.n_rxtns):

            self.rxtns[i].c -= self.stois[i]*r*dt
        for i in range(self.n_prdts):

            self.prdts[i].c += self.stois[self.n_rxtns + i]*r*dt

import __reactor__ as reactor

class detector:
    """
    Detector used to monitor reactions in reactor operando\n
    Parameters:
    ---
    attr: str
        short of attribute, define the type of property detected, current available options:
    >'c': concentration\n
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

