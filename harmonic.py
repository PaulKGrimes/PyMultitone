# harmonic.py
#------------------------
#
# Defines a class holding the properties of a single harmonic
#
# Class contains:
#    bessel b     : object for calculating bessel functions.  This is a very
#                    inefficient idea, as each harmonic will have it's own
#                    bessel engine.
#    list Anp     : list of Anp coefficients
#    int p        : order of harmonic
#    int totalJ   : number of Anp coefficients to calculate
#    complex Cvp  : Cvp coefficient of this class
#    double w_0  : angular frequency of harmonics fundamental
#    double w_Gap : angular frequency of junction gap (normalising frequency)
#    calc_Anp()   : method to cause recalculation of Anp coefficients
#    get_Anp()    : method to retrieve value of Anp coefficient
#
# Copyright (c) 2004, Paul Grimes
#
# Derived from code by Phichet Kittara
#
import bessel

class harmonic:
    """Class keeping all information about a single harmonic"""
    def __init__(self, totalJ, p, Vp, w0, w_Gap):
        """
        Constructor takes five values:
            totalJ : int number of Anp coefficients
            p : int order of harmonic
            Vp : complex coefficient Vp
            w_0 : double angular frequency of fundamental (harmonic where p=0)
            w_Gap : normalising angular frequency
        """
        self.totalJ = totalJ
        self.p = p
        self.Vp = Vp
        self.w0 = w0
        self.w_Gap = w_Gap
        
        self.Anp = range(totalJ*2 + 1)
        self.b = bessel.bessel()
        
        # Must do this last        
        self.calc_Anp()

        
    def calc_Anp(self):
        """
        Calculate Anp values.  These persist until this method is called again
        """
        # Set __ protected values to store state at time of calculation
        self.__last_totalJ__ = self.totalJ
        self.__last_p__ = self.p
        self.__last_Vp__ = self.Vp
        self.__last_w0__ = self.w0
        self.__last_w_Gap__ = self.w_Gap
        
        # Set length of Anp
        self.Anp = range(self.totalJ*2 + 1)
        
        # Get magnitude of voltage
        V = abs(self.Vp)
        # calculate junction drive level
        # For debug porpoises
        alpha = V*self.w_Gap/(self.p*self.w0)
        
        # Set coefficient on Bessel value in Anp
        if (V > 0.0): 
            unitV = self.Vp/V # i.e. unit vector in direction of Vp
        else:
            unitV = complex(1.0,0.0)
        
        real_unitV = pow(unitV.conjugate(), -self.totalJ)
        
        for i in range(-self.totalJ, self.totalJ+1):
            self.Anp[i+self.totalJ] = self.b.Jnx(i, alpha)
            self.Anp[i+self.totalJ] *= real_unitV
            # real_unitV *= unitV.conjugate()
            # self.Anp[i+totalJ] = pow(unitV.conjugate(), i)
            
            
    def get_Anp(self, n):
        """Returns value of Anp, given integer n in range
        -totalJ -> totalJ.  Only recalculates Anp if state has changed"""
        if (n < -self.totalJ) or (n > self.totalJ):
            return 0.0
        else:
            if (not self.valid()):
                self.calc_Anp()
        
            return self.Anp[n+self.totalJ]
        
        
    def valid(self):
        """Checks if calculated Anp values are valid"""
        if (self.__last_totalJ__ != self.totalJ or
            self.__last_p__ != self.p or
            self.__last_Vp__ != self.Vp or
            self.__last_w0__ != self.w0 or
            self.__last_w_Gap__ != self.w_Gap):
            return False
        else:
            return True
  