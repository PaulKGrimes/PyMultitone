# frequency.py
#------------------------
#
# Defines a frequency class that holds a collection of harmonics.
# We need things from TPropertyForm as well as One_Frequency
#
# Copyright (c) 2004, Paul Grimes
#
# Derived from code by Phichet Kittara
#
import harmonic
from constants import *

class frequency:
    """Defines a frequency object that combines a number of harmonics"""
    def __init__(self):
        """frequency objects take no parameters"""
        
        # Default set of required properties
        self.totalJ = 10
        self.Vs = [complex(0.0), complex(1.0)]
        self.Z = [complex(0.0), complex(1.0)]
        self.Vn = [complex(0.0), complex(1.0)]
        self.w0 = 0.5     # Fundamental frequency
        self.w_Gap = 1.0  # Junction gap (normalising) frequency
        self.Vph = self.w0/self.w_Gap
        
        # Some constants
        self.Delta_Ck = 1.0e-13
        self.Delta_Ck_0 = 1.0e-13
        
        # create list of harmonics with one harmonic in it
        h1 = harmonic.harmonic(self.totalJ, 1, self.Vs[1], self.w0,\
                                        self.w_Gap)
        
        self.harmonics = [0, h1]
        
        # Vectors used for storage
        self.DBCjk = range(self.totalJ*2+1)
        self.DBC2k = range(self.totalJ*2+1)
        for Cjk in self.DBCjk:
            Cjk = complex(0.0)
        for C2k in self.DBC2k:
            C2k = complex(0.0)
            
        # Private variables used for consistency checking
        self.__last_numHarmonics__ = 0
        self.__last_w0__ = 0.5
        self.__last_w_Gap__ = 1.0
        self.__last_Vs__ = self.Vs
        self.__last_Z__ = self.Z
        self.__last_Vn__ = self.Vn
        self.__last_Vph__ = self.Vph
        
        # For debug porpoises
        print "Calling f.__set_Ap__()..."
        self.__set_Ap__()
        
        print "Calling f.__set_Cjk__()..."
        self.__set_Cjk__()
        
        
    def getNumHarmonics(self):
        """Returns total number of harmonics in use"""
        return len(self.harmonics) - 1
    
    
    def addHarmonic(self, Vs=complex(0.0), Z=complex(0.0)):
        """
        Adds a harmonic to the frequency.  Returns total number of harmonics
        If Vs and Z are not set the new harmonic will be added with
        zero source voltage and embedding impedance.  These parameters
        must be set subsequently for the harmonic to have an effect
        """
        self.harmonics.append(harmonic.harmonic(self.totalJ, \
                                self.getNumHarmonics()+1, Vs, self.w0, \
                                self.w_Gap))
        self.Vs.append(Vs)
        self.Z.append(Z)
        self.Vn.append(complex(1.0))
        
        self.__set_Ap__()
        self.__set_C2k__()
        self.__set_Cjk__()
        
        return self.getNumHarmonics()
       
    
    def delHarmonic(self):
        """
        Removes the last harmonic from the frequency.  Returns number of
        harmonics.
        The last harmonic must be removed as there must be a continuous set of
        harmonics
        """
        self.harmonics.pop()
        self.Vs.pop()
        self.Z.pop()
        self.Vn.pop()
        
        return self.getNumHarmonics()
    
        
    def getNumBessel(self):
        """Returns number of Bessel coefficients per harmonic"""
        return self.totalJ
    
    def setNumBessel(self, totalJ):
        """
        Sets the number of Bessel coefficients per harmonic.  Does not
        cause any recalculation within the harmonics, as this will be done
        automatically when harmonic.get_Anp(n) is called
        """
        for h in harmonics:
            h.totalJ = totalJ
        
        self.totalJ = totalJ    
    
        
    def getFrequency(self):
        """Returns frequency of fundamental"""
        return self.w0
    
    def setFrequency(self, w0):
        """Set fundamental frequency"""
        self.w0 = w0
        self.Vph = w0/self.w_Gap
            
        
    def getVs(self, n):
        """Returns source voltage at harmonic n"""
        return self.Vs[n]
    
    def setVs(self, n, v):
        """Set Vs of harmonic n to v"""
        self.Vs[n] = complex(v)
        return 0
        
    
    def getZ(self, n):
        """Returns embedding impedance of harmonic n"""
        return self.Z[n]
   
    def setZ(self, n, z):
        """Set Z of harmonic n to z"""
        self.Z[n] = complex(z)
        return 0     
        
        
    def getVn(self, n):
        """Returns nonlinear voltage at harmonic n"""
        return self.Vn[n]
   
    def setVn(self, n, v):
        self.Vn[n] = complex(v)
        return 0
            
    
    def checkValid(self):
        """
        See if anything has changed since last calculation.
        Used in all private calculations that don't return anything
        """
        if (self.__last_numHarmonics__ == self.getNumHarmonics() ) and \
           (self.__last_w0__ == self.w0) and \
           (self.__last_w_Gap__ == self.w_Gap) and \
           (self.__last_Vs__ == self.Vs) and \
           (self.__last_Z__ == self.Z) and \
           (self.__last_Vn__ == self.Vn):
            return True
        else:
            return False
    
    
    def Cjk(self, j, k):
        """Returns value of Cjk"""
        return self.__get_each_Cjk__(j,k)
        
        
    def Ck(self, k):
        """Returns value of Cjk"""
        result = complex(0.0)
        if (abs(k) <= self.totalJ):
            result = self.DBCjk[k+self.totalJ]
        
        return result
        
        
    def Ip(self, responseFn, x, n):
        """
        Calculates and returns current flowing at harmonic n when biased
        at x.  1st argument is a response function object
        """
        rs_plus = complex(0.0)
        rs_minus = complex(0.0)
        
        for k in range(-self.totalJ, self.totalJ+1):
            Idc = responseFn.Idc(x + k*self.Vph)
            Ikk  = responseFn.Ikk(x + k*self.Vph)
            
            Ires = complex(Ikk, Idc)
            C0 = self.Ck(k)
            C_plus = self.Ck(k+n).conjugate()
            C_minus = self.Ck(k-n).conjugate()
            
            rs_plus += C0*C_plus*Ires
            rs_minus += C0*C_minus*Ires
            
        result = rs_minus - rs_plus.conjugate()
        
        if (n==0):
            result /= 2.0
            
        return complex(result.imag, result.real)
        
        
        
    def save_data(self, f):
        """Writes frequency data out to file object f"""
        pass
        
        
    def __set_Cjk__(self):
        """Private method to set Cjk values"""
        if self.checkValid():
            return 0
        
        # Get Cjk
        for k in range(-self.totalJ, self.totalJ+1):
            self.DBCjk[k+self.totalJ] = \
                self.__get_each_Cjk__(self.getNumHarmonics(), k)
                
        # check accuracy
        sum = complex(0.0)
        for k in range(-self.totalJ, self.totalJ+1):
            sum += self.DBCjk[k+self.totalJ] \
                    * self.DBCjk[k+self.totalJ].conjugate()
            sum += self.DBCjk[k+self.totalJ].conjugate() \
                    * self.DBCjk[k+self.totalJ]
                   
        sum *= 0.5
        
        self.Delta_Ck = abs(sum-1.0)
        
        # Accurate enough
        if (self.Delta_Ck > self.Delta_Ck_0) and (self.totalJ < 50):
            self.totalJ += 5
            self.DBCjk = range(self.totalJ*2+1) # magic to get list of right
            self.DBC2k = range(self.totalJ*2+1) # length and type
            for Cjk in self.DBCjk:
                Cjk = complex(0.0)
            for C2k in self.DBC2k:
                C2k = complex(0.0) # end of magic
            self.__set_Ap__()
            self.__set_Cjk__()
            
        # Too accurate
        if (self.Delta_Ck < 1.0e-5*self.Delta_Ck_0) and (self.totalJ > 15):
            self.totalJ -= 2
            self.DBCjk = range(self.totalJ*2+1) # magic to get list of right
            self.DBC2k = range(self.totalJ*2+1) # length and type
            for Cjk in self.DBCjk:
                Cjk = complex(0.0)
            for C2k in self.DBC2k:
                C2k = complex(0.0) # end of magic
            self.__set_Ap__()
            self.__set_Cjk__()
            
        
        
    def __set_C2k__(self):
        """Private method to set C2k values"""
        #if self.checkValid():
        #    return 0
        
        for k in range(-self.totalJ, self.totalJ+1):
            temp = complex(0.0)
            for m in range(-self.totalJ, self.totalJ+1):
                # This bit goes very wrong as k-2m runs from 10->-30
                # then 11->-29, which overflows h.Anp[k-2m]
                # Missing check on value of int passed to h.get_Anp
                temp +=self.harmonics[1].get_Anp(k-2*m) \
                     * self.harmonics[2].get_Anp(m)
                        
            self.DBC2k[k+self.totalJ] = temp
            
        
        
    def __set_Ap__(self):
        """Private method to set Ap values"""
        if self.checkValid():
            return 0
        
        for p in range(1, self.getNumHarmonics()):
            self.harmonics[p].totalJ = self.totalJ
            self.harmonics[p].p = p
            self.harmonics[p].Vp = self.Vn[p]
            self.harmonics[p].w0 = self.w0
            self.harmonics[p].w_Gap = self.w_Gap
            self.harmonics[p].Anp = range(self.totalJ*2 + 1)
            self.harmonics[p].calc_Anp()
            
        if (self.getNumHarmonics() > 1 ):
            self.__set_C2k__()
            
        
    def __get_each_Cjk__(self, j, k):
        """Private method to return individual Cjk values"""
        result = complex(0.0)
        if ( j<=0 ) or ( j > self.getNumHarmonics()):
            raise IndexError, "Requested non-existant harmonic"
        elif ( j==1 ):
            result = self.harmonics[j].get_Anp(k)
        elif ( j==2 ):
            if (abs(k) <= self.totalJ):
                result = self.DBC2k[k + self.totalJ]
        
        else: # j >=3, Only want to run this if harmonics exist
            for m in range(-self.totalJ, self.totalJ+1):
                result += self.__get_each_Cjk__(j-1, k-j*m) \
                            * self.harmonics[j].get_Anp(m)

        return result 
        
