# reponseFn.py
#------------------------
#
# Defines a response function class that holds an IV curve and its KK transform
#
# Copyright (c) 2004, Paul Grimes
#
# Derived from code by Phichet Kittara
#
import integrator, string
from constants import pi

class responseFn:
    """Class containing an SIS response function"""
    def __init__(self):
        """
        Constructs the SIS response function class
        """
        self.Igap = 1.0
        self.Vgap = 1.0
        self.Rn = 1.0
        self.yIntercept = 0.0
        self.noPoints = 201
        self.__bias__ = 0.0 # used in calculating the KK transform
        self.__ikk__ = range(self.noPoints)
        self.__idc__ = range(self.noPoints)
        self.__vdc__ = range(self.noPoints)
        self.__Int__ = integrator.integrator()
        self.__Int__.__jMax__ = 12
        self.__KK_vMax__ = 10.0
        self.separator = "\t" # data value separator in idc and ikk files
        
        
    def __bubbleFind__(self, testX, yList, xList, a, b):
        """
        Returns interpolated value of yList at x position given by testX.
        a and b are integers giving the extremes of the index range to be searched.
        For normal use, these should be the extents of the data lists
        """
        if ( abs(b-a)<=1 ): # We are within two points of the data we want.
            if( a==b ):     # Dead on x value
                return yList[a]
            else:           # Interpolate
                return yList[a] + (yList[b]-yList[a])*(testX-xList[a]) \
                                        /(xList[b]-xList[a])
                                        
        else: # Binary search until we just enclose the X value we want
            c = a+(b-a)/2
            if( xList[c] < testX ): # New index is below required index
                return self.__bubbleFind__(testX, yList, xList, c, b)
            else:
                return self.__bubbleFind__(testX, yList, xList, a, c)
        
        
    def Ikk(self, bias):
        """
        Returns the value of the KK transform at bias.
        If bias is outside of data range, returns final value of __ikk__
        """
        if (abs(bias) >= self.__vdc__[-1]): # Outside of bias data range
            return self.__ikk__[-1]
        else:
            return self.__bubbleFind__(bias, self.__ikk__, self.__vdc__, \
                                                0, len(self.__ikk__)-1)
        
        
    def Idc(self, bias):
        """Returns the DC current at bias"""
        if (abs(bias) > self.__vdc__[-1]): # Outside of bias data range
            result = self.yIntercept + self.Rn*abs(bias) # Linear extrapolation
        else:
            result = self.__bubbleFind__(bias, self.__idc__, self.__vdc__, \
                                                0, len(self.__idc__)-1)
        if (bias < 0):
            result = -result
            
        return result
        
             
    def Vdrive(self, n):
        """
        Returns the Bias Voltage with index n in the __vdc__ table.
        No checks are performed on range on n"""
        return self.__vdc__[n]
    
        
    def ReadData(self, idcFileName, ikkFileName):
        """
        Reads IV and KK data from the given filenames.
        Non numerical data is ignored and values are assumed to be tab
        separated.
        """
        # Get the data
        idcFile = open(idcFileName)
        idcLines = idcFile.readlines()
        idcFile.close()
        
        ikkFile = open(ikkFileName)
        ikkLines = ikkFile.readlines()
        ikkFile.close()
        
        # Process the IV data
        idcV = []
        idcI = []
        for line in idcLines:
            l = string.split(line, self.separator)
            try:
                idcV.append(float(l[0]))
                idcI.append(float(l[1]))
            except: # Something bad happened, and we don't care what
                continue
            
        self.__vdc__ = idcV
        self.__idc__ = idcI
        
        # Process KK data, which is probably shorter than IV data
        ikkV = []
        ikkI = []
        for line in ikkLines:
            l = string.split(line, self.separator)
            try:
                ikkV.append(float(l[0]))
                ikkI.append(float(l[1]))
            except: # Something bad happened, and we don't care what
                continue
            
        # Match KK data to IV data
        self.__ikk__ = range(len(self.__idc__))
        for n in range(len(self.__vdc__)):
            v = self.__vdc__[n]
            self.__ikk__[n] = self.__bubbleFind__(v, ikkI, ikkV, 0, len(ikkV)-1)
            
        self.noPoints = len(self.__vdc__)
        
        # Calculate Rn and yIntercept so we can extrapolate
        Rn = 0.0
        Sn = 0.0
        Vmid = 0.0
        Imid = 0.0
        
        # Find top vGap/5.0 of bias range
        start = len(idcV) - 1
        
        for i in range(len(idcV)):
            if ((idcV[-1] - idcV[start]) < self.Vgap/5.0):
                start -= 1
            else:
                break
        
        for i in range(start, len(idcV)):
            Rn += (idcV[i]-idcV[i-1])/(idcI[i]-idcI[i-1])
            Sn += (idcI[i]-idcI[i-1])/(idcV[i]-idcV[i-1])
            Vmid += idcV[i]
            Imid += idcI[i]
            
        Rn /= float(len(idcV)-start)
        Vmid /= float(len(idcV)-start)
        Imid /= float(len(idcV)-start)
        
        self.Rn = Rn
        self.yIntercept = Imid - Vmid/self.Rn
        
        
    def Kennedy(self, n, maxBias=2.0, points=201):
        """
        Generates a Kennedy fit approximation to an SIS IV curve.
        n is the order of the polynomial and points sets the resolution of the 
        data.
        """
        # Clear and set length of data lists
        self.__vdc__ = range(points)
        self.__idc__ = range(points)
        self.__ikk__ = range(points)
        
        for k in range(points):
            v = maxBias * k / (points-1.0)
            self.__vdc__[k] = v
            self.__idc__[k] = pow(v, 2*n+1)/(1.0+pow(v,2*n))
        
        self.noPoints = points
        self.Rn = 1.0
        self.yIntercept = 0.0
        self.__calc_Ikk__()
        
        
    def __calc_Ikk__(self):
        """Calculates the KK transform of the current IV data"""
        for i in range(self.noPoints):
            self.__bias__ = self.__vdc__[i]
            self.__ikk__[i] = self.__Int__.integrate(self.__KK_integrand__, \
                                    0.0, self.__KK_vMax__)/pi
            print self.__bias__, ":", self.__ikk__[i]

                        
    def __KK_integrand__(self, v):
        """Calculates the integrand used in the KK transform"""
        G1 = (self.Idc(self.__bias__+v) - (self.__bias__+v)) / v
        G2 = (self.Idc(self.__bias__-v) - (self.__bias__-v)) / (-v)
        return G1+G2