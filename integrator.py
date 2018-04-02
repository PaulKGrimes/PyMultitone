# integrator.py
#------------------------
#
# Defines an integrator for calculating KK transforms
# Uses technique from Numerical Recipes
#
# Copyright (c) 2004, Paul Grimes
#
# Derived from code by Phichet Kittara
#
class integrator:
    """Defines a integrator class for use in calculating KK transforms"""
    def __init__(self):
        """Constructor.  Sets numerical constants for integrator"""
        self.__jMax__ = 12
        self.__jMin__ = 5
        self.__eps__ = 1.0e-5
        self.s = 0.0
        
        
    def integrate(self, integrand, low, high):
        """Integration from Numerical Recipes in C 2nd Ed. pp 137 (qtrap)"""
        lastS = -1.0e30
        self.s = 0.0
        
        # This loop is a way of keeping the number of 
        # iterations as low as possible
        for j in range(1, self.__jMax__):
            self.s = self.__trapzd__(integrand, low, high, j)
            if (abs(self.s-lastS) < self.__eps__*abs(lastS)) \
                            and (j>=self.__jMin__):
                return self.s
            lastS = self.s
            
        # if accuracy isn't reached, return __jMax__'th s anyway
        return self.s
        
    def integrate2(self, integrand, low, high):
        """Integration from Numerical Recipes in C 2nd Ed. pp 139 (qsimp)"""
        self.s = 0.0
        st = 0.0
        ost = 0.0
        os = 0.0
        
        # This loop is a way of keeping the number of 
        # iterations as low as possible
        for j in range(1, self.__jMax__):
            st = self.__trapzd__(integrand, low, high, j)
            self.s = (4.0*st-ost)/3.0
            if (abs(self.s-os) < self.__eps__*abs(os)) \
                            and (j>=self.__jMin__):
                return self.s
            os = self.s
            ost = st
            
        # if accuracy isn't reached, return __jMax__'th s anyway
        return self.s
        
            
    def __trapzd__(self, integrand, low, high, n):
         """Integration from Numerical Recipes in C 2nd Ed. pp 137 (trapzd)"""
         
         if (n==1):
             return ((high-low)*integrand(0.5*(low+high)))
         else:
             iteration = 1
             for j in range(1, n-1):
                 iteration *=2
             tnm = iteration    
             delta = (high-low)/tnm
             x = low+0.5*delta
             sum = 0.0
             
             for j in range(iteration):
                 sum += integrand(x)
                 x += delta
                 
             s = 0.5*(self.s+(high-low)*sum/tnm)
             
             return s
                              
                 