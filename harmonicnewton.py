# harmonicnewton.py
#------------------------
#
# defines a hermonic balancer class 
#
# Copyright (c) 2004, Paul Grimes
#
# Derived from code by Phichet Kittara
#
from numarray import *
import numarray.linear_algebra as linear_algebra

class harmonicnewton:
    """Object that balances harmonics"""
    def __init__(self, HSize):
        """Constructs harmonicnewton object for HSize harmonics"""
        self.StopCheck = FALSE
        
        self.HarmonicSize = HSize
        
        # This needs to start out as the unit matrix
        self.Jacobian = identity(self.HarmonicSize, type=Complex)
        
        self.OldX = zeros(self.HarmonicSize, type=Complex)
        self.NewX = zeros(self.HarmonicSize, type=Complex)
        self.OldY = zeros(self.HarmonicSize, type=Complex)
        self.DeltaX = zeros(self.HarmonicSize, type=Complex)
        
        self.MIN_X = 1.0e-30
        self.MAX_X = 1.0e1
        self.stepFactor = 1.0e-3
        
        
    def clear(self):
        """Clears stored data"""
        self.OldX = zeros(self.HarmonicSize, type=Complex)
        self.NewX = zeros(self.HarmonicSize, type=Complex)
        self.OldY = zeros(self.HarmonicSize, type=Complex)
        self.DeltaX = zeros(self.HarmonicSize, type=Complex)
        
        # This actually needs to be the unit matrix
        self.Jacobian = identity(self.HarmonicSize, type=Complex)
            
        
    def setX(self, initX):
        """Set old and new x arrays to initX"""
        self.OldX = initX
        self.NewX = initX
        
        
    def setJacobian(self):
        """Calculate Jacobian"""
        TempX = zeros(self.HarmonicSize, type=Complex)
        TempY = zeros(self.HarmonicSize, type=Complex)
        
        for h in range(self.HarmonicSize):
            TempX = self.OldX
            
            dX = self.StepFactor*abs(self.OldX[i])
            
            if (dX < self.MIN_X):
                dX = MIN_X
                
            TempX[i] += dX
            FUNCTIONEVALUATION(TempX, TempY)
            
            for j in range(HarmonicSize):
                dY = TempY[j]-OldY[j]
                self.Jacobian[j][i] = dY/dX
                
                
    def setNewX(self):
        """Calculate new X values"""
        
        Jacobian = linear_algebra.inverse(self.Jacobian)
        
        self.DeltaX = matrixmultiply(Jacobian, self.OldY)
        
        self.NewX = self.OldX-self.DeltaX
        
        for h in range(self.HarmonicSize)
            if self.NewX > self.MAX_X:
                NewX[i] = complex(1.0)

                
    def Newton(self, max_it, tol)
        """Carry out the Harmonic Newton minimisation"""
        
        totalCheck = TRUE
        diff = 0.0
        
        count = 0
        
        while( (count < max_it) and totalCheck):
            self.OldX = self.NewX
            #FUNCTIONEVALUATION(self.OldX, self.OldY)
            
            self.setJacobian()
            self.setNewX()
            
            totalCheck = FALSE
            for h in range(self.HarmonicSize):
                diff = abs(self.NewX[i]-self.OldX[i])
                if (diff > epsilon*fabs(OldX[i]):
                    totalCheck = TRUE
                    break
                    
            count++
                     
        return count
            
