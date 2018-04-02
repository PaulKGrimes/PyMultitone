# multitone.py
#
# Defines class that carries out multitone analysis of SIS mixing
# This class contains everything required to analyse a SIS mixer
#
# This file also contains code to setup a basic simulation and carry
# out a basic analysis
#
# Copyright (C) 2004 Paul Grimes, Astrophysics Group, Cavendish Laboratory
#
# Developed from algortihms and C++ code by Phichet Kittara
#
# Updated 2016, Paul Grimes, Smithsonian Astrophysical Observatory

import scipy.constants as constants
import numpy as np
import numpy.linalg as la
import numpy.convolve as conv

class OneCk:
    """Object holding spectral coefficient and corresponding frequency"""
    def __init__(self):
        self.freq = 1.0
        self.Amp = complex(1.0)


        
class multitone:
    """Object that carries out multitone nonlinear analysis of SIS mixers"""
    def __init__(self):
        """Constructor.  Defines members of multitone object"""
        
        # list of all frequencies
        self.freqs = []
        
        # list of frequencies with Z != 0
        self.freqsToSolve = []
        
        # minimum Ck value to be kept
        self.min_Ck = 1.0e-18
        
        # Ck spectrum (list of OneCk objects)
        self.Ck = []
         
        
    def __setSpectrum__(self):
        """Calculate spectrum by convolving spectra of each frequency"""
        
        self.Ck = []
        TempCk = []
        TempOneCk = OneCk()
        
        # Process 1st frequency
        Tb = self.freqs[1].totalJ
        
        for b in range(-Tb, Tb, 1):
            TempOneCk.freq = b*self.freqs[1].Vph
            TempOneCk.Amp = self.freqs[1].Cjk(b)
            self.Ck.append(TempOneCk)
            
        # Process additional frequencies
        CkSize = len(self.Ck)
        Added = FALSE
        
        for f in range(2, len(self.freqs), 1):
            # Reset temporary variables
            Tb = self.freqs[f].totalJ
            TempCk = []
                
            # Calculate each Ck coefficient
            for b in range(-Tb, Tb, 1):
                for k in range(CkSize):
                    TempOneCk.Amp = Ck[k].Amp * self.freq[f].Cjk(b)
                    
                    # Check to see if Amp is big enough to keep
                    if( abs(TempOneCk.Amp) > self.min_Ck ):
                        Added = FALSE
                        TempOneCk.freq = self.Ck[k].freq + b*self.freqs.Vph
                        
                        # If freq is already in Ck, add new value to old,
                        # if not, add new value and freq to spectrum
                        for c in TempCk:
                            if abs(c.freq-TempOneCk.freq < DOUBLE_PRECISION):
                                c.Amp += TempOneCk.Amp
                                Added = TRUE
                                break
                        
                        if (not Added):
                            TempCk.append(TempOneCk)
                            
            self.Ck = TempCk
            CkSize = len(self.Ck)
            
            
    def __Delta__(self, inX, outY):
        """Calculate the Delta vector"""
        
        TempY = complex(1.0)
        TempI = complex(1.0)
        
        # Set Vn for each harmonic
        for i in range(0, self.HarmonicSize, 2):
            n = i/2 # We want this to truncate
            self.freqs[freqsToSolve[n]].Vn[harmToSolve[n]] = \
                                     complex(inX[i],inX[i+1])
                                     
        self.__setSpectrum__()
        
        # Calculate Delta
        for i in range(0, self.HarmonicSize, 2):
            n = i/2 # Again, we want truncation
            TempI = self.HarmonicCurrent(self.freqsToSolve[n], self.harmsToSolve[n]);
            TempY = self.freqs[self.freqsToSolve[n]].Vs[self.harmsToSolve[n]]
            TempY -=self.freqs[self.freqsToSolve[n]].Z[self.harmsToSolve[n]]*TempI
            TempY -=self.freqs[self.freqsToSolve[n]].Vn[self.harmsToSolve[n]]
            outY[i] = TempY.real
            outY[i+1] = TempY.imaginary
            

    def HarmonicCurrent(self, freq, harm):
        """Return harm'th harmonic of freq'th frequency"""
        f = harm*self.freqs[freq].Wph
        
        return self.Ip(f)
            

    def Ip(self, freq):
        """Return current flowing at frequency freq"""
        
        