# bessel.py
#------------------------
#
# Defines class for calculating Bessel functions
#
# So far, only simple Bessel values are calculated.  We will implement
# integration and arrays of Jnx later, if required
#
# Copyright (c) 2004, Paul Grimes
#
# Derived from code by Phichet Kittara
#

class bessel:
    def __init__(self):
        """Constructor for the bessel class"""
        
        
    def Jnx(self, n, x):
        """
        Takes integer n and double x.  If non-integer n is passed, it will be
        coerced to integer.  We should add proper exceptions and errors later
        Returns Jn(x) using recurrence relation from Ambramowitz & Stegun
        (1964).  Does this in a python list, which may not be the fastest
        method.
        """
        # Type checking
        n = int(n)
        
        # Beware of these zeros
        if (n==0 and x==0.0):
            return 1.0
        if (n!=0 and x==0.0):
            return 0.0
        
        # Get Bessel order
        abs_n = abs(n)
        
        max_n = 50 + 2*abs_n
       
        # Create array for holding the intermediate recurrence results
        J = range(max_n)
        
        # Start the recurrence relation with small values
        J[max_n-1] = 0.0
        J[max_n-2] = 1.0e-30
        
        # Carry out the recurrence
        
        for jj in range(max_n-2):
            # range runs from 0->max_n-3
            # This is a reverse recurrence from high j to low j
            j = max_n - 3 - jj
            J[j] = 2.*(j+1.)/x*J[j+1] - J[j+2]
            
            # Check for insanity
            if (J[j] > 1.0e10):
                for k in range(j, max_n, 1):
                    J[k] /= 1.0e15
                    
        # Calculate the normalisation factor for the result
        norm_factor = J[0]
        
        for jj in range(2, max_n, 2):
            norm_factor += 2*J[jj] 
        
        if (n>=0):
            return J[n]/norm_factor
        else:
            return pow(-1, (abs_n%2)+2)*J[abs_n]/norm_factor


