# harmonicbalance.py
#------------------------
#
# defines a hermonic balancer class 
#
# Copyright (c) 2004, Paul Grimes
#
# Derived from code by Phichet Kittara
#

class harmonicbalance:
    """Object that balances harmonics"""
    def __init__(self):
        self.StopCheck = FALSE
        
        # Only keep those with Z !=0
        self.HW.clear()
        self.HH.clear()
        
        # Set sizes of data stores and fill with zeros
        self.Old_VN = range
        
        
    def