import numpy as np

"""
Basic File to Keep Track of Financial Vehicles available

All Durations are modified durations, not Macauley Durations.

#TODO:Check all convexity formulas. 

"""

# Present Values 
def present_value_level_annuity_immediate(n:int,
                                          i:float):
    """
    Present Value of a level annuity immediate

    Parameters:
    -----------
    n: int
        Number of payments in the annuity.

    i: float
        Interest rate

    Returns:
    -------
    PV : float
        Present value of the annuity

    """
    nu = 1/(1+i)
    PV = (1-np.power(nu,n))/i
    return PV

def present_value_arithmetic_annuity_immediate(n:int,
                                               i:float):
    """
    Present Value of a arithmetically increasing
    annuity immediate with the first value of 1 
    and increasing by 1 every step. 

    Parameters:
    -----------
    n: int
        Number of payments in the annuity.

    i: float
        Interest rate

    Returns:
    -------
    PV : float
        Present value of the annuity

    """
    an = present_value_level_annuity_immediate(n=n,i=i)
    nu = 1/(1+i)
    PV = ((1+i)*an - n*np.power(nu,n))/i
    return PV

def present_value_quadratic_annuity_immediate(n:int,
                                              i:float):
    """
    Computes the present value of an asset that grows 
    quadratically starting at 1 i.e. 
        \sum_{m=1}^n m(m+1) nu^m

    Parameters:
    -----------
    n: int
        Number of payments in the annuity.

    i: float
        Interest rate

    Returns:
    -------
    PV : float
        Present value of the annuity
    """
    #TODO: Test This is correct
    first = (2*np.power(i+1,2))/np.power(i,3)
    second = (np.power(i,2)*(n+1)*(n+2) + i*2(n+2) + 2)/(np.power(1+i,n)*np.power(i,3))
    PV = first+second
    return PV

def present_value_perpetuity_immediate(i:float):
    """
    Present value of a perpetuity immediate

    Parameters:
    -----------
    n: int
        Number of payments in the annuity.

    i: float
        Interest rate

    Returns:
    -------
    PV : float
        Present value of the annuity
    """
    PV = 1/i
    return PV

def present_value_arithmetic_perpetuity_immediate(i:float):
    """
    Present value of an arithmetically increasing
    perpetuity

    Parameters:
    -----------
    i: float
        Interest rate

    Returns:
    -------
    PV : float
        Present value of the annuity
    
    """
    PV = (1+i)/np.power(i,2)
    return PV

def present_value_quadratic_perpetuity_immediate(i:float):
    """
    Present value of a quadratically increasing perpetuity 
    immediate i.e.
        \sum_{m=1}^\infty m(m+1)nu^m
    where nu is the discount factor

    Parameters:
    -----------
    i: float
        Interest rate

    Returns:
    -------
    PV : float
        Present value of the annuity
    """
    PV = 2*np.power(1+i,2)/np.power(i,3)
    return PV


def present_value_lump(n:int,
                       i:float):
    """
    Present value of a lump sum

    Parameters:
    -----------
    n: int
        Number of payments in the annuity.

    i: float
        Interest rate

    Returns:
    -------
    PV : float
        Present value of the annuity
    """
    PV  =  np.power(1+i,-n)
    return PV

class Abstract_Financial_Instrument:
    def __init__(self) -> None:
        pass
    def compute_PV(self,
                   i:float):
        pass
    def compute_duration(self,
                         i:float):
        pass
    def compute_convexity(self,
                          i:float):
        pass


class ZeroCouponBond(Abstract_Financial_Instrument):
    def __init__(self,
                 n:int,
                 C:float):
        self.C = C
        self.n = n
    def compute_PV(self,
                   i:float):
        return self.C*present_value_lump(n=self.n,
                                         i=i)
    def compute_duration(self,
                         i:float):
        return self.n/(1+i)
    def compute_convexity(self,
                          i:float):
        return  self.n*(self.n+1)/np.power(1+i,2)
    
class CouponBond(Abstract_Financial_Instrument):
    def __init__(self,
                 n:int,
                 r:int,
                 F:float,
                 C:float):
        self.n = n
        self.r = r
        self.F = F
        self.C = C
    def compute_PV(self,
                   i:float):
        an = present_value_level_annuity_immediate(n=self.n,
                                                   i=i)
        p = present_value_lump(n=self.n,
                               i=i)
        PV = self.F*self.r*an + self.C*p
        return PV
    def compute_duration(self,
                         i:float):
        num = self.F*self.r*present_value_arithmetic_annuity_immediate(n=self.n,
                                                                       i=i)
        num += self.C*self.n*present_value_lump(n=self.n,i=i)
        denom =self.compute_PV(i=i)
        D = num/((1+i)*denom)
        return D
    def compute_convexity(self,
                          i: float):
        num = self.F*self.r*present_value_quadratic_annuity_immediate(n=self.n,
                                                                      i=i)
        num+= self.C*(self.n)*(self.n+1)*present_value_lump(n=self.n,
                                                            i=i)
        denom = self.compute_PV(i=i)
        return num/(denom*np.power(1+i,2))
    
class Pepetuity(Abstract_Financial_Instrument):
    def __init__(self,
                 P:float):
        self.P = P
    
    def compute_PV(self, i: float):
        return super().compute_PV(i)
    
    def compute_duration(self, i: float):
        return super().compute_duration(i)
    
    def compute_convexity(self, i: float):
        return super().compute_convexity(i)
        

        
        

    

