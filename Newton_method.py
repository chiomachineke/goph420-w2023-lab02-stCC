import numpy as np
import matplotlib.pyplot as plt
import math
def function(rho_ratio,Beta,phi,zeta):
    term1=(rho_ratio*np.sqrt(Beta-zeta**2))/zeta
    term2=np.tan(phi*zeta)
    return term1-term2


def dfdx(rho_ratio,Beta,phi,zeta):
    term1=-(rho_ratio/(zeta*np.sqrt(Beta-zeta**2)))
    term2=phi*math.cos(phi*zeta) ** -2
    return term1-term2


def newton_raphson_method(zeta0,f,dfdx,eps_s=1*10^-8):
    """finds the root of the equation
    zeta : the initial guess   
    f     : the function we are solving
    dfdx : the second derivative of the function
    """
    eps_a = 2*eps_s
    k=0
    zeta_k=zeta0
    while eps_a > eps_s:
        dx=-f(zeta_k)/dfdx(zeta_k)
        zeta_k  = zeta_k + dx
        k  = k + 1
        eps_a=np.abs(dx/zeta_k)
        
    return zeta_k, k, eps_a

