import numpy as np
import matplotlib.pyplot as plt
def function(rho_ratio,Beta,om,zeta):
    term1=(rho_ratio*np.sqrt(Beta-zeta**2))/zeta
    term2=np.tan(om*zeta)
    return term1-term2




def secant_method(x1,x2,f,eps_s=1*10^-8):
    """finds the root of the equation
    x1,x2 : the interval used while searching for the root;
            f(x1)*f(x2)<0   
    f     : the function we are solving
    """
    eps_a = 2*eps_s
    k=0
    while eps_a > eps_s:
        fx1=f(x1)
        fx2=f(x2)
        xk=x2-((x1-x2)*fx1)/(fx1-fx2)
        dx=xk-x2
        #k + = 1     