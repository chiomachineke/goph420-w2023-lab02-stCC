import numpy as np
import matplotlib.pyplot as plt
#g= lambda zeta:rho_rat*np.sqrt(zeta_beta_2/(zeta**2) - 1) - np.tan(om*zeta)
#g = lambda zeta: (rho_ratio*np.sqrt(Beta-zeta**2))/(zeta) - np.tan(om*zeta)
def calcuate(rho_ratio,Beta,om,zeta):
    term1=(rho_ratio*np.sqrt(Beta-zeta**2))/zeta
    term2=np.tan(om*f)
    return term1-term2

H=4000
B1=1900
B2=3200
P1=1800
P2=2500
f=[0,0.25,0.5,0.75, 1,1.25,1.5]
rho_ratio= P2/P1
Beta=H**2*(B1**(-2)-B2**(-2))
zeta=np.arange(1,Beta,Beta/1000)
om=2*np.pi*0.5
#y=g(zeta)
#plt.plot(zeta,y)
#plt.ylim(-5,5);
result= calcuate(rho_ratio=10,Beta=20,om=2,zeta=5)