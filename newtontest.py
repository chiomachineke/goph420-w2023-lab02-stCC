import numpy as np              #importing the numpy and math libraries
import matplotlib.pyplot as plt
import math
def newton_raphson_method(zeta0,f,dfdx):                    #defining the root finding method
    """finds the root of the equation
    zeta0 : the initial guess
            float   
    f    : the function we are solving
    dfdx : the first derivative of the function

    pseudo code
    1: evaluate f and dfdx at zeta0
    2: evaluate zeta1 from f and dfdx
    3: the iteration stops if eps_a < eps_s, if not, the iterations continues until eps_a< eps_s
        """
    
    zeta_k = zeta0          #setting the initial true root equal to the initial guess
    eps_s=1e-8              #stopping criterium
    eps_a = 2*eps_s
    k=0                     #iteration counter set initially at zero
    eps_k= []

    while eps_a > eps_s:        #the algorithm stops when relative error is less than the stopping criterium
        k  = k + 1              # iteration counter increases by 1 after each iteration
        dx=-f(zeta_k)/dfdx(zeta_k)      
        zeta_k  = zeta_k + dx        
        eps_a=np.abs(dx/zeta_k)
        #eps_a=np.abs(dx)
  
        eps_k.append(eps_a)
    return zeta_k, k, eps_a

#%% Question 2

# function parameters
P1=1800
P2=2500
H=4000
B1=1900
B2=3200

zeta_max = H * np.sqrt(((B1**-2)-(B2**-2)))

#print(zeta_max)
f_list=[]           #an empty list for frequency values
zeta_r_lists=[]     #an empty list for zeta values

#estimating the frequency values that would give zeta values within the range of zeta_max
for n in range(10):
    f= 0.25*(2*n+1)/zeta_max
    f_list.append(f)
    zeta_r_list=[]

    #defining the given function
    def function(zeta0):
        Beta=(H**2*((B1**-2)-(B2**-2)))

        rho_ratio= P2/P1
        Beta=(H**2)*((B1**(-2))-(B2**(-2)))

        phi=2*np.pi*f

        term1 = (rho_ratio * np.sqrt(Beta - (zeta0 ** 2))) / zeta0
        term2 = np.tan(phi * zeta0)
        result = term1 - term2
        
        return result

    #first derivative of the given function
    def dfdx(zeta0):
    
        dz = 1e-10
        result2 = (function(zeta0+dz) - function(zeta0))/(dz)
        
        return result2

    for k in range(n+1):        # estimating the zeta values 
        zeta_k=0.25*(2*k+1)/f
        zeta0=zeta_k-1e-3
        zeta_r=newton_raphson_method(zeta0,function,dfdx)[0]
        zeta_r_list.append(zeta_r)
    zeta_r_lists.append(np.array(zeta_r_list))
    
#creating empty lists to store zeta_modes,frequency_modes,velocity_modesa and wavelength_modes
zeta_modes=[]
f_modes=[]
vel_modes=[]
wav_l_modes=[]
for m in range(4):
    zeta_mode=[]
    f_mode=[]

    for k in range(m,10):
        zeta_mode.append(zeta_r_lists[k][m])
        f_mode.append(f_list[k])

    zeta_mode=np.array(zeta_mode)
    f_mode=np.array(f_mode)
    vel_mode=1/np.sqrt(B1**-2-zeta_mode**2/H**2)
    wav_l_mode=vel_mode/f_mode
    vel_modes.append(vel_mode)
    wav_l_modes.append(wav_l_mode)
    


    zeta_modes.append(zeta_mode)
    f_modes.append(f_mode)



    

print(f_modes)
plt.figure(figsize=(6,8))
plt.subplot(3,1,1)
for f,z in zip(f_modes,zeta_modes):
    plt.plot(f,z)
plt.xlabel('f[Hz]')
plt.ylabel('zeta[s]')
plt.legend([f'mode {k}' for k in range(4)])



plt.subplot(3,1,2)
for f,v in zip(f_modes,vel_modes):
    plt.plot(f,v)
plt.xlabel('f[Hz]')
plt.ylabel('C_L[m/s]')
plt.legend([f'mode {k}' for k in range(4)])

plt.subplot(3,1,3)
for f,wavl in zip(f_modes,wav_l_modes):
    plt.plot(f,wavl)
plt.xlabel('f[Hz]')
plt.ylabel('lamda_l[m]')
plt.legend([f'mode {k}' for k in range(4)])



plt.show()
plt.savefig('modes.png')


