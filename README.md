def f(x):
    f=(p2/p1)((H^2(B1^-2-B2^-2)-x^2)^-1/2)/x




def secant_method(x1,x2,f,dfdx,eps_s=1*10^-8):
"""finds the root of the equation
    x1,x2 : the interval used while searching for the root;
            f(x1)*f(x2)<0   
    f     : the function we are solving

"""
    eps_a = 2*eps_s
    dx=xk-x2
    while eps_a > eps_s:
        fx1=f(x1)
        fx2=f(x2)
        xk=x2-((x1-x2)*fx1)/(fx1-fx2)
        eps_a = np.abs(dx/xk)
    return xk, k, eps_a