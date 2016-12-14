import numpy as np
import matplotlib.pyplot as plt

def TDMAsolver(a,b,c,d):    # solver of tridiagnal matrix equations
    nf = len(d)
    ac,bc,cc,dc = map(np.array,(a,b,c,d))
    for it in xrange(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
        	    
    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in xrange(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc

def heat_diffusion(X,T,Nx,Nt,D,q0):
    # 1d grid
    dx = X/Nx
    dt = T/Nt

    x_grid = np.array(range(Nx))*dx
    t_grid = np.array(range(Nt))*dt

    alpha = D*dt/(2*dx*dx)
    mu = 5.67e-10

    u = np.zeros(Nx,float)
    q = np.zeros(Nx,float)
    q[10:20] = q0

    
    # matrix for Crank Nicolson iterations
    A = np.diagflat([-alpha for i in range(Nx-1)], -1) +\
          np.diagflat([1.+alpha]+[1.+2.*alpha for i in range(Nx-2)]+[1.+alpha]) +\
          np.diagflat([-alpha for i in range(Nx-1)], 1)
            
    B = np.diagflat([alpha for i in range(Nx-1)], -1) +\
          np.diagflat([1.-alpha]+[1.-2.*alpha for i in range(Nx-2)]+[1.-alpha]) +\
          np.diagflat([alpha for i in range(Nx-1)], 1)

    u_tevol = []
    u_tevol.append(u)

    a_u = [-alpha for i in range(Nx-1)]
    b_u = [1.+alpha]+[1.+2.*alpha for i in range(Nx-2)]+[1.+alpha]
    c_u = [-alpha for i in range(Nx-1)]

    fig = plt.figure()
    
    for ti in range(1,Nt):
        d_u = np.dot(B,u) + dt*q - dt*mu*((u-293.0)**4-293.0**4)
        u = TDMAsolver(a_u, b_u, c_u, d_u)
        u_tevol.append(u)
        #plt.plot(x_grid,u,'b',alpha = float(ti)/Nt)
        #fig.show()

    return u_tevol
