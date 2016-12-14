import numpy as np
import matplotlib.pyplot as plt

def update_ghost(u,ghost,Nx,Ny):  # update the Neumann BC ghost zones
    u_upd = np.copy(u)
    for i in range(3,Nx-3):
        for j in range(3,Ny-3):
            if ghost[i,j] > 0:  # find the ghost zones and neighbors
                u_upd[i,j] = (u[i-1,j]+u[i,j-1]+u[i+1,j]+u[i,j+1])/ghost[i,j]
    return u_upd

def heat_diffusion(X,Y,T,Nx,Ny,Nt,D,q0):
    # grid of the disk
    dx = X/Nx
    dy = Y/Ny
    dt = T/Nt
    
    # center of the disk
    ctr = float(Nx)/2
    r = ctr*0.8

    # center of the light spot
    ctr_lig = float(Nx)/4
    r_lig = 0.05*ctr

    # use matrices to get the shape and positions of the disk and ghost zones
    disk = np.zeros([Nx,Ny])
    ghost1 = np.zeros([Nx,Ny])
    ghost2 = np.zeros([Nx,Ny])

    # u is temperature and q is the source function
    u = np.zeros([Nx,Ny,Nt],float)
    q = np.zeros([Nx,Ny],float)
    for i in range(Nx-1):
        for j in range(Ny-1):
            if ((float(i)-ctr_lig)**2+(float(j)-ctr_lig)**2) < r_lig**2:
                q[i,j] = q0
    mu = 5.67e-10
                
    # obtain ghost zones around the disk along x and y directions
    for i in range(Nx-1):
        find = 0
        for j in range(Ny-1):
            if ((float(i)-ctr)**2+(float(j)-ctr)**2) < r**2:
                disk[i,j] = 1.0
                if find == 0:
                    ghost1[i,j-1] = 1
                find = 1
            else:
                if find == 1:
                    ghost1[i,j] = 1
                find = 0

    for j in range(Ny-1):
        find = 0
        for i in range(Nx-1):
            if ((float(i)-ctr)**2+(float(j)-ctr)**2) < r**2:
                disk[i,j] = 1.0
                if find == 0:
                    ghost2[i-1,j] = 1
                find = 1
            else:
                if find == 1:
                    ghost2[i,j] = 1
                find = 0
    # add two ghost zones of x and y direction up
    ghost = ghost1 + ghost2

    # do the time evolution by explicit iterations
    for i in range(1,Nt):
        u[1:Nx-1,1:Ny-1,i] = (u[1:Nx-1,1:Ny-1,i-1] +\
                             dt*D*((u[2:Nx,1:Ny-1,i-1] -\
                             2*u[1:Nx-1,1:Ny-1,i-1] +\
                             u[0:Nx-2,1:Ny-1,i-1])/(dx*dx) +\
                             (u[1:Nx-1,2:Ny,i-1]-2*u[1:Nx-1,1:Ny-1,i-1] +\
                             u[1:Nx-1,0:Ny-2,i-1])/(dy*dy)) +\
                             -dt*mu*((u[1:Nx-1,1:Ny-1,i-1]+293.0)**4-293.0**4) +\
                             dt*q[1:Nx-1,1:Ny-1])*disk[1:Nx-1,1:Ny-1];
        u[1:Nx-1,1:Ny-1,i] = update_ghost(u[1:Nx-1,1:Ny-1,i],ghost[1:Nx-1,1:Ny-1],Nx,Ny)

    return u

    
