import numpy as np
import matplotlib.pyplot as plt
import diffusion_1d as d2d

X = 0.01
Nx = 100
T = 1.0
Nt = 200
D = 1.22e-3
q0 = 6.78e2

u = d2d.heat_diffusion(X,T,Nx,Nt,D,q0)

print np.amax(u)
print np.amin(u)

dx = X/Nx
dt = T/Nt
x_grid = np.array(range(Nx))*dx

fig = plt.figure()
for i in range(12):
    plt.subplot(3,4,i+1)
    plt.plot(x_grid,u[int(i*(Nt-1)/11)])
    fig.text(0.5, 0.04, r'$x/m$', ha='center', va='center')
    fig.text(0.06, 0.5, r'$T-293K/K$', ha='center', va='center', rotation='vertical')

#plt.show()

def maglevmotion(u):   
    Bd = 0.01
    rho = 2.27e3
    g = 10
    h = 1e-4
    m = rho*h*X**2
    mu0 = 1e-7
    chi0 = 1.72e-5*np.log(293.0) + 1.386e-4
    z0 = rho*g*h*mu0/(chi0*Bd*Bd)
    B0 = Bd*z0
    B = B0
    vz = 0.0
    vx = 0.0
    omega = 0.0
    I = m*X*X/12

    z = z0*np.ones(Nt,float)
    theta = np.zeros(Nt,float)
    x = np.zeros(Nt,float)

    for i in range(Nt-1):
        chi = -1.72e-5*np.log(u[i]+293.0) + 1.386e-4
        B = Bd*(z[i] + np.sin(theta[i])*(x_grid-X/2+dx/2))
        
        f = sum(chi*B)*X*X*Bd/mu0
        a = (rho*g*h*X*X-f*np.cos(theta[i]))/m
        z[i+1] = vz*dt + 0.5*a*dt*dt + z[i]
        vz = a*dt

        tau = sum(chi*B*(x_grid-X/2+dx/2))*X*X*Bd/mu0
        a_omega = tau/I
        theta[i+1] = omega*dt + 0.5*a_omega*dt*dt + theta[i]
        omega = a_omega*dt

        fx = sum(chi*B)*X*X*Bd/mu0*np.sin(theta[i])
        ax = fx/m
        x[i+1] = vx*dt + 0.5*ax*dt*dt + x[i]
        vx = ax*dt

    return x

x = maglevmotion(u)
    
plt.figure()
plt.plot(dt*np.array(range(Nt)),x)
plt.xlabel(r't/s')
plt.ylabel(r'x/m')

plt.show()
