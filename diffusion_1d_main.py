import numpy as np
import matplotlib.pyplot as plt
import diffusion_1d as d2d

X = 0.01
Nx = 100
T = 0.1*10
Nt = 2000
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
    #print np.amax(u[int(i*(Nt-1)/11)])-np.amin(u[int(i*(Nt-1)/11)])

#plt.show()

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
#a = np.zeros(Nt,float)

print x_grid-X/2+dx/2

for i in range(Nt-1):
    chi = 1.72e-5*np.log(u[i]+293.0) + 1.386e-4
    B = Bd*(z[i] + np.sin(theta[i])*(x_grid-X/2+dx/2))
    
    #print B
    f = sum(chi*B)*X*X*Bd/mu0
    #print f
    #print rho*g*h*X*X
    a = (rho*g*h*X*X-f*np.cos(theta[i]))/m
    #print a
    z[i+1] = vz*dt + 0.5*a*dt*dt + z[i]
    vz = a*dt
    #print z[i+1]

    tao = sum(chi*B*(x_grid-X/2+dx/2))*X*X*Bd/mu0
    #print tao
    #print I
    a_omega = tao/I
    theta[i+1] = omega*dt + 0.5*a_omega*dt*dt + theta[i]
    omega = a_omega*dt

    #print B
    fx = sum(chi*B)*X*X*Bd/mu0*np.sin(theta[i])
    #print f
    #print rho*g*h*X*X
    ax = fx/m
    #print a
    x[i+1] = vx*dt + 0.5*ax*dt*dt + x[i]
    vx = ax*dt
    #print z[i+1]
    

plt.figure()
plt.plot(range(Nt),z)
plt.figure()
plt.plot(range(Nt),theta)
plt.ylabel(r'$\theta$')
plt.figure()
plt.plot(range(Nt),x)
plt.ylabel(r'x')


plt.show()
    
    
    
