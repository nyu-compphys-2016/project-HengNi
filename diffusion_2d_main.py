import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import diffusion_2d as d2d

X = 0.01
Y = 0.01
Nx = 40
Ny = 40
T = 0.1
Nt = 20000
D = 1.22e-3
q0 = 6.78e5

u = d2d.heat_diffusion(X,Y,T,Nx,Ny,Nt,D,q0)

print np.amax(u)
print np.amin(u)

fig, axes = plt.subplots(nrows=3, ncols=4)
i = 0
for ax in axes.flat:
    im = ax.imshow(u[:,:,round(i*(Nt-1)/11)],vmin=np.amin(u),vmax=np.amax(u),cmap='jet')
    print np.amax(u[:,:,round(i*(Nt-1)/11)])-u[10,10,round(i*(Nt-1)/11)]
    i = i+1
    ax.axis('off')
    
cax,kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
plt.colorbar(im, cax=cax, **kw)
plt.show()
