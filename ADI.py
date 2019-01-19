import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plotIt(x, y, v, n, saveIt=False):
    Y, X = np.meshgrid(y, x)
    fig = plt.figure(figsize=(12, 8))
    imAx = fig.add_subplot(121)
    imAx.imshow(v, cmap='Blues', vmin=0., vmax=5., interpolation='none')
    threeDAx = fig.add_subplot(122, projection='3d')
    threeDAx.plot_surface(X, Y, v, cmap='Reds', edgecolor='none',
                          cstride=1, rstride=1, vmin=0., vmax=1.)
    threeDAx.set_zlim(0., 5.)
    fig.suptitle('n = '+str(n))
    if saveIt:
        s = str(n)
        for i in range(5-len(s)):
            s = '0'+s
        s += '.png'
        s = 'plotIt_n'+s
        print('Saving: '+s)
        plt.savefig(s)
        plt.close()
    else:
        return fig
    #this function plots the graph and saves it if you choose
    #x and y are spatial dimensions, v is the array of values for each grid point (temp, density, etc.)
    #n is the time step

def thomasAlgorithm(a, b, c, d):
    y = np.empty(a.shape)
    x = y.copy()
    z = 1./b[0]
    y[0] = z*c[0]
    x[0] = z*d[0]
    for i in range(1, a.size):
        z = 1./(b[i]-a[i]*y[i-1])
        y[i] = z*c[i]
        x[i] = z*(d[i]+a[i]*x[i-1])
    li = range(a.size-1)
    li.reverse()
    for i in li:
        x[i] += y[i]*x[i+1]
    return x
    #algorithm for solving tri diagonal matrices

def cyclicThomasAlgorithm(a, b, c, d):
    alpha, beta, gamma = -c[-1], a[0], b[0]
    aa = a.copy()
    aa[0] = 0.
    bb = b.copy()
    bb[0]  -= gamma
    bb[-1] -= alpha*beta/gamma
    cc = c.copy()
    cc[-1] = 0.
    x = thomasAlgorithm(aa, bb, cc, d)
    e = np.zeros(d.size)
    e[0]  = gamma
    e[-1] = alpha
    y = thomasAlgorithm(aa, bb, cc, e)
    f = (x[0]+beta*x[-1]/gamma)/(1.0+y[0]+beta*y[-1]/gamma)
    x -= f*y
    return x
    #same but cyclic

def capped(z):
    try:
        out = np.zeros(z.size)
    except:
        out = 0.
    return out
    #for the Nuemann boundary conditions

def lump(x, y):
    return 4*np.cos(0.5*np.pi*x)*np.cos(0.5*np.pi*y)/1.62197252869
    #makes a smooth central peak
    
def slosh(x, y):
    if x > 0:
        return 0
    else:
        return 4*np.sqrt(np.sin(0.5*np.pi*x)**2)*np.cos(0.5*np.pi*y)/0.810986264347
    #makes a peak at the center of one edge
        
    #above are initial conditions for v
    #weird denominators make the total mass an integer (more or less).

def twoD_ADI_NeumannBoundaries(Nx=40, Ny=40, Nt=1500,
                               theta=0.5, alpha=0.1, deltat=0.001,
                               xmin=-1., xmax=1., ymin=-1., ymax=1.,
                               initialConditions=slosh,
                               lowerXBoundaryDerivative=capped,
                               lowerYBoundaryDerivative=capped,
                               upperXBoundaryDerivative=capped,
                               upperYBoundaryDerivative=capped,
                               inc=100, saveIt=False):
    """
           D
        -------
       |       |
     A |       | C
       |       |
        -------
           B
    """
    def getDx():
        if j==0:
            dx =  (1.0-(1.0-theta)*muy)  *v[:, 0]\
                 +(1.0-theta)*muy        *v[:, 1]\
                 -(1.0-theta)*muy*deltay *B
        elif j==Ny-1:
            dx =  (1.0-theta)*muy        *v[:, -2]\
                 +(1.0-(1.0-theta)*muy)  *v[:, -1]\
                 +(1.0-theta)*muy*deltay *D
        else:
            dx =  (1.0-theta)*muy           *v[:, j-1]\
                 +(1.0-2.0*(1.0-theta)*muy) *v[:, j  ]\
                 +(1.0-theta)*muy           *v[:, j+1]
        dx[0]  -= theta*mux*deltax *A[j]
        dx[-1] += theta*mux*deltax *C[j]
        return dx
    
    def getDy():
        if i==0:
            dy =  (1.0-(1.0-theta)*mux)  *v[0, :]\
                 +(1.0-theta)*mux        *v[1, :]\
                 -(1.0-theta)*mux*deltax *A
        elif i==Nx-1:
            dy =  (1.0-theta)*mux        *v[-2, :]\
                 +(1.0-(1.0-theta)*mux)  *v[-1, :]\
                 +(1.0-theta)*mux*deltax *C
        else:
            dy =  (1.0-theta)*mux           *v[i-1, :]\
                 +(1.0-2.0*(1.0-theta)*mux) *v[i  , :]\
                 +(1.0-theta)*mux           *v[i+1, :]
        dy[0]  -= theta*muy*deltay *B[i]
        dy[-1] += theta*muy*deltay *D[i]
        return dy
        #magic maths stuff
    
    x = np.linspace(xmin, xmax, Nx+1)[:-1]
    deltax = np.diff(x)[0]
    x += 0.5*deltax
    #x grid
    
    y = np.linspace(ymin, ymax, Ny+1)[:-1]
    deltay = np.diff(y)[0]
    y += 0.5*deltay
    #y grid
    
    #t = np.arange(0, Nt, deltat)
    #time array, not actually neccessary
    
    mux = alpha*deltat/(deltax**2.)
    muy = alpha*deltat/(deltay**2.)
    
    ax = np.ones(Nx)*theta*mux
    ax[0] = 0.
    bx = np.ones(Nx)*(1.+theta*mux)
    bx[1:-1] += theta*mux
    cx = np.ones(Nx)*theta*mux
    cx[-1] = 0.
    
    ay = np.ones(Ny)*theta*muy
    ay[0] = 0.
    by = np.ones(Ny)*(1.+theta*muy)
    by[1:-1] += theta*muy
    cy = np.ones(Ny)*theta*muy
    cy[-1] = 0.
    #defining terms for the Thomas algorithm
    
    A = lowerXBoundaryDerivative(y)
    B = lowerYBoundaryDerivative(x)
    C = upperXBoundaryDerivative(y)
    D = upperYBoundaryDerivative(x)
    #Neumann boundaries
    
    v = np.empty([Nx, Ny])
    for i in range(Nx):
        for j in range(Ny):
            v[i, j] = initialConditions(x[i], y[j])
    #v is the parameter to vary i.e. density, which is initialised here at t = 0
    
    M = deltax*deltay*np.sum(v)
    print 'Total mass =', M
    #finds the total mass if v is density by approximating volume integral    
    
    vstar = v.copy()
    tolerance = 10**-5
    #the maximum variation in v between time steps
    #which indicates a steady state has been reached
    
    for n in range(Nt-1):
        
        for j in range(Ny):
            vstar[:, j] = thomasAlgorithm(ax, bx, cx, getDx())
            
        if np.all(v - vstar < tolerance):
            print 'Steady state reached at n =', n
            break
        #stops the loop when a steady state is reached for a given tolerance
                  
        v = vstar.copy()
        
        for i in range(Nx):
            vstar[i, :] = thomasAlgorithm(ay, by, cy, getDy())
            
        if np.all(v - vstar < tolerance):
            print 'Steady state reached at n =', n
            break
        
        v = vstar.copy()
        #finds v at each time step
              
        
        if n%inc==0:
            plotIt(x, y, v, n, saveIt=saveIt)
            #plots v at each time step which is a multiple of inc
            
        
    
    #plotIt(x, y, v, Nt, saveIt=saveIt)
    #plots the final time step. not neccessary.
    
    if not saveIt:
        plotIt(x, y, v, Nt)
    plt.show()
    #shows the plot if you choose not to save it

if __name__=='__main__':
    import os
    
    twoD_ADI_NeumannBoundaries(Nt=15000, inc=100, saveIt=True)
    #runs the simulation
    
    cwd = os.getcwd()
    dest = os.path.join(cwd, 'ImagesFromADI')
    for fileName in os.listdir(cwd):
        if 'plotIt' in fileName:
            os.rename(os.path.join(cwd, fileName),
                      os.path.join(dest, fileName))
    #saves it to the same folder as the .py file, if saveIt is true