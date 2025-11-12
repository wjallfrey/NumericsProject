import numpy as np
import matplotlib.pyplot as plt

def modA2(b, c, kdx): #modAsquared in Von Neumann analysis (linearised)
    return (1 + b * (np.cos(kdx) - 1)) ** 2 + (c * np.sin(kdx)) ** 2

def maxmodA2(b,c): #calculates the maximum value of modA2 over all kdx
    kdx = np.arange(1, 314) / 100
    return np.max(modA2(b, c, kdx))

def FTCS(nu,T,NT,N,plotting,stabilityTest,multiplot,massplot):
    mass=[]
    U = 1
    L = np.pi
    deltat=T/NT
    xgrid = np.linspace(0, L, N + 1)    
    deltax = L / N
    D = 2 * nu * deltat / deltax ** 2
    C = U * deltat / deltax
    initialphi = np.sin(2 * xgrid) ** 2
    currentphi = initialphi
    energy0=deltax*sum(currentphi**2)
    energy=[energy0]
    for k in range(1, NT + 1):
        newphi = np.zeros_like(currentphi)
        newphi[0] = currentphi[0] + nu * deltat / deltax ** 2 * (currentphi[1] - 2 * currentphi[0] + currentphi[N]) - 0.5*deltat / deltax * currentphi[0] * (currentphi[1] - currentphi[N])
        newphi[N] = currentphi[N] + nu * deltat / deltax ** 2 * (currentphi[0] - 2 * currentphi[N] + currentphi[N - 1]) - 0.5*deltat / deltax * currentphi[N] * (currentphi[0] - currentphi[N - 1])
        for j in range(1, N):
            newphi[j] = currentphi[j] + nu * deltat / deltax ** 2 * (currentphi[j + 1] - 2 * currentphi[j] + currentphi[j - 1]) - 0.5*deltat / deltax * currentphi[j] * (currentphi[j+1] - currentphi[j - 1])
        currentphi = newphi
        mass.append(deltax*sum(currentphi))
        energy.append(deltax*sum(currentphi**2))
        if stabilityTest and max(energy)>energy0*1.1: #if energy has increased then exit the loop and record instability
            return(C,D,True)
            break        
        if multiplot and k % (NT // 10) == 0: #plot the solution at 10 evenly spaced times.
            plt.plot(xgrid, newphi, color='blue')
        if massplot:
            if sum(currentphi**2)>100000000000000:
                break
    if plotting:
        plt.figure()
        plt.plot(xgrid,currentphi)
    if massplot:
        plt.figure()
        plt.plot(np.linspace(0,T,k),mass) #plot of mass over time
        plt.xlabel("time")
        plt.ylabel("mass")
        #plt.title("Plot of solution at t="+str(T)+" with varying spatial resolution, keeping c fixed")
        if stabilityTest:
            plt.figure()
            plt.plot(np.linspace(0,T,k+1),energy) # plot of energy over time
            plt.xlabel("time")
            plt.ylabel("energy")
    if stabilityTest:
        return(C,D,False)
    else:
        return(currentphi)


def FFT(nu,T,NT,N,plotting): #FFT scheme
    L = np.pi
    deltat=T/NT
    xgrid = np.linspace(0, L, N + 1)
    initialphi = np.sin(2 * xgrid) ** 2
    deltax = L / N
    mass=[]
    currentphifft=initialphi
    k=2*np.concatenate((np.linspace(0,0,1),np.linspace(1,N//2,N//2),np.linspace(-N//2,-1,N//2)))
    ik=1j*k
    k2=k**2
    for i in range(1, NT + 1): # for each time step
        phihat=np.fft.fft(currentphifft)
        phisquaredhat=np.fft.fft(currentphifft**2)
        newphihat=phihat-nu*deltat*k2*phihat-deltat/2*ik*phisquaredhat
        newphi=np.fft.ifft(newphihat).real #Calculate phi at next timestep - note imaginary parts small
        currentphifft=newphi
        mass.append(deltax*sum(currentphifft))
        if plotting and i % (NT // 10) == 0: #plot the solution at 10 evenly spaced times.
            plt.plot(xgrid, newphi, color='blue')
    if plotting:
        plt.figure()
        plt.plot(xgrid, initialphi, color='blue')
        plt.plot(xgrid,currentphifft)
        plt.xlabel("x")
        plt.ylabel("phi")
        plt.figure()
        plt.plot(np.linspace(0,T,NT),mass)
        plt.xlabel("time")
        plt.ylabel("mass")
    return(currentphifft)

def VNstability(resolution): #Calculates the value of modA^2 for a grid of points of d and c. resolution 400 is good
    BN = int(1.25*resolution)
    CN = int(1.25*resolution)
    bvec = np.arange(BN + 1) / resolution
    cvec = np.arange(CN + 1) / resolution
    heatmat = np.zeros((BN + 1, CN + 1))
    plt.figure()
    for i in range(BN + 1):
        for j in range(CN + 1):
            heatmat[i, j] = maxmodA2(bvec[i], cvec[j])
            #heatmat[i, j] = modA2(bvec[i], cvec[j],np.pi/2)
    plt.contour(bvec, cvec, heatmat.T, levels=np.arange(0, 2), cmap='Greys')
    plt.xlabel('d')
    plt.ylabel('c')
    plt.title('Plot of max|A|^2=1')
    plt.colorbar()

def interpolateanalytic(analyticSol,maxN,N): #interpolates analyticSol from an xgrid of maxN to a courser xgrid of N
    interpolatedSol=np.linspace(0,0,N+1)
    for i in range(0,N):
        interpolatedSol[i]=analyticSol[int((i*maxN/N)//1)]*(1-(i*maxN/N)%1)+analyticSol[int((i*maxN/N)//1)+1]*((i*maxN/N)%1)
    interpolatedSol[N]=analyticSol[maxN]
    return(interpolatedSol)

def analyticSol(nu,T,NT,maxN):
    return(FFT(nu,T,NT,maxN,False))

def errors(nu,T,maxN,invc,L,analyticSol,Nrange,scheme,plotting): #gives the errors between interpolated analyticSol and the numerical method using scheme (FFT or FTCS) over a range of NX 
    errorsl1=[]
    errorsl2=[]
    errorslinf=[]
    if plotting:
        plt.figure()
    for N in Nrange:
        N=int(N)
        if scheme=="FFT":
            t1solution=FFT(nu,T,invc*N,N,False)
        elif scheme=="FTCS":
            t1solution=FTCS(nu,T,invc*N,N,False,False,False,False)
        else:
            break
        if plotting:
            plt.plot(np.linspace(0,L,N+1),t1solution)
        interpolatedSol=interpolateanalytic(analyticSol,maxN,N)
        errorsl1.append(sum(abs(t1solution-interpolatedSol))/sum(abs(interpolatedSol)))
        errorsl2.append((sum((t1solution-interpolatedSol)**2)/sum((interpolatedSol)**2))**(1/2))
        errorslinf.append(max(abs(t1solution-interpolatedSol))/max(abs(interpolatedSol)))
    deltaxs=L*Nrange**(-1)
    if plotting:
        plt.figure()
        plt.plot(deltaxs,errorsl1,'-x')
        plt.plot(deltaxs,errorsl2,'-x')
        plt.xlabel("Delta x")
        plt.ylabel("Error")
        plt.title("plot of errors (l1 and l2) against grid scale deltax")
        plt.xlim(0,max(deltaxs)*1.05)
        plt.figure()
        plt.plot(np.log(deltaxs),np.log(errorsl1),'-x')
        plt.plot(np.log(deltaxs),np.log(errorsl2),'-x')
        plt.plot(np.log(deltaxs),np.log(errorsl2[-1])+2*(np.log(deltaxs)-np.log(deltaxs[-1])))
        plt.xlabel("log(Delta x)")
        plt.ylabel("log(Error)")
        plt.title("log plot of errors (l1 and l2) and grid scale. Line of gradient 2 shown in green")

    return(errorsl1,errorsl2,errorslinf,deltaxs)

def stability(resolution): #checks for stability via energy increase over a grid of c and d values. 0<c<1.4, 0<d<1.2.
    N=100
    Kc=int(1.4*resolution)
    Kd=int(1.2*resolution)
    outputs=np.zeros([Kc,Kd])
    testcrange=np.linspace(0.001,1.4,Kc)
    for i in range(Kc):
        testc=testcrange[i]
        nurange=np.linspace(0,1.2*np.pi/(2*N*testc),Kd)
        T=10*(np.pi*testc)
        NT=round(N*T/(np.pi*testc)) #fixed at 1000
        for j in range(Kd):
            outputs[i,j]=FTCS(nurange[j],T,NT,N,False,True,False,False)[2]
    plt.figure()
    plt.contour(np.linspace(0,1.2,Kd),testcrange,outputs,levels=np.arange(0,2))
    plt.plot(np.linspace(0.001,1,30)**2,np.linspace(0.001,1,30))
    plt.title("plot of the boundary of stability in the c,d parameter space - stable inside the region")

#plots: 
FTCS(0.1,1,1000,100,False,False,True,False) #F1
FTCS(0.004,0,100,150,True,False,False,False) #F4
FTCS(0.004,1,100,150,True,False,False,False)
FTCS(0.004,2,200,150,True,False,False,False)
FTCS(0.004,3,300,150,True,False,False,False)
FTCS(0.004,7,700,150,True,False,False,False)
FTCS(0.004,10,1000,150,True,False,False,False)

VNstability(400) #F2
stability(25) #F3 NOTE: For the plot in report resolution of 100 was used. Here we use 25 to save time.

analytic=analyticSol(0.1,1,100000,800)
errors(0.1,1,800,30,np.pi,analytic,np.linspace(100,800,8),"FTCS",True) #F5,6,7

FTCS(0.1,1,1000,100,False,False,False,True) #massplot1 F8
FTCS(0.1,1,1000,400,False,False,False,True) #massplot2 F9
plt.show()  
