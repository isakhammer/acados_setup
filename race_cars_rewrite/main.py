#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#

# author: Daniel Kloeser

import time, os
import numpy as np
from acados_settings import *
from pathlib import Path
from matplotlib import cm
import matplotlib.pyplot as plt
import centerline

def transformProj2Orig(si,ni,alpha,v,filename='LMS_Track.txt'):
    [sref,xref,yref,psiref,_]=getTrack(filename=filename)
    tracklength=sref[-1]
    si=si%tracklength
    idxmindist=findClosestS(si,sref)
    idxmindist2=findSecondClosestS(si,sref,idxmindist)
    t=(si-sref[idxmindist])/(sref[idxmindist2]-sref[idxmindist])
    x0=(1-t)*xref[idxmindist]+t*xref[idxmindist2]
    y0=(1-t)*yref[idxmindist]+t*yref[idxmindist2]
    psi0=(1-t)*psiref[idxmindist]+t*psiref[idxmindist2]

    x=x0-ni*np.sin(psi0)
    y=y0+ni*np.cos(psi0)
    psi=psi0+alpha
    v=v
    return x,y,psi,v


def findClosestS(si,sref):
    # Get number of elements
    if(np.isscalar(si)):
        N=1
    else:
        N=np.array(si).shape[0]
    mindist=100000*np.ones(N)
    idxmindist=np.zeros(N)
    for i in range(sref.size):
        di=abs(si-sref[i])
        idxmindist = np.where(di < mindist,i, idxmindist)
        mindist = np.where(di < mindist, di, mindist)
    idxmindist = np.where(idxmindist==sref.size,1,idxmindist)
    idxmindist = np.where(idxmindist<1,sref.size-1,idxmindist)
    return idxmindist.astype(int)


def findSecondClosestS(si,sref,idxmindist):
    d1=abs(si-sref[idxmindist-1])               # distance to node before
    d2=abs(si-sref[(idxmindist+1)%sref.size])   # distance to node after
    idxmindist2 = np.where(d1>d2,idxmindist+1,idxmindist-1) # decide which node is closer
    idxmindist2 = np.where(idxmindist2==sref.size,0,idxmindist2)    # if chosen node is too large
    idxmindist2 = np.where(idxmindist2<0,sref.size-1,idxmindist2)   # if chosen node is too small

    return idxmindist2

def transformOrig2Proj(x,y,psi,v,filename='LMS_Track.txt'):
    [sref,xref,yref,psiref,_]=getTrack(filename=filename)
    idxmindist=findClosestPoint(x,y,xref,yref)
    idxmindist2=findClosestNeighbour(x,y,xref,yref,idxmindist)
    t=findProjection(x,y,xref,yref,sref,idxmindist,idxmindist2)
    s0=(1-t)*sref[idxmindist]+t*sref[idxmindist2]
    x0=(1-t)*xref[idxmindist]+t*xref[idxmindist2]
    y0=(1-t)*yref[idxmindist]+t*yref[idxmindist2]
    psi0=(1-t)*psiref[idxmindist]+t*psiref[idxmindist2]

    s=s0
    n=np.cos(psi0)*(y-y0)-np.sin(psi0)*(x-x0)
    alpha=psi-psi0
    v=v
    return s,n,alpha,v

def findProjection(x,y,xref,yref,sref,idxmindist,idxmindist2):
    vabs=abs(sref[idxmindist]-sref[idxmindist2])
    vl=np.empty(2)
    u=np.empty(2)
    vl[0]=xref[idxmindist2]-xref[idxmindist]
    vl[1]=yref[idxmindist2]-yref[idxmindist]
    u[0]=x-xref[idxmindist]
    u[1]=y-yref[idxmindist]
    t=(vl[0]*u[0]+vl[1]*u[1])/vabs/vabs
    return t

def findClosestPoint(x,y,xref,yref):
    mindist=1
    idxmindist=0
    for i in range(xref.size):
        dist=dist2D(x,xref[i],y,yref[i])
        if dist<mindist:
            mindist=dist
            idxmindist=i
    return idxmindist

def findClosestNeighbour(x,y,xref,yref,idxmindist):
    distBefore=dist2D(x,xref[idxmindist-1],y,yref[idxmindist-1])
    distAfter=dist2D(x,xref[idxmindist+1],y,yref[idxmindist+1])
    if(distBefore<distAfter):
        idxmindist2=idxmindist-1
    else:
        idxmindist2=idxmindist+1
    if(idxmindist2<0):
        idxmindist2=xref.size-1
    elif(idxmindist==xref.size):
        idxmindist2=0
    return idxmindist2


def dist2D(x1,x2,y1,y2):
    return np.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))



def plotTrackProj(simX,filename='LMS_Track.txt', T_opt=None):
    # load track
    s=simX[:,0]
    n=simX[:,1]
    alpha=simX[:,2]
    v=simX[:,3]
    distance=0.12
    # transform data
    [x, y, _, _] = transformProj2Orig(s, n, alpha, v,filename)
    # plot racetrack map

    #Setup plot
    plt.figure()
    plt.ylim(bottom=-1.75,top=0.35)
    plt.xlim(left=-1.1,right=1.6)
    plt.ylabel('y[m]')
    plt.xlabel('x[m]')

    # Plot center line
    [Sref,Xref,Yref,Psiref,_]=getTrack(filename)
    plt.plot(Xref,Yref,'--',color='k')

    # Draw Trackboundaries
    Xboundleft=Xref-distance*np.sin(Psiref)
    Yboundleft=Yref+distance*np.cos(Psiref)
    Xboundright=Xref+distance*np.sin(Psiref)
    Yboundright=Yref-distance*np.cos(Psiref)
    plt.plot(Xboundleft,Yboundleft,color='k',linewidth=1)
    plt.plot(Xboundright,Yboundright,color='k',linewidth=1)
    plt.plot(x,y, '-b')

    # Draw driven trajectory
    heatmap = plt.scatter(x,y, c=v, cmap=cm.rainbow, edgecolor='none', marker='o')
    cbar = plt.colorbar(heatmap, fraction=0.035)
    cbar.set_label("velocity in [m/s]")
    ax = plt.gca()
    ax.set_aspect('equal', 'box')

    # Put markers for s values
    xi=np.zeros(9)
    yi=np.zeros(9)
    xi1=np.zeros(9)
    yi1=np.zeros(9)
    xi2=np.zeros(9)
    yi2=np.zeros(9)
    for i in range(int(Sref[-1]) + 1):
        try:
            k = list(Sref).index(i + min(abs(Sref - i)))
        except:
            k = list(Sref).index(i - min(abs(Sref - i)))
        [_,nrefi,_,_]=transformOrig2Proj(Xref[k],Yref[k],Psiref[k],0)
        [xi[i],yi[i],_,_]=transformProj2Orig(Sref[k],nrefi+0.24,0,0)
        # plt.text(xi[i], yi[i], f'{i}m', fontsize=12,horizontalalignment='center',verticalalignment='center')
        plt.text(xi[i], yi[i], '{}m'.format(i), fontsize=12,horizontalalignment='center',verticalalignment='center')
        [xi1[i],yi1[i],_,_]=transformProj2Orig(Sref[k],nrefi+0.12,0,0)
        [xi2[i],yi2[i],_,_]=transformProj2Orig(Sref[k],nrefi+0.15,0,0)
        plt.plot([xi1[i],xi2[i]],[yi1[i],yi2[i]],color='black')

def plotRes(simX,simU,t):
    # plot results
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.step(t, simU[:,0], color='r')
    plt.step(t, simU[:,1], color='g')
    plt.title('closed-loop simulation')
    plt.legend(['dD','ddelta'])
    plt.ylabel('u')
    plt.xlabel('t')
    plt.grid(True)
    plt.subplot(2, 1, 2)
    plt.plot(t, simX[:,:])
    plt.ylabel('x')
    plt.xlabel('t')
    plt.legend(['s','n','alpha','v','D','delta'])
    plt.grid(True)

def plotalat(simX,simU,constraint,t):
    Nsim=t.shape[0]
    plt.figure()
    alat=np.zeros(Nsim)
    for i in range(Nsim):
        alat[i]=constraint.alat(simX[i,:],simU[i,:])
    plt.plot(t,alat)
    plt.plot([t[0],t[-1]],[constraint.alat_min, constraint.alat_min],'k--')
    plt.plot([t[0],t[-1]],[constraint.alat_max, constraint.alat_max],'k--')
    plt.legend(['alat','alat_min/max'])
    plt.xlabel('t')
    plt.ylabel('alat[m/s^2]')

def getTrack(filename):
    track_file = os.path.join(str(Path(__file__).parent),"tracks", filename)
    array=np.loadtxt(track_file)
    sref=array[:,0]
    xref=array[:,1]
    yref=array[:,2]
    psiref=array[:,3]
    kapparef=array[:,4]
    return sref,xref,yref,psiref,kapparef

"""
Example of the frc_racecars in simulation without obstacle avoidance:
This example is for the optimal racing of the frc race cars. The model is a simple bicycle model and the lateral acceleration is constraint in order to validate the model assumptions.
The simulation starts at s=-2m until one round is completed(s=8.71m). The beginning is cut in the final plots to simulate a 'warm start'.
"""
def generate_circle():
    n = 10000
    r = 1
    w = 3
    th = np.linspace(0, 2*np.pi, n)
    w = r*np.ones(n)
    p = r*np.array([np.cos(th), np.sin(th)]).T
    return p, w


p,w = generate_circle()
reftrack = np.column_stack(( p, w, w ))
cl = centerline.Centerline(reftrack)
s, reftrack, kappa, normvec = cl.discretize(0, cl.end(), 100)

track = "LMS_Track.txt"
[s_imp, _, _, _, kappa_imp] = getTrack(track)

Sref        = s_imp
s0          = s_imp
kapparef    = kappa_imp

# s0          = s
# kapparef    = kappa

plot = False
if plot:
    plt.subplot(211)
    plt.title("ori")
    plt.plot(s0)
    plt.plot(kapparef)
    plt.subplot(212)
    plt.title("disc")
    plt.plot(s)
    plt.plot(kappa)
    plt.show()

Sref        = s0
Tf = 1.0  # prediction horizon
N = 50  # number of discretization steps
T = 10.00  # maximum simulation time[s]
sref_N = 3  # reference for final reference progress

# load model
constraint, model, acados_solver = acados_settings(Tf, N, s0, kapparef)

# dimensions
nx = model.x.size()[0]
nu = model.u.size()[0]
ny = nx + nu
Nsim = int(T * N / Tf)

# initialize data structs
simX = np.ndarray((Nsim, nx))
simU = np.ndarray((Nsim, nu))
s0 = model.x0[0]
tcomp_sum = 0
tcomp_max = 0

# simulate
for i in range(Nsim):
    # update reference
    sref = s0 + sref_N
    for j in range(N):
        yref = np.array([s0 + (sref - s0) * j / N, 0, 0, 0, 0, 0, 0, 0])
        # yref=np.array([1,0,0,1,0,0,0,0])
        acados_solver.set(j, "yref", yref)
    yref_N = np.array([sref, 0, 0, 0, 0, 0])
    # yref_N=np.array([0,0,0,0,0,0])
    acados_solver.set(N, "yref", yref_N)

    # solve ocp
    t = time.time()

    status = acados_solver.solve()
    if status != 0:
        print("acados returned status {} in closed loop iteration {}/{}.".format(status, i, Nsim))

    elapsed = time.time() - t

    # manage timings
    tcomp_sum += elapsed
    if elapsed > tcomp_max:
        tcomp_max = elapsed

    # get solution
    x0 = acados_solver.get(0, "x")
    u0 = acados_solver.get(0, "u")
    for j in range(nx):
        simX[i, j] = x0[j]
    for j in range(nu):
        simU[i, j] = u0[j]

    # update initial condition
    x0 = acados_solver.get(1, "x")
    acados_solver.set(0, "lbx", x0)
    acados_solver.set(0, "ubx", x0)
    s0 = x0[0]

    # check if one lap is done and break and remove entries beyond
    if x0[0] > Sref[-1] + 0.1:
        # find where vehicle first crosses start line
        N0 = np.where(np.diff(np.sign(simX[:, 0])))[0][0]
        Nsim = i - N0  # correct to final number of simulation steps for plotting
        simX = simX[N0:i, :]
        simU = simU[N0:i, :]
        break

# Plot Results
t = np.linspace(0.0, Nsim * Tf / N, Nsim)
plotRes(simX, simU, t)
plotTrackProj(simX, track)
plotalat(simX, simU, constraint, t)

plt.show()
