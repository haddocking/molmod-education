#!/usr/local/bin/python 
#
###########################################################################
###########################################################################
# Simple EM of Lennard Jones charged or uncharges particles
# Alexandre Bonvin, Utrecht University
#
# adapted from a script from Patrick Fuchs, Uni. Paris VI
#
###########################################################################
###########################################################################

##################
# import modules #
##################
from math import sqrt,exp,log,sin,cos
from random import random,randint,seed

###########################################################################
###########################################################################
### define parameters #####################################################
###########################################################################
###########################################################################

nAtoms = 20        # number of atoms
Radius = 25.0      # beware that Radius must be in a good range (according to nAtoms)
                   # in order to be able to place all atoms
Rmin = 2.24 * Radius  # distance at which rmin is mini
BoxDim = [500,500] # box dimension
Atom_Coord = []    # list of the form : [NMAX][2]
Epsilon = 25.0     # well depth
Dielec = 1.0       # dielectric constant
qat = Radius       # Atom absolute charge
frac_neg = 0.5     # Fraction negative charges
OverlapFr = 0.0    # fraction of overlap allowed
CutOff = 250       # non-bonded cutoff
CutOffSquare = CutOff**2
speed = 20         # canvas update speed
cstboltz = 0.00198722         # Boltzmann's constant in kcal/mol/K
cstboltz = 1000*cstboltz/4.18 #in J/mol/K
drinit = 1.00	   # dr from EM
drmin  = 0.00001   # minimum dr value to step EM
drmax  = 5.00      # maximum dr
alpha  = 1.05      # scaling factor for dr if Enew < Eold  
beta   = 0.90      # scaling factor for dr if Enew > Eold
deltaE = 0.001     # energy difference threshold to stop EM
normFmin = 0.001   # minimum force norm to step EM
Seed   = 100       # random number seed




###########################################################################
###########################################################################
# Steepest descent minimizer ##############################################
###########################################################################
###########################################################################

def Steepest_descent(atom_coord,drstep,force):
    #the first step for Verlet
    list=[]

# This function gets as input parameters:
# - atom_coord, a vector containing the x and y position and the charge of the i atoms
# - drstep, the displacement for the minimizer
# - force, a vector containing the x and y components of the force on the atoms
#
# The function return a list array (vector containing the new positions)
# 
# Implement in the following loop over all atoms the steepest descent algorithm
#
# A few hints:
# - powers in python are given by **, e.g.: x to the square is x**2
# - squared root x: sqrt(x)
# - avoid dividing by zero
#
# 1) First calculate the norm of the total force vector
#
    normf = 0.0
    for i in range(len(atom_coord)):
        normf=normf+force[i][0]**2.0+force[i][1]**2.0
    normf=sqrt(normf)
# 
# 2) Then move the particles
#
    for i in range(len(atom_coord)):
        q=atom_coord[i][2]
        r0x=atom_coord[i][0]		#coordinates
        r0y=atom_coord[i][1]
        if (normf > 0):
#
# Insert below the lines defining the new coordinates based on the old ones + forces + drstep
#
# forces are contained in force[i][0] for the x force component and force[i][1] for the y force component
# the step size for the move is given by drstep
#
# 
# ====>>>>>

           r0xnew=r0x
           r0ynew=r0y
# <<<<<====
        r0x=r0xnew
        r0y=r0ynew
        list.append([r0x,r0y,q])
    return list,normf



###########################################################################
###########################################################################
# move particules with EM  ################################################
###########################################################################
###########################################################################

def Go(*args):
  import sys,time
  global Atom_Coord,Radius,BoxDim,Epsilon,Rmin,Cutoff,CutOffSquare,Iterations,Ene,Ene_prev,Accepted
  global drstep,drmax,drmin,deltaE,hwtextene1,hwtextene2,alpha,beta,normFmin
  global Color,sttext0,ptext2,paccept,Dielec,root,canevas,speed,nout
  hwtextene1.destroy()
  hwtextene2.destroy()
  sttext0.destroy()

  if Iterations==0:
    nout = 0
    drstep = drinit
    Ene,EneLJ,EneCoul = Calc_Ene2(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim)
    Ene_prev=Ene
    outtext="Iteration: %8d Epot: %6.1f Elj: %6.1f Ecoul: %6.1f" % (Iterations,Ene,EneLJ,EneCoul)
    print(outtext)
  nout = nout + 1
  Force = Calc_Force2(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim)
  Atom_Coord, normF=Steepest_descent(Atom_Coord,drstep,Force)
  Ene,EneLJ,EneCoul = Calc_Ene2(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim)
  Ene_diff= Ene - Ene_prev  
  if (Ene_diff < 0.0):
    drstep = drstep * alpha
    drstep = min(drmax,drstep)
  else:
    drstep = drstep * beta
  Ene_prev=Ene

  #update energies
  mynewtext="step: %d" % (Iterations)
  sttext0=Label(top1,text=mynewtext)
  sttext0.pack(side='left')
  mynewtext="Epot: %6.1f deltaE: %10.6f dr: %8.6f" % (Ene,Ene_diff,drstep)
  hwtextene1=Label(top2,text=mynewtext)
  hwtextene1.pack(side='left')
  mynewtext="Elj: %6.1f Ecoul: %6.1f" % (EneLJ,EneCoul)
  hwtextene2=Label(top3,text=mynewtext)
  hwtextene2.pack(side='left')
  
  #apply boudary conditions
  for pp in range(len(Atom_Coord)):
   for i in range(2): # i=0 -> case x coordinate ; i=1 -> case y coordinate
    if Atom_Coord[pp][i] < 0:
     Atom_Coord[pp][i] += BoxDim[i]
    if Atom_Coord[pp][i] > BoxDim[i]:
     Atom_Coord[pp][i] -= BoxDim[i]

  #draw new canvas coordinates
  for i in range(len(Atom_Coord)):  
    x1 = Atom_Coord[i][0] + Radius
    y1 = Atom_Coord[i][1] + Radius
    x2 = Atom_Coord[i][0] - Radius
    y2 = Atom_Coord[i][1] - Radius
    canevas.coords(ATOM[i],x1,y1,x2,y2)
  Iterations=Iterations+1

  #print to terminal window
  normF = normF/len(Atom_Coord)
  if nout == 20:
    nout = 0
    outtext="Iteration: %8d Epot: %6.1f Elj: %6.1f Ecoul: %6.1f deltaE: %10.6f <normF>: %8.6f dr: %8.6f" % (Iterations,Ene,EneLJ,EneCoul,Ene_diff,normF,drstep)
    print(outtext)

  if (abs(Ene_diff) < deltaE or drstep < drmin or normF < normFmin ):
    print("STOPPING... deltaE<",deltaE,", or drstep<",drmin,", or normF<",normFmin)
    outtext="Iteration: %8d Epot: %6.1f Elj: %6.1f Ecoul: %6.1f deltaE: %10.6f <normF>: %8.6f dr: %8.6f" % (Iterations,Ene,EneLJ,EneCoul,Ene_diff,normF,drstep)
    print(outtext)
  else:
    canevas.after(speed,Go)    

def reset(*args):
  import sys,time
  global Atom_Coord,ATOM,Radius,BoxDim,Epsilon,Rmin,Cutoff,CutOffSquare,Iterations,Ene,Ene_prev,Accepted
  global drstep,drmax,drmin,deltaE,hwtextene1,hwtextene2,alpha,beta,normFmin
  global Color,sttext0,ptext2,paccept,Dielec,root,canevas,speed,nout
  hwtextene1.destroy()
  hwtextene2.destroy()
  sttext0.destroy()
  canevas.destroy()
  canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
  canevas.bind("<Button-1>",Go)
  Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color)
  update_ene()
  canevas.pack()
  Iterations=0


###########################################################################
###########################################################################
### energy functions  #####################################################
###########################################################################
###########################################################################

# calculate LJ from the squared distance
def LJ2(distsquare, epsilon, rmin_exp6):
  Z = (1/distsquare)**3 * rmin_exp6
  return epsilon * Z * (Z-1)

# classical Coulomb from the squared distance
def Coulomb2(r,dielec,qa,qb):
  return qa*qb/(dielec*sqrt(r))

# Calculate energy Evdw + Ecoulomb (used squared distance)
# version with boundary conditions
def Calc_Ene2(coord,epsilon,rmin,dielec,cutoffsquare,boxdim,elec=1):
  Ene = 0.0 ; distsquare = 0
  ELJ = 0.0; ECoul=0.0
  rmin_exp6 = rmin**6
  # doubly nested loop over all particule pairs
  for i in range(len(coord)-1):
    for j in range(i+1,len(coord)):
      # calculate the squared atomic distance
      distsquare = 0
      for k in range(2):
        tmp = coord[j][k] - coord[i][k]
        # chooses the nearest image
        halfbox = boxdim[k]/2 
        tmp = tmp - SignR(halfbox,tmp-halfbox) - SignR(halfbox,tmp+halfbox)
        distsquare += tmp**2
      # compute vdw and Coulomb energy
      if distsquare < cutoffsquare:
        qa = coord[i][2]
        qb = coord[j][2]
        vdw  = LJ2(distsquare, epsilon, rmin_exp6)
        Ene += vdw
        ELJ += vdw
        if (elec):
          CC = Coulomb2(distsquare,dielec,qa,qb)
          Ene+=CC
          ECoul+=CC
  return Ene,ELJ,ECoul

###########################################################################
###########################################################################
### force functions  ######################################################
###########################################################################
###########################################################################

# force LJ (use squared distance)
def ForceLJ2(distsquare, epsilon, rmin_exp6,xi):
  rij=sqrt(distsquare)
  Z = (1/distsquare)**3 * rmin_exp6
  dedz=epsilon*(2*Z-1)
  dzdr=rmin_exp6*(-6.0/rij**(7.0))
  drdx=xi/rij
  return dedz*dzdr*drdx

# Force Coulomb (use squared distance)
def ForceCoulomb2(distsquare,dielec,qa,qb,xi):
  rij=sqrt(distsquare)
  dedr=-1.0*(qa*qb/dielec)*(1/distsquare)
  drdx=xi/rij
  return dedr*drdx

# Calculate force from Evdw + Ecoulomb (uses squared distance)
def Calc_Force2(coord,epsilon,rmin,dielec,cutoffsquare,boxdim):
  Force=[] ; distsquare = 0
  rmin_exp6 = rmin**6
  # doubly nested loop over all particle pairs
  for i in range(len(coord)):
    tmpforce=[0.0,0.0]
    for j in range(len(coord)):
     if not (i==j):
      # calculate the squared atomic distance
      distsquare = 0
      for k in range(2):
        tmp = coord[j][k] - coord[i][k]
        # chooses the nearest image
        halfbox = boxdim[k]/2 
        tmp = tmp - SignR(halfbox,tmp-halfbox) - SignR(halfbox,tmp+halfbox)
        distsquare += tmp**2
      # compute vdw force
      if distsquare < cutoffsquare:
       qa = coord[i][2]
       qb = coord[j][2]
       fflist=[]
       for k in range(2):
        tmp = coord[j][k] - coord[i][k]
        ff = ForceLJ2(distsquare, epsilon, rmin_exp6,tmp)
        ff += ForceCoulomb2(distsquare,dielec,qa,qb,tmp)
        fflist.append(ff)
       for k in range(2):
        tmpforce[k]=tmpforce[k]+fflist[k]
    Force.append(tmpforce)
  return Force


###########################################################################
###########################################################################
### other functions  ######################################################
###########################################################################
###########################################################################


### distance ###
def dist(A,B):
  return sqrt((A[0]-B[0])**2+(A[1]-B[1])**2)

### squared distance ###
def dist2(A,B):
  return (A[0]-B[0])**2+(A[1]-B[1])**2

### change sign ###
def SignR(a,b):
  if b > 0:
    return a
  else:
    return -a

### color particules based on charge ###
def charge_color(charge,qat):
  tmp = "#111111"
  if charge == qat:
    tmp = "#FFFFFF"
  else:
    tmp = "#333333"
  return tmp

def die(event=0):
  import sys
  sys.exit()
   
###########################################################################
###########################################################################
### initialization  #######################################################
###########################################################################
###########################################################################

### generates random coordinates ###
def InitConf(n,dim,radius,qat,frac_neg):
 seed(Seed)
 print("Initializing box, please wait...")
 # generate a list of random positions
 tmp_coord = []
 i = 0
 nneg = 0
 ntrial = 0
 # fix first atom
 x = random()*(dim[0]-2*radius)+radius#dim[0]
 y = random()*(dim[1]-2*radius)+radius#dim[1]
 nneg = int(float(n) * frac_neg)
 npos = n - nneg
 charge = -qat
 if (npos == n): charge = qat
 i += 1
 if (n==2):
   tmp_coord.append([175,300,charge])
 else:
   tmp_coord.append([x,y,charge])
 while(i < nneg):
    x = random()*(dim[0]-2*radius)+radius#dim[0]
    y = random()*(dim[1]-2*radius)+radius#dim[1]
    # check wether the new particule ovelap an existing one
    OVERLAP = 1
    for j in range(i):
      if dist(tmp_coord[j],[x,y]) < (1-OverlapFr)*2*radius:
        OVERLAP = 0
    if OVERLAP:
      charge = -qat
      if (n==2):
        tmp_coord.append([325,300,charge])
      else:
        tmp_coord.append([x,y,charge])
      i += 1
    ntrial = ntrial + 1
    if ntrial > 100000: 
      print("initialisation failed")
      print("==> reduce radius or number of atoms")
      sys.exit()
 while(i < n):
    x = random()*(dim[0]-2*radius)+radius#dim[0]
    y = random()*(dim[1]-2*radius)+radius#dim[1]
    # check wether the new particule overlap an existing one
    OVERLAP = 1
    for j in range(i):
      if dist(tmp_coord[j],[x,y]) < (1-OverlapFr)*2*radius:
        OVERLAP = 0
    if OVERLAP:
      charge = qat
      if (n==2):
        tmp_coord.append([325,300,charge])
      else:
        tmp_coord.append([x,y,charge])
      i += 1
    ntrial = ntrial + 1
    if ntrial > 10**10: 
      print("initialisation failed")
      print("==> reduce radius or number of atoms")
      sys.exit()
 return tmp_coord

### generates random charges ###
def InitCharge(n,dim,qat,frac_neg):
  global Atom_Coord
  print("Initializing charges, please wait...")
  i = 0
  nneg = 0
  nneg = int(float(n) * frac_neg)
  npos = n - nneg
  charge = -qat
  if (npos == n): charge = qat
  Atom_Coord[i][2]=charge
  i += 1
  while(i < nneg):
      charge = -qat
      Atom_Coord[i][2]=charge
      i += 1
  while(i < n):
      charge = qat
      Atom_Coord[i][2]=charge
      i += 1


###########################################################################
###########################################################################
### various functions for input + layout ##################################
###########################################################################
###########################################################################

def stop():
 global Atom_Coord,ATOM,Color,canevas
 global r,nAtoms,size,Radius,drinit
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 update_ene()
 canevas.pack()
 update_ene()

### setup system ###
def setupall(Atom_Coord,ATOM,Color,Repack=1):
 global Iterations,drinit,drstep
 ATOM = [] # liste contenant les widgets
 if (Repack==1): 
     Atom_Coord = InitConf(nAtoms,BoxDim,Radius,qat,frac_neg)
     Iterations = 0
     drstep = drinit
     Color = []
     for i in range(nAtoms):
         Color.append(charge_color(Atom_Coord[i][2],qat))
 if (Repack==2): 
     InitCharge(nAtoms,BoxDim,qat,frac_neg)
     Color = []
     for i in range(nAtoms):
         Color.append(charge_color(Atom_Coord[i][2],qat))
 for i in range(len(Atom_Coord)):
  x1 = Atom_Coord[i][0] + Radius
  y1 = Atom_Coord[i][1] + Radius
  x2 = Atom_Coord[i][0] - Radius
  y2 = Atom_Coord[i][1] - Radius
  ATOM.append(Canvas.create_oval(canevas,x1,y1,x2,y2,fill=Color[i]))
  update_ene()
 return Atom_Coord,ATOM,Color


### set number of particules ###
def set_r(event):
 global Atom_Coord,ATOM,Color,canevas
 global r,nAtoms,Iterations
 nAtoms=int(r.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color)
 update_ene()
 canevas.pack()
 Iterations=0
 

### set atom Radius ###
def set_size(event):
 global Atom_Coord,ATOM,Color,canevas,Radius,Rmin
 global r,nAtoms,size
 Radius=int(size.get())
 Rmin = 2 * Radius
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color)
 update_ene()
 canevas.pack()
 
### set epsilon for Lennard-Jones ###
def set_vdw1(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,nAtoms,vdw1,Epsilon
 Epsilon=int(vdw1.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 update_ene()
 canevas.pack()

### set sigma for Lennard-Jones ###
def set_vdw2(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,nAtoms,vdw2,Rmin
 Rmin=int(vdw2.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 update_ene()
 canevas.pack()

### set charge fraction ###
def set_frac(event):
 global Atom_Coord,ATOM,Color,canevas,frac_neg
 global r,nAtoms
 frac_neg=float(frac.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,2)
 update_ene()
 canevas.pack()

### set particule charge ###
def set_q(event):
 global Atom_Coord,ATOM,Color,canevas,qat
 global r,nAtoms
 qat=float(q.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,2)
 update_ene()
 canevas.pack()

### set dielectric constant ###
def set_diel(event):
 global Atom_Coord,ATOM,Color,canevas,Dielec
 global r,nAtoms
 Dielec=float(diel.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 update_ene()
 canevas.pack()

### set minimum Force norm difference for stop ###
def set_dFmin(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,nAtoms,Iterations,normFmin
 normFmin=float(Fmin.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 update_ene()
 canevas.pack()

### set minimum Energy difference for stop ###
def set_deltaE(event):
 global Atom_Coord,ATOM,Color,canevas
 global r,nAtoms,Iterations,deltaE
 deltaE=float(Emin.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 update_ene()
 canevas.pack()

### set refresh delay for graphical update ###
def set_spd(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,nAtoms,spd,speed
 speed=max(1,int(spd.get()))
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 update_ene()
 canevas.pack()

### set initial displacement for minimizer ###
def set_dxstep(event):
 global Atom_Coord,ATOM,Color,canevas,Radius,Iterations
 global r,nAtoms,drstep,drinit
 drinit=float(dxstep.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 update_ene()
 canevas.pack()

### set alpha factor for increasing dr ###
def set_alpha(event):
 global Atom_Coord,ATOM,Color,canevas,Radius,Iterations
 global r,nAtoms,drstep,drinit,alpha
 alpha=float(alphafactor.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 update_ene()
 canevas.pack()

### set beta factor for decreasing dr ###
def set_beta(event):
 global Atom_Coord,ATOM,Color,canevas,Radius,Iterations
 global r,nAtoms,drstep,drinit,beta
 beta=float(betafactor.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 update_ene()
 canevas.pack()

### update energy ###
def update_ene():
  global Atom_Coord,Radius,BoxDim,Epsilon,Rmin,CutOff,CutOffSquare,Iterations,Ene,Ene_prev,Accepted
  global drstep,drmax,deltaE,hwtextene1,hwtextene2
  global Color,sttext0,ptext2,paccept,Dielec,root,canevas,speed
#  Ene,EneLJ,EneCoul = Calc_Ene(Atom_Coord,Epsilon,Rmin,Dielec,CutOff,BoxDim)
  Ene,EneLJ,EneCoul = Calc_Ene2(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim)
  hwtextene1.destroy()
  hwtextene2.destroy()

  mynewtext="Epot: %6.1f deltaE: %10.6f dr: %8.6f" % (Ene,Ene_diff,drstep)
  hwtextene1=Label(top1,text=mynewtext)
  hwtextene1.pack(side='left')

  mynewtext="Elj: %6.1f Ecoul: %6.1f" % (EneLJ,EneCoul)
  hwtextene2=Label(top2,text=mynewtext)
  hwtextene2.pack(side='left')
 

###########################################################################
###########################################################################
### MAIN PROGRAM ##########################################################
###########################################################################
###########################################################################

from tkinter import *
from tkinter import Canvas

root = Tk() #root (main) window

canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
canevas.bind("<Button-1>",Go)
root.bind("<Escape>",die)
top=Frame(root)
top.pack(side='top')
top1=Frame(root)
top1.pack(side='top')
top2=Frame(root)
top2.pack(side='top')
top3=Frame(root)
top3.pack(side='top')
low=Frame(root)
low.pack(side='bottom')
low1=Frame(root)
low1.pack(side='bottom')
low2=Frame(root)
low2.pack(side='bottom')
low3=Frame(root)
low3.pack(side='bottom')
mynewtext="Charged particles Energy Minimization"
hwtext=Label(top,text=mynewtext,foreground='red',font='times 18 bold')
hwtext.pack(side='left')
mynewtext="nAtoms = " 
hwtext=Label(low1,text=mynewtext)
hwtext.pack(side='left')

r = DoubleVar()
size = DoubleVar()
vdw1=DoubleVar()
vdw2=DoubleVar()
frac=DoubleVar()
diel=DoubleVar()
Emin=DoubleVar()
Fmin=DoubleVar()
alphafactor=DoubleVar()
betafactor=DoubleVar()
q=DoubleVar()
temp=DoubleVar()
spd=DoubleVar()
dxstep=DoubleVar()
swp=DoubleVar()
drstep=drinit

Color = []
ATOM = [] # liste contenant les widgets
Atom_Coord = InitConf(nAtoms,BoxDim,Radius,qat,frac_neg)
Iterations = 0

r.set(nAtoms)
r_entry=Entry(low1,width=6,textvariable=r)
r_entry.pack(side='left')
r_entry.bind('<Return>', set_r)

mynewtext2="VDW param: Radius = " 
hwtext2=Label(low3,text=mynewtext2)
hwtext2.pack(side='left')
size.set(Radius)
size_entry=Entry(low3,width=6,textvariable=size)
size_entry.pack(side='left')
size_entry.bind('<Return>', set_size)

mynewtext3="Epsilon = " 
hwtext3=Label(low3,text=mynewtext3)
hwtext3.pack(side='left')
vdw1.set(Epsilon)
vdw1_entry=Entry(low3,width=6,textvariable=vdw1)
vdw1_entry.pack(side='left')
vdw1_entry.bind('<Return>', set_vdw1)

mynewtext2="Coulomb param: frac neg charge = " 
hwtext2=Label(low2,text=mynewtext2)
hwtext2.pack(side='left')
frac.set(frac_neg)
frac_entry=Entry(low2,width=6,textvariable=frac)
frac_entry.pack(side='left')
frac_entry.bind('<Return>', set_frac)

mynewtext3="abs charge = " 
hwtext3=Label(low2,text=mynewtext3)
hwtext3.pack(side='left')
q.set(qat)
q_entry=Entry(low2,width=6,textvariable=q)
q_entry.pack(side='left')
q_entry.bind('<Return>', set_q)

mynewtext3="Dielec = " 
hwtext3=Label(low2,text=mynewtext3)
hwtext3.pack(side='left')
diel.set(Dielec)
diel_entry=Entry(low2,width=6,textvariable=diel)
diel_entry.pack(side='left')
diel_entry.bind('<Return>', set_diel)

mynewtext3="DeltaE threshold = " 
hwtext3=Label(low1,text=mynewtext3)
hwtext3.pack(side='left')
Emin.set(deltaE)
Emin_entry=Entry(low1,width=6,textvariable=Emin)
Emin_entry.pack(side='left')
Emin_entry.bind('<Return>', set_deltaE)

mynewtext3="dr = " 
hwtextdr=Label(low,text=mynewtext3)
hwtextdr.pack(side='left')
dxstep.set(drinit)
dxstep_entry=Entry(low,width=6,textvariable=dxstep)
dxstep_entry.pack(side='left')
dxstep_entry.bind('<Return>', set_dxstep)

mynewtext3="alpha = " 
hwtext3=Label(low,text=mynewtext3)
hwtext3.pack(side='left')
alphafactor.set(alpha)
alphafactor_entry=Entry(low,width=6,textvariable=alphafactor)
alphafactor_entry.pack(side='left')
alphafactor_entry.bind('<Return>', set_alpha)

mynewtext3="beta = " 
hwtext3=Label(low,text=mynewtext3)
hwtext3.pack(side='left')
betafactor.set(beta)
betafactor_entry=Entry(low,width=6,textvariable=betafactor)
betafactor_entry.pack(side='left')
betafactor_entry.bind('<Return>', set_beta)

resetbutton=Button(top3,text=' reset ',command=reset,background='yellow',foreground='green',relief='flat')
resetbutton.pack(side='left',fill='x')
startbutton=Button(top3,text=' start ',command=Go,background='yellow',foreground='blue',relief='flat')
startbutton.pack(side='left',fill='x')
stopbutton=Button(top3,text=' stop ',command=stop,background='yellow', foreground='red',relief='flat')
stopbutton.pack(side='left')

mynewtext3="Delay = " 
hwtext3=Label(low1,text=mynewtext3)
hwtext3.pack(side='left')
spd.set(speed)
spd_entry=Entry(low1,width=6,textvariable=spd)
spd_entry.pack(side='left')
spd_entry.bind('<Return>', set_spd)

# set up a global var with the nb of iterations
Iterations = 0 

# calculate first energy
Ene,EneLJ,EneCoul = Calc_Ene2(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim)
Ene_prev=Ene
Ene_diff=abs(Ene_prev-Ene)
mynewtext="Epot: %6.1f deltaE: %10.6f dr: %8.6f" % (Ene,Ene_diff,drstep)
hwtextene1=Label(top1,text=mynewtext)
hwtextene1.pack(side='left')

mynewtext="Elj: %6.1f Ecoul: %6.1f" % (EneLJ,EneCoul)
hwtextene2=Label(top2,text=mynewtext)
hwtextene2.pack(side='left')

mynewtext="step: %d" % (Iterations)
sttext0=Label(top3,text=mynewtext)
sttext0.pack(side='left')


Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color)
update_ene()
drstep = drinit

canevas.pack()

print("Click on mouse button within the box to go ahead !!!")
print("<ESC> to quit")

# conserver la fenetre ouverte (inutile en interactif)
root.mainloop()
