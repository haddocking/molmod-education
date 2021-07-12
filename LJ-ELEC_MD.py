#!/usr/local/bin/python 
#
# Simple MD of Lennard Jones charged or uncharges particles
# Alexandre Bonvin, Aalt Jan van Dijk, Utrecht University
#
# adapted from a script from Patrick Fuchs, Uni. Paris VI
#
##################
# import modules #
##################
from math import sqrt,exp,log,sin,cos
from random import random,randint

####################
# define constants #
####################

nAtoms = 5         # number of atoms
Radius = 25.0      # beware that Radius must be in a good range (according to nAtoms)
                   # in order to be able to place all atoms
Rmin = 2 * Radius  # distance at which rmin is mini
BoxDim = [500,500] # box dimension
Atom_Coord = []    # list of the form : [NMAX][2]
Epsilon = 20.0     # well depth
Dielec = 1.0       # dielectric constant
qat = Radius       # Atom absolute charge
frac_neg = 0.5     # Fraction negative charges
OverlapFr = 0.0    # fraction of overlap allowed
CutOff = 250       # non-bonded cutoff
CutOffSquare = CutOff**2
speed = 20         # canvas update speed
cstboltz = 0.00198722         # Boltzmann's constant in kcal/mol/K
cstboltz = 1000*cstboltz/4.18 #in J/mol/K
Temperature = 300.0           # temperature in K
timestep=1.0E-2	              # MD time step

##################
# some functions #
##################

def dist(A,B):
  return sqrt((A[0]-B[0])**2+(A[1]-B[1])**2)

# change sign
def SignR(a,b):
  if b > 0:
    return a
  else:
    return -a

def charge_color(charge,qat):
  tmp = "#111111"
  if charge == qat:
    tmp = "#FFFFFF"
  else:
    tmp = "#333333"
  return tmp

def calc_temp(vel,nat,k):
   mass=1.0				#mass is set to 1.0
   v2=0.0 
   for i in range(len(vel)):
        vx=vel[i][0]			#velocity
        vy=vel[i][1]
   v2=v2+vx**2+vy**2
   nkt=v2*0.5*mass		#kinetic energy equals 0.5*m*v**2
   temp=nkt/(nat*k)			#N*k*T=Kinetic Energy
   return temp
   
#########################
# initialize parameters #
#########################

# generates random coordinates
def InitConf(n,dim,radius,qat,frac_neg):
 print("Initializing box, please wait...")
 # generate a list of random positions
 tmp_coord = []
 i = 0
 nneg = 0
 ntrial = 0
 #for testing purposes:
 # fix first atom
 x = random()*(dim[0]-radius)+radius#dim[0]
 y = random()*(dim[1]-radius)+radius#dim[1]
 nneg = int(float(n) * frac_neg)
 npos = n - nneg
 charge = -qat
 if (npos == n): charge = qat
 i += 1
 if (n==2):
   tmp_coord.append([175,300,charge])
 else:
   tmp_coord.append([x,y,charge])
#  print "atom ",i, charge, frac_neg, float(nneg)/float(i)
 while(i < nneg):
    x = random()*(dim[0]-radius)+radius#dim[0]
    y = random()*(dim[1]-radius)+radius#dim[1]
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
#      print "atom ",i, charge, frac_neg, float(nneg)/float(i)
    ntrial = ntrial + 1
    if ntrial > 100000: 
      print("initialisation failed")
      print("==> reduce radius or number of atoms")
      sys.exit()
 while(i < n):
    x = random()*(dim[0]-radius)+radius#dim[0]
    y = random()*(dim[1]-radius)+radius#dim[1]
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
#      print "atom ",i, charge, frac_neg, float(nneg)/float(i)
    ntrial = ntrial + 1
    if ntrial > 10**10: 
      print("initialisation failed")
      print("==> reduce radius or number of atoms")
      sys.exit()
 return tmp_coord

# generates random charges
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

#generates initial velocities according to Maxwell distribution
def InitVel(n,temperature):
   mass=1.0
   stdev=sqrt(cstboltz*temperature/mass)
   print("initializing velocities, please wait..")
   tmp_vel=[]
   #for testing purposes:
   if (n==2):
    vel=[[0,0],[0,0]]
   else:
    i=0
    while (i<n):
      r1=random() 
      r2=random() 
      #generate random numbers according to Gaussian:
      x1=sqrt(-2.0*log(r1))*cos(r2)
      x2=sqrt(-2.0*log(r1))*sin(0.5*r2)
      x1=x1*stdev
      x2=x2*stdev
      tmp_vel.append([x1,x2])
      i+=1
    #remove overall motion
    i=0
    vxt=0.0
    vyt=0.0
    while (i<n):
     vxt+=tmp_vel[i][0]
     vyt+=tmp_vel[i][1]
     i+=1
    i=0
    while (i<n):
     tmp_vel[i][0]=tmp_vel[i][0]-vxt/float(n)
     tmp_vel[i][1]=tmp_vel[i][1]-vyt/float(n)
     i+=1
    #scaling factor is used to get temperature exactly equal to desired temperature
    tt=calc_temp(tmp_vel,n,cstboltz)
    scaling=sqrt(temperature/tt)
    vel=[]
    i=0
    while (i<n):
     vx,vy=tmp_vel[i][0],tmp_vel[i][1]
     vx=vx*scaling
     vy=vy*scaling
     vel.append([vx,vy])
     i+=1
   return vel

####################
# calculate energy #
####################
# classical LJ
def LJ(r,epsilon,rmin):
  return epsilon*((rmin/r)**12-(rmin/r)**6)

# classical Coulomb
def Coulomb(r,dielec,qa,qb):
  return qa*qb/(dielec*r)

# classical Coulomb2
# argument is squared distance
def Coulomb2(r,dielec,qa,qb):
  return qa*qb/(dielec*sqrt(r))

# version without boundary conditions
def Calc_Ene2(coord,epsilon,rmin,dielec,cutoffsquare,boxdim,elec=1):
  Ene = 0.0 ; distsquare = 0
  ELJ = 0.0; ECoul=0.0
  rmin_exp6 = rmin**6
  # doubly nested loop over all particule pairs
  for i in range(len(coord)-1):
    for j in range(i+1,len(coord)):
      # calculate the squared atomic distance
      for k in range(2):
        tmp = coord[j][k] - coord[i][k]
        distsquare += tmp**2
      qa = coord[i][2]
      qb = coord[j][2]
      LJ  = LJ2(distsquare, epsilon, rmin_exp6)
      Ene += LJ
      ELJ += LJ
      if (elec):
          CC = Coulomb2(distsquare,dielec,qa,qb)
          Ene+=CC
          ECoul+=CC
  return Ene,ELJ,ECoul

# calculate LJ from the squared distance
def LJ2(distsquare, epsilon, rmin_exp6):
  Z = (1/distsquare)**3 * rmin_exp6
  return epsilon * Z * (Z-1)

def ForceLJ2(distsquare, epsilon, rmin_exp6,xi):
  rij=sqrt(distsquare)
  Z = (1/distsquare)**3 * rmin_exp6
  dedz=epsilon*(2*Z-1)
  dzdr=rmin_exp6*(-6.0/rij**(7.0))
  drdx=xi/rij
  return dedz*dzdr*drdx

def ForceCoulomb(distsquare,dielec,qa,qb,xi):
  rij=sqrt(distsquare)
  dedr=-1.0*(qa*qb/dielec)*(1/distsquare)
  drdx=xi/rij
  return dedr*drdx

# version with boundary conditions
def Calc_Ene(coord,epsilon,rmin,dielec,cutoffsquare,boxdim,elec=1):
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
        LJ  = LJ2(distsquare, epsilon, rmin_exp6)
        Ene += LJ
        ELJ += LJ
        if (elec):
          CC = Coulomb2(distsquare,dielec,qa,qb)
          Ene+=CC
          ECoul+=CC
  return Ene,ELJ,ECoul

def Calc_Kin(vel):
   mass=1.0
   v2=0.0 
   for i in range(len(vel)):
        vx=vel[i][0]			#velocity
        vy=vel[i][1]
   v2=v2+vx**2+vy**2
   kin=v2*0.5*mass
   return kin

def Calc_Force(coord,epsilon,rmin,dielec,cutoffsquare,boxdim):
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
        ff += ForceCoulomb(distsquare,dielec,qa,qb,tmp)
        fflist.append(ff)
       for k in range(2):
        tmpforce[k]=tmpforce[k]+fflist[k]
    Force.append(tmpforce)
  return Force


#####################
# move particles   #  
#####################
def Step1(atom_coord,velocity,force,h):
    #the first step for Verlet
    list=[]
    for i in range(len(atom_coord)):
        q=atom_coord[i][2]
        r0x=atom_coord[i][0]			#coordinates
        r0y=atom_coord[i][1]
        v0x=velocity[i][0]			#velocity
        v0y=velocity[i][1]
        r0x=r0x+h*v0x+0.5*h**2*force[i][0]
        r0y=r0y+h*v0y+0.5*h**2*force[i][1]
        list.append([r0x,r0y,q])
    return list

def Verlet(atom_coord,velocity,force,h,old_atom_coord):
    #the Verlet algorithm
    list=[]
    #NB: put mass m to 1
    for i in range(len(atom_coord)):
        q=atom_coord[i][2]
        r0x=atom_coord[i][0]			#coordinates
        r0y=atom_coord[i][1]
        oldr0x=old_atom_coord[i][0]		#old coordinates
        oldr0y=old_atom_coord[i][1]
        v0x=velocity[i][0]			    #velocity
        v0y=velocity[i][1]
        r0x=2*r0x-oldr0x+h**2*force[i][0]	#verlet
        r0y=2*r0y-oldr0y+h**2*force[i][1]
        list.append([r0x,r0y,q])
    return list

def CalcVel(old_atom_coord,atom_coord,h):
    #calculate velocities based on old and new positions
    list=[]
    #NB: put mass m to 1
    for i in range(len(atom_coord)):
        r0x=atom_coord[i][0]
        r0y=atom_coord[i][1]
        oldr0x=old_atom_coord[i][0]
        oldr0y=old_atom_coord[i][1]
        v0x=(r0x-oldr0x)/(2*h)
        v0y=(r0y-oldr0y)/(2*h)
        list.append([v0x,v0y])
    return list


############################
# move particules with MD  #  
############################
def Go(*args):
  global Atom_Coord,Radius,BoxDim,Epsilon,Rmin,CutOffSquare,Iterations,Ene,Accepted,Old_Atom_Coord
  global Velocity,timestep,hwtextene1,hwtextene2
  global Color,sttext0,ptext2,paccept,Dielec,root,canevas,speed
  hwtextene1.destroy()
  hwtextene2.destroy()
  sttext0.destroy()

  Ene,EneLJ,EneCoul = Calc_Ene(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim)
  Kin = Calc_Kin(Velocity)
  Force = Calc_Force(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim)
  if Iterations==0:
       Old_Atom_Coord=Atom_Coord
       Atom_Coord=Step1(Atom_Coord,Velocity,Force,timestep)
  else:
   tmp=Atom_Coord
   Atom_Coord=Verlet(Atom_Coord,Velocity,Force,timestep,Old_Atom_Coord)
   Velocity=CalcVel(Old_Atom_Coord,Atom_Coord,timestep)
   Old_Atom_Coord=tmp
  
  mynewtext="step: %d Time: %8.3f" % (Iterations,float(Iterations)*timestep)
  sttext0=Label(top1,text=mynewtext)
  sttext0.pack(side='left')

  mynewtext="Etot: %6.1f Ekin: %6.1f Epot: %6.1f" % (Ene+Kin,Kin,Ene)
  hwtextene1=Label(top2,text=mynewtext)
  hwtextene1.pack(side='left')

  temperature=calc_temp(Velocity,nAtoms,cstboltz)
  mynewtext="Elj: %6.1f Ecoul: %6.1f Temp: %6.1f" % (EneLJ,EneCoul,temperature)
  hwtextene2=Label(top3,text=mynewtext)
  hwtextene2.pack(side='left')
  
  if temperature > 1000000: 
    print("The system is exploding !!!")
    print("step: %d Time: %8.3f" % (Iterations,float(Iterations)*timestep))
    print("Etot: %6.1f Ekin: %6.1f Epot: %6.1f" % (Ene+Kin,Kin,Ene))
    print("Elj: %6.1f Ecoul: %6.1f Temp: %6.1f" % (EneLJ,EneCoul,temperature))
    print("Emergency stop")
    sys.exit()

#  #apply boudary conditions
#  for pp in range(len(Atom_Coord)):
#   for i in range(2): # i=0 -> case x coordinate ; i=1 -> case y coordinate
#     if Atom_Coord[pp][i] < 0:
#       Atom_Coord[pp][i] -= BoxDim[i]*(-1+int(Atom_Coord[pp][i]/BoxDim[i]))
#     if Atom_Coord[pp][i] > BoxDim[i]:
#       Atom_Coord[pp][i] -= BoxDim[i]*(int(Atom_Coord[pp][i]/BoxDim[i]))

  #apply boudary conditions
  for pp in range(len(Atom_Coord)):
   for i in range(2): # i=0 -> case x coordinate ; i=1 -> case y coordinate
    if Atom_Coord[pp][i] < 0:
     Atom_Coord[pp][i] += BoxDim[i]
     Old_Atom_Coord[pp][i] += BoxDim[i]
    if Atom_Coord[pp][i] > BoxDim[i]:
     Atom_Coord[pp][i] -= BoxDim[i]
     Old_Atom_Coord[pp][i] -= BoxDim[i]

  # draw new canvas coordinates
  for i in range(len(Atom_Coord)):  
    x1 = Atom_Coord[i][0] + Radius
    y1 = Atom_Coord[i][1] + Radius
    x2 = Atom_Coord[i][0] - Radius
    y2 = Atom_Coord[i][1] - Radius
    canevas.coords(ATOM[i],x1,y1,x2,y2)
  Iterations=Iterations+1

  canevas.after(speed,Go)

######################## 
# some other functions #
########################
def die(event=0):
  import sys
  sys.exit()

###########################################################################
###########################################################################

#some functions added for extra buttons etc.
def stop():
 global Atom_Coord,ATOM,Color,canevas
 global r,nAtoms,size,Radius
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 canevas.pack()

def setupall(Atom_Coord,ATOM,Color,Repack=1):
 global Iterations,Velocity,Temperature
 ATOM = [] # liste contenant les widgets
 if (Repack==1): 
     Atom_Coord = InitConf(nAtoms,BoxDim,Radius,qat,frac_neg)
     Color = []
     for i in range(nAtoms):
         Color.append(charge_color(Atom_Coord[i][2],qat))
     Velocity=InitVel(nAtoms,Temperature) 
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
 return Atom_Coord,ATOM,Color

def set_r(event):
 global Atom_Coord,ATOM,Color,canevas
 global r,nAtoms,Velocity,Iterations
 nAtoms=int(r.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color)
 canevas.pack()
 Iterations=0
 
def set_size(event):
 global Atom_Coord,ATOM,Color,canevas,Radius,Rmin
 global r,nAtoms,size
 Radius=int(size.get())
 Rmin = 2 * Radius
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color)
 canevas.pack()
 
def set_vdw1(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,nAtoms,vdw1,Epsilon
 Epsilon=int(vdw1.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 canevas.pack()
 
def set_vdw2(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,nAtoms,vdw2,Rmin
 Rmin=int(vdw2.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 canevas.pack()

def set_frac(event):
 global Atom_Coord,ATOM,Color,canevas,frac_neg
 global r,nAtoms,temp,Temperature
 frac_neg=float(frac.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,2)
 canevas.pack()

def set_q(event):
 global Atom_Coord,ATOM,Color,canevas,qat
 global r,nAtoms,temp,Temperature
 qat=float(q.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,2)
 canevas.pack()

def set_diel(event):
 global Atom_Coord,ATOM,Color,canevas,Dielec
 global r,nAtoms,temp,Temperature
 Dielec=float(diel.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 canevas.pack()

def set_temp(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,nAtoms,temp,Temperature,Velocity,Iterations
 Temperature=float(temp.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 Velocity=InitVel(nAtoms,Temperature)
 canevas.pack()
 Iterations=0

def set_spd(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,nAtoms,spd,speed
 speed=max(1,int(spd.get()))
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 canevas.pack()

def set_tstep(event):
 global Atom_Coord,ATOM,Color,canevas,Radius,Iterations,timestep,Velocity
 global r,nAtoms,tstep,Temperature
 timestep=float(tstep.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 Velocity=InitVel(nAtoms,Temperature)
 canevas.pack()
 Iterations=0

###########################################################################
###########################################################################

################
# MAIN PROGRAM #
################

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
low2=Frame(root)
low2.pack(side='bottom')
low3=Frame(root)
low3.pack(side='bottom')
mynewtext="Charged particles Molecular Dynamics"
hwtext=Label(top,text=mynewtext,foreground='red',font='times 18 bold')
hwtext.pack(side='left')
mynewtext="nAtoms = " 
hwtext=Label(low,text=mynewtext)
hwtext.pack(side='left')

r = DoubleVar()
size = DoubleVar()
vdw1=DoubleVar()
vdw2=DoubleVar()
frac=DoubleVar()
diel=DoubleVar()
q=DoubleVar()
temp=DoubleVar()
spd=DoubleVar()
tstep=DoubleVar()
swp=DoubleVar()

Color = []
ATOM = [] # liste contenant les widgets
Atom_Coord = InitConf(nAtoms,BoxDim,Radius,qat,frac_neg)
Velocity=InitVel(nAtoms,Temperature)

r.set(nAtoms)
r_entry=Entry(low,width=6,textvariable=r)
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

mynewtext3="T[K] = " 
hwtext3=Label(low,text=mynewtext3)
hwtext3.pack(side='left')
temp.set(Temperature)
temp_entry=Entry(low,width=6,textvariable=temp)
temp_entry.pack(side='left')
temp_entry.bind('<Return>', set_temp)

mynewtext3="Timestep = " 
hwtext3=Label(low,text=mynewtext3)
hwtext3.pack(side='left')
tstep.set(timestep)
tstep_entry=Entry(low,width=6,textvariable=tstep)
tstep_entry.pack(side='left')
tstep_entry.bind('<Return>', set_tstep)

startbutton=Button(top3,text=' start ',command=Go,background='yellow',foreground='blue',relief='flat')
startbutton.pack(side='left',fill='x')
stopbutton=Button(top3,text=' stop ',command=stop,background='yellow', foreground='red',relief='flat')
stopbutton.pack(side='left')

mynewtext3="Delay = " 
hwtext3=Label(low,text=mynewtext3)
hwtext3.pack(side='left')
spd.set(speed)
spd_entry=Entry(low,width=6,textvariable=spd)
spd_entry.pack(side='left')
spd_entry.bind('<Return>', set_spd)

# set up a global var with the nb of iterations
Iterations = 0 

# calculate first energy
Ene,EneLJ,EneCoul = Calc_Ene(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim)
Kin = Calc_Kin(Velocity)

mynewtext="Etot: %6.1f Ekin: %6.1f Epot: %6.1f" % (Ene+Kin,Kin,Ene)
hwtextene1=Label(top1,text=mynewtext)
hwtextene1.pack(side='left')

temperature=calc_temp(Velocity,nAtoms,cstboltz)
mynewtext="Elj: %6.1f Ecoul: %6.1f Temp: %6.1f" % (EneLJ,EneCoul,temperature)
hwtextene2=Label(top2,text=mynewtext)
hwtextene2.pack(side='left')

mynewtext="step: %d Time: %8.3f" % (Iterations,float(Iterations)*timestep)
sttext0=Label(top3,text=mynewtext)
sttext0.pack(side='left')

Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color)

canevas.pack()

print("Click on mouse button within the box to go ahead !!!")
print("<ESC> to quit")

# conserver la fenetre ouverte (inutile en interactif)
root.mainloop()
