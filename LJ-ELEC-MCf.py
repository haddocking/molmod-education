#!/usr/local/bin/python 

##################
# import modules #
##################
from math import sqrt,exp
from random import random,randint,seed

####################
# define constants #
####################

nAtoms = 40         # nb of atoms
Radius = 25.0      # beware that Radius must be in a good range (according to nAtoms)
                   # in order to be able to place all atoms
Rmin = 2 * Radius  # distance at which rmin is mini
BoxDim = [500,500] # box dimension
Atom_Coord = []    # list of the form : [NMAX][2]
deltaRmax = 50.0   # step up to which particules move
Epsilon = 20.0     # well depth
Dielec = 1.0       # dielectric constant
qat = Radius       # Atom absolute charge
frac_neg = 0.5     # Fraction negative charges
frac_swap = 0.2    # Fraction of MC steps with charge swapping
OverlapFr = 0.0    # fraction of overlap allowed
CutOff = 250       # non-bonded cutoff
CutOffSquare = CutOff**2
speed = 50         # canvas update speed
cstboltz = 8.3502E-03  # Boltzmann's constant in in J/mol/K
#cstboltz = 0.00198722 # Boltzmann's contstant in cal/mol/K
Temperature = 300.0 # temperature in K
Seed = 100         # random seed

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

# generate a random rgb color like  #xxxxxx (xx should be an hexadecimal nb)
def random_color():
  tmp = "#"
  for i in range(6):
    rdm = randint(0,15)
    tmp += hex(rdm)[-1]
  return tmp
# generate a rgb color based on charge like  #xxxxxx (xx should be an hexadecimal nb)
def charge_color(charge,qat):
  tmp = "#111111"
  if charge == qat:
    tmp = "#FFFFFF"
  else:
    tmp = "#333333"
  return tmp
  

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
  # fix first atom
  x = random()*(dim[0]-radius)+radius#dim[0]
  y = random()*(dim[1]-radius)+radius#dim[1]
  nneg = int(float(n) * frac_neg)
  npos = n - nneg
  charge = -qat
  if (npos == n): charge = qat
  i += 1
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
    # check wether the new particule ovelap an existing one
    OVERLAP = 1
    for j in range(i):
      if dist(tmp_coord[j],[x,y]) < (1-OverlapFr)*2*radius:
        OVERLAP = 0
    if OVERLAP:
      charge = qat
      tmp_coord.append([x,y,charge])
      i += 1
#      print "atom ",i, charge, frac_neg, float(nneg)/float(i)
    ntrial = ntrial + 1
    if ntrial > 100000: 
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
def Coulomb2(r,dielec,qa,qb):
  return qa*qb/(dielec*sqrt(r))

# version without boundary conditions
def Calc_Ene2(coord,epsilon,rmin,dielec,cutoffsquare,boxdim,elec=1):
  Ene = 0.0; distsquare = 0
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
      Ene += LJ2(distsquare, epsilon, rmin_exp6)
      if (elec): Ene += Coulomb2(distsquare,dielec,qa,qb)
  return Ene

# version without boundary conditions, single particles
def Calc_Ene2Single(nat,coord,epsilon,rmin,dielec,cutoffsquare,boxdim,elec=1):
  Ene = 0.0; distsquare = 0
  rmin_exp6 = rmin**6
  # doubly nested loop over all particule pairs
  i = nat
  for j in range(0,len(coord)-1):
    if (j != i):
      # calculate the squared atomic distance
      for k in range(2):
        tmp = coord[j][k] - coord[i][k]
        distsquare += tmp**2
      qa = coord[i][2]
      qb = coord[j][2]
      Ene += LJ2(distsquare, epsilon, rmin_exp6)
      if (elec): Ene += Coulomb2(distsquare,dielec,qa,qb)
  return Ene

# calculate LJ from the squared distance
def LJ2(distsquare, epsilon, rmin_exp6):
  Z = (1/distsquare)**3 * rmin_exp6
  return epsilon * Z * (Z-1)

# version with boundary conditions
def Calc_Ene(coord,epsilon,rmin,dielec,cutoffsquare,boxdim,elec=1):
  Ene = 0.0 ; distsquare = 0
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
        Ene += LJ2(distsquare, epsilon, rmin_exp6)
        if (elec): Ene += Coulomb2(distsquare,dielec,qa,qb)
  return Ene

# version with boundary conditions, single particle
def Calc_EneSingle(nat,coord,epsilon,rmin,dielec,cutoffsquare,boxdim,elec=1):
  Ene = 0.0 ; distsquare = 0
  rmin_exp6 = rmin**6
  # doubly nested loop over all particule pairs
  i = nat
  for j in range(0,len(coord)-1):
    if j != i:
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
        Ene += LJ2(distsquare, epsilon, rmin_exp6)
        if (elec): Ene += Coulomb2(distsquare,dielec,qa,qb)
  return Ene

############################
# move particules in a MC  #  
############################
def Go(*args):
  global Atom_Coord,Radius,deltaRmax,BoxDim,Epsilon,Rmin,CutOffSquare,Iterations,Ene,Accepted
  global Color,sttext0,entext1,ptext2,tmove,paccept,Dielec,root,canevas,proba_boltzmann,xx
  if Iterations > 0: 
    sttext0.destroy()
    entext1.destroy()
    ptext2.destroy()
    tmove.destroy()
  # calculate Energy
  Ene = Calc_Ene(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim)
  mynewtext="step %d" % Iterations
  sttext0=Label(top,text=mynewtext)
  sttext0.pack(side='left')
  mynewtext="E=%.3f" % Ene
  entext1=Label(top,text=mynewtext)
  entext1.pack(side='left')
  mynewtext="Paccept= %.2f" % paccept
  ptext2=Label(top,text=mynewtext)
  ptext2.pack(side='left')
  mynewtext="   "
  tmove=Label(top,text=mynewtext)
  tmove.pack(side='left')
  ACCEPTED = 0
  frac_simple_move=1-frac_swap
  calc_elec=1
  if qat == 0: 
    frac_simple_move=1.0
    calc_elec=0
  while (not ACCEPTED):
   Iterations += 1
   # select a coordinate randomly
   RANDOM_atom  = randint(0,len(Atom_Coord)-1)
   RANDOM_coord = randint(0,1)
   EneSingle = Calc_EneSingle(RANDOM_atom,Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim,calc_elec)
   # save old coordinates of that atom
   Xold = Atom_Coord[RANDOM_atom][0]
   Yold = Atom_Coord[RANDOM_atom][1]
   xx=random()
   if xx<frac_simple_move:
    type = "move"
    # move the particule
    factor = ((2*random() - 1) * deltaRmax)
    Atom_Coord[RANDOM_atom][RANDOM_coord] += factor
    # apply boudary conditions
    for i in range(2): # i=0 -> case x coordinate ; i=1 -> case y coordinate
      halfbox = BoxDim[i]/2
      Z = Atom_Coord[RANDOM_atom][RANDOM_coord]
      Atom_Coord[RANDOM_atom][RANDOM_coord] = Z - SignR(halfbox,Z) - SignR(halfbox,Z-BoxDim[i])
    EneSingle_new = Calc_EneSingle(RANDOM_atom,Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim,calc_elec)
    deltaE = EneSingle_new - EneSingle
    if deltaE < 0.0:
      proba_boltzmann = 1.0
      ACCEPTED = 1 ; Accepted += 1 ; Ene = Ene + deltaE
    else:
      proba_boltzmann = exp(-deltaE/(cstboltz*Temperature))
      xx = random()
      if proba_boltzmann > xx:
        ACCEPTED = 1 ; Accepted += 1 ; Ene = Ene + deltaE
      else:
        ACCEPTED = 0
        # get back the orginal coordinate
        Atom_Coord[RANDOM_atom][RANDOM_coord] -= factor
    paccept = 0.01*ACCEPTED+0.99*paccept;
    if Iterations > 0: 
      sttext0.destroy()
      entext1.destroy()
      ptext2.destroy()
      tmove.destroy()
    mynewtext="step %d" % Iterations
    sttext0=Label(top,text=mynewtext)
    sttext0.pack(side='left')
    mynewtext="E= %.3f" % Ene
    entext1=Label(top,text=mynewtext)
    entext1.pack(side='left')
    mynewtext="Paccept= %.2f" % paccept
    ptext2=Label(top,text=mynewtext)
    ptext2.pack(side='left')
    mynewtext="move"
    tmove=Label(top,text=mynewtext)
    tmove.pack(side='left')
   else:
    type="swap"
    RANDOM_atom2  = randint(0,len(Atom_Coord)-1)
    Atom_Coord[RANDOM_atom][0] = Atom_Coord[RANDOM_atom2][0]
    Atom_Coord[RANDOM_atom][1] = Atom_Coord[RANDOM_atom2][1]
    Atom_Coord[RANDOM_atom2][0] = Xold
    Atom_Coord[RANDOM_atom2][1] = Yold
    Ene_new = Calc_Ene(Atom_Coord,Epsilon,Rmin,Dielec,CutOffSquare,BoxDim,calc_elec)
    if Ene_new < Ene:
      proba_boltzmann = 1.0
      ACCEPTED = 1 ; Accepted += 1 ; Ene = Ene_new
    else:
      deltaE = Ene_new - Ene
      proba_boltzmann = exp(-deltaE/(cstboltz*Temperature))
      xx = random()
      if proba_boltzmann > xx:
        ACCEPTED = 1 ; Accepted += 1 ; Ene = Ene_new
      else:
        ACCEPTED = 0
        # get back the orginal coordinate
        Atom_Coord[RANDOM_atom2][0] = Atom_Coord[RANDOM_atom][0]
        Atom_Coord[RANDOM_atom2][1] = Atom_Coord[RANDOM_atom][1]
        Atom_Coord[RANDOM_atom][0] = Xold
        Atom_Coord[RANDOM_atom][1] = Yold
    paccept = 0.01*ACCEPTED+0.99*paccept;
    if Iterations > 0: 
      sttext0.destroy()
      entext1.destroy()
      ptext2.destroy()
      tmove.destroy()
    mynewtext="step %d" % Iterations
    sttext0=Label(top,text=mynewtext)
    sttext0.pack(side='left')
    mynewtext="E= %.3f" % Ene
    entext1=Label(top,text=mynewtext)
    entext1.pack(side='left')
    mynewtext="Paccept= %.2f" % paccept
    ptext2=Label(top,text=mynewtext)
    ptext2.pack(side='left')
    mynewtext="swap"
    tmove=Label(top,text=mynewtext)
    tmove.pack(side='left')
   if (float(Iterations)%50.0) == 0.0:  
     print("step %6i E= %f  Paccepted= %5.3f Boltzmann-factor= %10.8f random= %6.4f" % (Iterations,Ene,(Accepted/float(Iterations)),proba_boltzmann,xx))
   elif (proba_boltzmann < 1.0 and proba_boltzmann > xx):
     print("step %6i E= %f  Paccepted= %5.3f Boltzmann-factor= %10.8f random= %6.4f" % (Iterations,Ene,(Accepted/float(Iterations)),proba_boltzmann,xx))
  if type=="move":
   atoms_to_do=[RANDOM_atom]
  elif type=="swap":
   atoms_to_do=[RANDOM_atom,RANDOM_atom2]
  # update in the canvas
  # delete old
  for ii in atoms_to_do:
   canevas.delete(ATOM[ii])
  # draw new
  # facteur de conversion
   x1 = Atom_Coord[ii][0] + Radius
   y1 = Atom_Coord[ii][1] + Radius
   x2 = Atom_Coord[ii][0] - Radius
   y2 = Atom_Coord[ii][1] - Radius
   ATOM[ii] = Canvas.create_oval(canevas,x1,y1,x2,y2,fill=Color[ii])
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
 global Iterations
 ATOM = [] # liste contenant les widgets

 if (Repack==1): 
     Atom_Coord = InitConf(nAtoms,BoxDim,Radius,qat,frac_neg)
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
 return Atom_Coord,ATOM,Color

def set_r(event):
 global Atom_Coord,ATOM,Color,canevas
 global r,nAtoms
 nAtoms=int(r.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color)
 canevas.pack()
 
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
 global r,nAtoms,temp,Temperature
 Temperature=float(temp.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 canevas.pack()

def set_dmv(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,nAtoms,dmv,deltaRmax
 deltaRmax =int(dmv.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 canevas.pack()

def set_swp(event):
 global Atom_Coord,ATOM,Color,canevas,frac_swap
 global r,nAtoms,dmv,deltaRmax
 frac_swap=float(swp.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 canevas.pack()

def set_spd(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,nAtoms,spd,speed
 speed=max(1,int(spd.get()))
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color,0)
 canevas.pack()
###########################################################################
###########################################################################

################
# MAIN PROGRAM #
################

from tkinter import *
from tkinter import Canvas

seed(Seed)
root = Tk() #root (main) window

canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
canevas.bind("<Button-1>",Go)
root.bind("<Escape>",die)
top2=Frame(root)
top2.pack(side='top')
top=Frame(root)
top.pack(side='top')
low=Frame(root)
low.pack(side='bottom')
low2=Frame(root)
low2.pack(side='bottom')
low3=Frame(root)
low3.pack(side='bottom')
mynewtext="Charged particles Monte Carlo"
hwtext=Label(top2,text=mynewtext,foreground='red',font='times 18 bold')
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
dmv=DoubleVar()
swp=DoubleVar()

Color = []
ATOM = [] # liste contenant les widgets
Atom_Coord = InitConf(nAtoms,BoxDim,Radius,qat,frac_neg)

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

#mynewtext3="Rmin = " 
#hwtext3=Label(low3,text=mynewtext3)
#hwtext3.pack(side='left')
#vdw2.set(Rmin)
#vdw2_entry=Entry(low3,width=6,textvariable=vdw2)
#vdw2_entry.pack(side='left')
#vdw2_entry.bind('<Return>', set_vdw2)


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

mynewtext3="Dmove = " 
hwtext3=Label(low,text=mynewtext3)
hwtext3.pack(side='left')
dmv.set(deltaRmax)
dmv_entry=Entry(low,width=6,textvariable=dmv)
dmv_entry.pack(side='left')
dmv_entry.bind('<Return>', set_dmv)

mynewtext3="frac charge swap = " 
hwtext3=Label(low,text=mynewtext3)
hwtext3.pack(side='left')
swp.set(frac_swap)
swp_entry=Entry(low,width=6,textvariable=swp)
swp_entry.pack(side='left')
swp_entry.bind('<Return>', set_swp)

startbutton=Button(top,text=' start ',command=Go,background='yellow',foreground='blue',relief='flat')
startbutton.pack(side='left',fill='x')
stopbutton=Button(top,text=' stop ',command=stop,background='yellow', foreground='red',relief='flat')
stopbutton.pack(side='left')

mynewtext3="Delay = " 
hwtext3=Label(low,text=mynewtext3)
hwtext3.pack(side='left')
spd.set(speed)
spd_entry=Entry(low,width=6,textvariable=spd)
spd_entry.pack(side='left')
spd_entry.bind('<Return>', set_spd)

Atom_Coord,ATOM,Color=setupall(Atom_Coord,ATOM,Color)

# set up a global var with the nb of iterations and nb of accepted conf
Iterations = 0 ; Accepted = 0; paccept = 1
# set up Ene as a global variable
Ene = 0.0
counter=0
# adapter a la dimension de la fenetre
canevas.pack()

# calculate first energy
#Calc_Ene(Atom_Coord)

print("Click on mouse button within the box to go ahead !!!")
print("<ESC> to quit")

# conserver la fenetre ouverte (inutile en interactif)
root.mainloop()
