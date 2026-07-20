#!/usr/local/bin/python 

##################
# import modules #
##################
from math import sqrt,exp
from random import random,randint

####################
# define constants #
####################

Radius = 25.0      	# radius used to draw
BoxDim = [500,500] 	# box dimension
Atom_Coord = []    	# list of the form : [NMAX][2]
deltaRmax = 500.0   	# step up to which particules move
b1 = 100.0     		# well depth
b2 = 120.0     		# well depth
a1 = 100.0     		# well length
a2 = 50.0     		# well length
l = 100.0     		#distance used in drawing in-between the two wells
speed = 50         	# canvas update speed
cstboltz = 8.34E-3 	# Boltzmann's contstant in kcal/mol/K
Temperature = 300.0 	# temperature in K
startx,starty=170,300 	#position to start
scaley=10		#scalefactor

#########################
# initialize parameters #
#########################

# define the double well
def InitConf(dim,a1,a2,b1,b2,l):
  print("Initializing box, please wait...")
  tmp_coord = [[startx,starty]]
  x0=startx-a1
  y0=10
  x1=startx-a1
  y1=starty
  x2=startx+a1
  y2=starty
  x3=startx+a1
  y3=10
  x4=startx+a1+l
  y4=10
  x5=startx+a1+l
  y5=starty+scaley*(-b1+b2)
  x6=startx+a1+l+2*a2
  y6=starty+scaley*(-b1+b2)
  x7=startx+a1+l+2*a2
  y7=10
  pot_coord=[[x0,y0],[x1,y1],[x2,y2],[x3,y3],[x4,y4],[x5,y5],[x6,y6],[x7,y7]]
  return tmp_coord,pot_coord

####################
# calculate energy #
####################
def Calc_Ene(coord,a1,a2,b1,b2):
  Ene = 0.0 
  big=1E99
  #big is meant to represent infinity
  x=coord[0][0]
  if (x>(startx-a1)) and (x<(startx+a1)):
   #within well 1
   Ene=-1.0*b1
  elif (x>(startx+a1+l)) and (x<(startx+a1+l+2*a2)):
   #within well 2
   Ene=-1.0*b2
  else:
   Ene=big
   
  return Ene

############################
# move particules in a MC  #  
############################
def Go(*args):
  global Atom_Coord,Radius,deltaRmax,Iterations,Ene,Accepted
  global Color,sttext0,entext1,ptext2,ptext3,paccept,canevas,nr1,nr2

  RANDOM_atom  = 0
  RANDOM_coord = 0

  # ONE Metropolis trial move per step (proper residence-time weighting: a
  # rejected move keeps the current configuration, which is still counted below).
  Iterations += 1
  # save the old coordinates of that atom
  Xold = Atom_Coord[RANDOM_atom][0]
  Yold = Atom_Coord[RANDOM_atom][1]
  # energy of the current position (depends only on x -> recompute, no stale state)
  Eold = Calc_Ene(Atom_Coord,a1,a2,b1,b2)
  # propose a random displacement
  factor = ((2*random() - 1) * deltaRmax)
  Atom_Coord[RANDOM_atom][RANDOM_coord] += factor
  # set the drawing height according to which well the particle now sits in
  if ((Atom_Coord[RANDOM_atom][RANDOM_coord])> (startx-a1)) and ((Atom_Coord[RANDOM_atom][RANDOM_coord])< (startx+a1)):
    Atom_Coord[RANDOM_atom][1]=starty
  elif ((Atom_Coord[RANDOM_atom][RANDOM_coord])> (startx+a1+l)) and ((Atom_Coord[RANDOM_atom][RANDOM_coord])< (startx+a1+l+2*a2)):
    Atom_Coord[RANDOM_atom][1]=starty+scaley*(-b1+b2)
  Ene_new = Calc_Ene(Atom_Coord,a1,a2,b1,b2)
  deltaE = Ene_new - Eold
  # Metropolis acceptance criterion
  if deltaE < 0.0:
    ACCEPTED = 1
  else:
    test = -deltaE/(cstboltz*Temperature)
    proba_boltzmann = exp(test) if test > -100 else 0.0
    ACCEPTED = 1 if proba_boltzmann > random() else 0
  if ACCEPTED:
    Accepted += 1
    Ene = Ene_new
  else:
    # reject: restore the original coordinates (x and drawing height)
    Atom_Coord[RANDOM_atom][0] = Xold
    Atom_Coord[RANDOM_atom][1] = Yold
    Ene = Eold
  paccept = 0.01*ACCEPTED+0.99*paccept

  # count the current well every step (accepted or rejected) -> Boltzmann populations
  xnow = Atom_Coord[RANDOM_atom][RANDOM_coord]
  if (xnow > (startx-a1)) and (xnow < (startx+a1)):
    nr1=nr1+1
  elif (xnow > (startx+a1+l)) and (xnow < (startx+a1+l+2*a2)):
    nr2=nr2+1

  # update the on-screen labels
  if Iterations > 1:
    sttext0.destroy()
    entext1.destroy()
    ptext2.destroy()
    ptext3.destroy()
  sttext0=Label(top,text="step %d" % Iterations)
  sttext0.pack(side='left')
  entext1=Label(top,text="E= %.3f" % Ene)
  entext1.pack(side='left')
  ptext2=Label(top,text="Paccept= %.2f" % paccept)
  ptext2.pack(side='left')
  frac1=float(nr1)/(float(nr1)+float(nr2))
  ptext3=Label(top,text="fraction in well1= %.2f" % frac1)
  ptext3.pack(side='left')

  # update the atom in the canvas
  canevas.delete(ATOM[RANDOM_atom])
  x1 = Atom_Coord[RANDOM_atom][0] + Radius
  y1 = Atom_Coord[RANDOM_atom][1] + Radius
  x2 = Atom_Coord[RANDOM_atom][0] - Radius
  y2 = Atom_Coord[RANDOM_atom][1] - Radius
  ATOM[RANDOM_atom] = canevas.create_oval(x1,y1,x2,y2,fill=Color[RANDOM_atom])
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
 global Atom_Coord,Pot_Coord,ATOM,Color,canevas,POT
 global r,Radius
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color,POT=setupall(Atom_Coord,ATOM,Color,0)
 canevas.pack()

def setupall(Atom_Coord,ATOM,Color,Repack=1):
 global Iterations,Pot_Coord,nr1,nr2
 ATOM = [] # liste contenant les widgets
 POT=[]

 if (Repack==1): 
     nr1=1
     nr2=0
     Atom_Coord, Pot_Coord = InitConf(BoxDim,a1,a2,b1,b2,l)
     Color = ["orange"]
 for i in range(len(Atom_Coord)):
  x1 = Atom_Coord[i][0] + Radius
  y1 = Atom_Coord[i][1] + Radius
  x2 = Atom_Coord[i][0] - Radius
  y2 = Atom_Coord[i][1] - Radius
  ATOM.append(canevas.create_oval(x1,y1,x2,y2,fill=Color[i]))
 for i in range(len(Pot_Coord)-1):
  x1 = Pot_Coord[i][0]
  y1 = Pot_Coord[i][1]
  x2 = Pot_Coord[i+1][0]
  y2 = Pot_Coord[i+1][1]
  POT.append(canevas.create_line(x1,y1,x2,y2,fill="red"))
 
 return Atom_Coord,ATOM,Color,POT

def set_a1(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,a1
 a1=int(pota1.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color,POT=setupall(Atom_Coord,ATOM,Color)
 canevas.pack()
 
def set_a2(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,a2
 a2=int(pota2.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color,POT=setupall(Atom_Coord,ATOM,Color)
 canevas.pack()
 
def set_b1(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,b1
 b1=int(potb1.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color,POT=setupall(Atom_Coord,ATOM,Color)
 canevas.pack()
 
def set_b2(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,b2
 b2=int(potb2.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color,POT=setupall(Atom_Coord,ATOM,Color)
 canevas.pack()

def set_temp(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,temp,Temperature
 Temperature=float(temp.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color,POT=setupall(Atom_Coord,ATOM,Color)
 canevas.pack()

def set_dmv(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,dmv,deltaRmax
 deltaRmax =int(dmv.get())
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color,POT=setupall(Atom_Coord,ATOM,Color,0)
 canevas.pack()

def set_spd(event):
 global Atom_Coord,ATOM,Color,canevas,Radius
 global r,spd,speed
 speed=max(1,int(spd.get()))
 canevas.destroy()
 canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
 canevas.bind("<Button-1>",Go)
 Atom_Coord,ATOM,Color,POT=setupall(Atom_Coord,ATOM,Color,0)
 canevas.pack()
###########################################################################
###########################################################################

################
# MAIN PROGRAM #
################

from tkinter import *
#from Canvas import *

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
mynewtext="Double well Monte Carlo"
hwtext=Label(top2,text=mynewtext,foreground='red',font='times 18 bold')
hwtext.pack(side='left')

r = DoubleVar()
pota1=DoubleVar()
pota2=DoubleVar()
potb1=DoubleVar()
potb2=DoubleVar()
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
Atom_Coord, Pot_Coord = InitConf(BoxDim,a1,a2,b1,b2,l)

mynewtext3="a1 = " 
hwtext3=Label(low3,text=mynewtext3)
hwtext3.pack(side='left')
pota1.set(a1)
pota1_entry=Entry(low3,width=6,textvariable=pota1)
pota1_entry.pack(side='left')
pota1_entry.bind('<Return>', set_a1)

mynewtext3="a2 = " 
hwtext3=Label(low3,text=mynewtext3)
hwtext3.pack(side='left')
pota2.set(a2)
pota2_entry=Entry(low3,width=6,textvariable=pota2)
pota2_entry.pack(side='left')
pota2_entry.bind('<Return>', set_a2)

mynewtext3="b1 = " 
hwtext3=Label(low3,text=mynewtext3)
hwtext3.pack(side='left')
potb1.set(b1)
potb1_entry=Entry(low3,width=6,textvariable=potb1)
potb1_entry.pack(side='left')
potb1_entry.bind('<Return>', set_b1)

mynewtext3="b2 = " 
hwtext3=Label(low3,text=mynewtext3)
hwtext3.pack(side='left')
potb2.set(b2)
potb2_entry=Entry(low3,width=6,textvariable=potb2)
potb2_entry.pack(side='left')
potb2_entry.bind('<Return>', set_b2)

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

Atom_Coord,ATOM,Color,POT=setupall(Atom_Coord,ATOM,Color)

# set up a global var with the nb of iterations and nb of accepted conf
Iterations = 0 ; Accepted = 0; paccept = 1 
nr1=1; nr2=0 #number of times in well1 and in well2
# set up Ene as a global variable
Ene = 0.0
counter=0

canevas.pack()

print("Click on mouse button within the box to go ahead !!!")
print("<ESC> to quit")

# start the main loop
root.mainloop()
