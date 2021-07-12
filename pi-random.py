#!/usr/local/bin/python 
#
# Monte Carlo estimation of pi by counting the number of randomly generated
# points within a circle
# Alexandre Bonvin, Aalt Jan van Dijk, Utrecht University
#
##################
# import modules #
##################
from math import sqrt
from random import random,randint

####################
# define constants #
####################

BoxDim = [500.0,500.0] # box dimension
Center = [250.0,250.0] # coordinates of center
R2 = 250.0**2          # radius squared
Radius = 0             # radius of points in pixel
speed = 1              # canvas update speed

##################
# some functions #
##################

def dist2(x,y,B):
  return (x-B[0])**2+(y-B[1])**2

#########################
# initialize parameters #
#########################

#####################
# move particules   #  
#####################
def Go(*args):
  global Radius,sttext0,counter,lcounter,ncircle,pi,pi_old,tpi,frac_change,error
  global sttext1,pi_bailey,pi_bailey_old,wfrac_change,werror
  counter=counter+1.0
  lcounter=lcounter+1.0

  x = int(random()*BoxDim[0])
  y = int(random()*BoxDim[1])
  
  pi_old = pi
  pi_bailey_old = pi_bailey
  
  if dist2(x,y,Center)<=R2:
    ncircle = ncircle + 1.0
    pi = 4*ncircle/counter
    frac_change=100.0*abs(pi-pi_old)/pi
    error=100.0*abs(tpi-pi)/tpi
  
# and now calculate it using bailey's sum:
# Pi = SUMk=0 to infinity 16-k  [ 4/(8k+1) - 2/(8k+4) - 1/(8k+5) - 1/(8k+6) ].
  
  pi_bailey = pi_bailey + 16**(-counter) * (4.0/(8.0*counter+1.0) - 2.0/(8*counter+4.0) - 1.0/(8.0*counter+5.0) - 1.0/(8.0*counter+6.0))
  if counter < 50.0: 
    print("Step ",counter," PI-Bailey= ",pi_bailey,"PI-MC= ",pi," Error= ",error)

  wfrac_change=100.0*abs(pi_bailey-pi_bailey_old)/pi_bailey
  werror=100.0*abs(tpi-pi_bailey)/tpi
  
  x1 = x + Radius
  y1 = y + Radius
  x2 = x - Radius
  y2 = y - Radius
  xr = Canvas.create_oval(canevas,x1,y1,x2,y2,fill="#FFFFFF")
  if lcounter==50.0:
    print("Step ",counter," PI-Bailey= ",pi_bailey,"PI-MC= ",pi," Error= ",error)
    if counter>50.0:
      sttext0.destroy()
      sttext1.destroy()
    lcounter=0
    mynewtext="step: %10d PI= %11.8f %%error %8.5f %%change %8.5f" % (counter,pi,error,frac_change)
    sttext0=Label(low,text=mynewtext)
    sttext0.pack(side='left')
    mynewtext="Bailey's sum PI= %11.8f %%error %8.5f %%change %8.5f" % (pi_bailey,werror,wfrac_change)
    sttext1=Label(low2,text=mynewtext)
    sttext1.pack(side='left')
  canevas.coords(xr,x1,y1,x2,y2)
  canevas.after(speed,Go)
  
###########################################################################
# some functions added for extra buttons etc.
###########################################################################
def die(event=0):
  import sys
  sys.exit()

################
# MAIN PROGRAM #
################

from tkinter import *
from tkinter import Canvas
from math import acos

root = Tk() #root (main) window
spd=DoubleVar()

canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
canevas.bind("<Button-1>",Go)
root.bind("<Escape>",die)
top=Frame(root)
top.pack(side='top')
low=Frame(root)
low.pack(side='bottom')
low2=Frame(root)
low2.pack(side='bottom')

tpi=acos(-1)

mynewtext="Monte Carlo estimation of PI \n     %10.8f" % (tpi)
hwtext=Label(top,text=mynewtext,foreground='red',font='times 18 bold')
hwtext.pack(side='left')

counter=0.0
ncircle=0.0
pi=0.0
pi_old=0.0
pi_bailey = (4.0 - 0.5 - 1.0/5.0 - 1.0/6.0)

frac_change=0.0
error=0.0
lcounter=0
canevas.pack()

x1 = BoxDim[0]
y1 = BoxDim[1]
x2 = 0
y2 = 0
xr = Canvas.create_oval(canevas,x1,y1,x2,y2,fill="#FFFFFF")
canevas.coords(xr,x1,y1,x2,y2)


print("Click on mouse button within the box to go ahead !!!")
print("<ESC> to quit")

root.mainloop()
