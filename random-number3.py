#!/usr/local/bin/python 

##################
# import modules #
##################
from math import sqrt
from random import random,randint

####################
# define constants #
####################

BoxDim = [500.0,500.0] # box dimension
Radius = 1         # radius of points in pixel
speed = 1          # canvas update speed
xn = 11            # seed for random number generator

##################
# some functions #
##################

def dist(A,B):
  return sqrt((A[0]-B[0])**2+(A[1]-B[1])**2)
  
# Random number generator based on:
# http://crypto.mat.sbg.ac.at/results/karl/server/node3.html

def myrandom():
 global xn
 #LCG: x_n+1=a*x_n + b (mod m); seed x0
 #LCG(m,a,b,x0)
 #bad: LCG(2^32,477211307,0,1)
 a=156778911
 b=3
 m=2E28
 xn=a*xn+b
 xn=xn%m
 val=xn/m
 return val

#########################
# initialize parameters #
#########################

#####################
# move particules   #  
#####################
def Go(*args):
  global Radius,sttext0,counter,lcounter
  counter=counter+1
  lcounter=lcounter+1

  #x = int(random()*BoxDim[0])
  #y = int(random()*BoxDim[1])
  x = int(myrandom()*BoxDim[0])
  y = int(myrandom()*BoxDim[1])
  x1 = x + Radius
  y1 = y + Radius
  x2 = x - Radius
  y2 = y - Radius
  xr = Canvas.create_oval(canevas,x1,y1,x2,y2,fill="#FFFFFF")
  if lcounter==50:
    if counter>50:
      sttext0.destroy()
    lcounter=0
    mynewtext="step: %10d" % (counter)
    sttext0=Label(low,text=mynewtext)
    sttext0.pack(side='left')
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

root = Tk() #root (main) window
spd=DoubleVar()

canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1],bg="#ccddff")
canevas.bind("<Button-1>",Go)
root.bind("<Escape>",die)
top=Frame(root)
top.pack(side='top')
low=Frame(root)
low.pack(side='bottom')

mynewtext="How random is your random number generator?"
hwtext=Label(top,text=mynewtext,foreground='red',font='times 18 bold')
hwtext.pack(side='left')

counter=0
lcounter=0
canevas.pack()

print("Click on mouse button within the box to go ahead !!!")
print("<ESC> to quit")


root.mainloop()
