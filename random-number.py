#!/usr/local/bin/python 
#
# Check your random number generator by plotting points within a square
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
Radius = 1         # radius of points in pixel
speed = 1          # canvas update speed

##################
# some functions #
##################

def dist(A,B):
  return sqrt((A[0]-B[0])**2+(A[1]-B[1])**2)

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

  x = int(random()*BoxDim[0])
  y = int(random()*BoxDim[1])
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
