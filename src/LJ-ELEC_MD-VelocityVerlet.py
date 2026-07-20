#!/usr/local/bin/python
#
###########################################################################
###########################################################################
# Simple MD of Lennard Jones charged or uncharges particles
# Velocity-Verlet integrator + Berendsen weak-coupling thermostat
# Alexandre Bonvin, Aalt Jan van Dijk, Utrecht University
#
# adapted from a script from Patrick Fuchs, Uni. Paris VI
###########################################################################
###########################################################################
#
##################
# import modules #
##################
from math import sqrt, exp, log, sin, cos
from random import random, randint, seed

###########################################################################
###########################################################################
### define parameters #####################################################
###########################################################################
###########################################################################

nAtoms = 20  # number of atoms
Radius = 25.0  # beware that Radius must be in a good range (according to nAtoms)
# in order to be able to place all atoms
Mass = 10.0  # atom mass
Rmin = 2.24 * Radius  # distance at which rmin is mini
BoxDim = [500, 500]  # box dimension
Atom_Coord = []  # list of the form : [NMAX][2]
Epsilon = 25.0  # well depth
Dielec = 1.0  # dielectric constant
qat = Radius  # Atom absolute charge
frac_neg = 0.5  # Fraction negative charges
OverlapFr = 0.0  # fraction of overlap allowed
CutOff = 250  # non-bonded cutoff
CutOffSquare = CutOff**2
# soft core: over a long run an oppositely-charged pair can drift into hard contact and
# the r^-12 core blows up at this timestep. Clamping the pair distance to >= SoftCore in
# the energy/force evaluation caps the repulsion for deeply overlapping pairs and keeps
# the integration stable (a rare collision becomes a recoverable bounce, not an explosion).
SoftCore = 0.8 * Rmin  # soft-core distance
SoftCoreSquare = SoftCore**2
speed = 20  # canvas update speed
cstboltz = 0.00198722  # Boltzmann's constant in kcal/mol/K
cstboltz = 1000 * cstboltz / 4.18  # in kJ/mol/K
Temperature = 300.0  # target temperature in K
timestep = 1.0e-3  # MD time step
tauT = 0.1  # thermostat coupling time (same units as timestep); larger = weaker coupling
Seed = 100  # random number seed


###########################################################################
###########################################################################
# Velocity-Verlet integrator  #############################################
###########################################################################
###########################################################################
#
# The velocity-Verlet scheme propagates BOTH positions and velocities and is
# self-starting (no need for a special first step or for the previous positions).
# Writing a = F/m for the acceleration, one MD step from time t to t+h reads:
#
#   (1) r(t+h) = r(t) + v(t) h + 1/2 a(t) h^2
#   (2) evaluate the new forces F(t+h) at the new positions r(t+h)
#   (3) v(t+h) = v(t) + 1/2 [ a(t) + a(t+h) ] h
#
# So each step needs the OLD forces (from the previous step) to move the positions,
# and the NEW forces to finish the velocity update. We therefore keep the current
# forces around between steps (global Force).


def VelVerlet_Pos(atom_coord, velocity, force, h, mass):
    # step (1): new positions from current positions, velocities and forces
    #
    # A few hints:
    # - powers in python are given by **, e.g.: x to the square is x**2
    # - indents are important in python
    list = []
    for i in range(len(atom_coord)):
        q = atom_coord[i][2]
        r0x = atom_coord[i][0]  # coordinates
        r0y = atom_coord[i][1]
        v0x = velocity[i][0]  # velocities
        v0y = velocity[i][1]

        # Insert below the lines defining the new x and y positions using the
        # velocity-Verlet position update r + v*h + 0.5*(F/m)*h**2
        # >>>>>
        r0xnew = r0x + v0x * h + 0.5 * force[i][0] / mass * h**2
        r0ynew = r0y + v0y * h + 0.5 * force[i][1] / mass * h**2
        # <<<<<
        list.append([r0xnew, r0ynew, q])
    return list


def VelVerlet_Vel(velocity, force_old, force_new, h, mass):
    # step (3): new velocities from the average of the old and new forces
    list = []
    for i in range(len(velocity)):
        v0x = velocity[i][0]
        v0y = velocity[i][1]

        # Insert below the lines defining the new x and y velocities using the
        # velocity-Verlet velocity update v + 0.5*(F_old + F_new)/m*h
        # >>>>>
        v0xnew = v0x + 0.5 * (force_old[i][0] + force_new[i][0]) / mass * h
        v0ynew = v0y + 0.5 * (force_old[i][1] + force_new[i][1]) / mass * h
        # <<<<<
        list.append([v0xnew, v0ynew])
    return list


###########################################################################
###########################################################################
# Berendsen weak-coupling thermostat  #####################################
###########################################################################
###########################################################################
#
# The system is coupled to an external heat bath at the target temperature T0.
# Each step the velocities are rescaled by a factor lambda so that the temperature
# relaxes towards T0 with a time constant tau:
#
#   lambda = sqrt( 1 + (h/tau) ( T0/T - 1 ) )
#
# - tau = h  -> lambda = sqrt(T0/T): the temperature is reset to T0 every step
#               (strong coupling / simple velocity rescaling)
# - tau >> h -> lambda ~ 1: very weak coupling, the temperature drifts slowly to T0
# The kinetic energy (hence T) is only gently nudged, so the dynamics stay smooth.


def Berendsen(velocity, temp_current, temp_target, h, tau):
    if temp_current <= 0.0:
        return velocity
    # weak-coupling velocity scaling factor
    fac = 1.0 + (h / tau) * (temp_target / temp_current - 1.0)
    if fac < 0.0:
        fac = 0.0
    lam = sqrt(fac)
    return [[v[0] * lam, v[1] * lam] for v in velocity]


###########################################################################
###########################################################################
# move particules with MD  ################################################
###########################################################################
###########################################################################


def Go(*args):
    global \
        Atom_Coord, \
        Radius, \
        Mass, \
        BoxDim, \
        Epsilon, \
        Rmin, \
        CutOffSquare, \
        Iterations, \
        Ene, \
        Force
    global Velocity, timestep, tauT, hwtextene1, hwtextene2
    global Color, sttext0, ptext2, paccept, Dielec, root, canevas, speed, nout
    hwtextene1.destroy()
    hwtextene2.destroy()
    sttext0.destroy()

    # On the very first step (or after a reset) compute the forces at the current
    # positions -- velocity-Verlet needs a(t) before it can move the atoms.
    if Iterations == 0:
        nout = 0
        Force = Calc_Force(Atom_Coord, Epsilon, Rmin, Dielec, CutOffSquare, BoxDim)

    # --- velocity-Verlet step ---
    # (1) new positions from the current forces
    Atom_Coord = VelVerlet_Pos(Atom_Coord, Velocity, Force, timestep, Mass)
    # (2) forces at the new positions
    Force_new = Calc_Force(Atom_Coord, Epsilon, Rmin, Dielec, CutOffSquare, BoxDim)
    # (3) new velocities from the average of old and new forces
    Velocity = VelVerlet_Vel(Velocity, Force, Force_new, timestep, Mass)
    Force = Force_new
    nout = nout + 1

    # --- temperature control by weak coupling (Berendsen thermostat) ---
    Kin, temperature = calc_temp(Velocity, nAtoms, cstboltz, Mass)
    Velocity = Berendsen(Velocity, temperature, Temperature, timestep, tauT)
    Kin, temperature = calc_temp(Velocity, nAtoms, cstboltz, Mass)

    # potential energy at the new positions
    Ene, EneLJ, EneCoul = Calc_Ene2(
        Atom_Coord, Epsilon, Rmin, Dielec, CutOffSquare, BoxDim
    )

    mynewtext = "step: %d Time: %8.3f" % (Iterations, float(Iterations) * timestep)
    sttext0 = Label(top1, text=mynewtext)
    sttext0.pack(side="left")

    mynewtext = "Etot: %6.1f Ekin: %6.1f Epot: %6.1f" % (Ene + Kin, Kin, Ene)
    hwtextene1 = Label(top2, text=mynewtext)
    hwtextene1.pack(side="left")

    mynewtext = "Elj: %6.1f Ecoul: %6.1f Temp: %6.1f" % (EneLJ, EneCoul, temperature)
    hwtextene2 = Label(top3, text=mynewtext)
    hwtextene2.pack(side="left")

    if temperature > 1000000:
        print("The system is exploding !!!")
        print("step: %d Time: %8.3f" % (Iterations, float(Iterations) * timestep))
        print("Etot: %6.1f Ekin: %6.1f Epot: %6.1f" % (Ene + Kin, Kin, Ene))
        print("Elj: %6.1f Ecoul: %6.1f Temp: %6.1f" % (EneLJ, EneCoul, temperature))
        print("Emergency stop")
        sys.exit()

    # apply boudary conditions (velocity-Verlet stores velocities explicitly, so
    # wrapping the positions does not corrupt anything -- no predecessor to shift)
    for pp in range(len(Atom_Coord)):
        for i in range(2):  # i=0 -> case x coordinate ; i=1 -> case y coordinate
            if Atom_Coord[pp][i] < 0:
                Atom_Coord[pp][i] += BoxDim[i]
            if Atom_Coord[pp][i] > BoxDim[i]:
                Atom_Coord[pp][i] -= BoxDim[i]

    # draw new canvas coordinates
    for i in range(len(Atom_Coord)):
        x1 = Atom_Coord[i][0] + Radius
        y1 = Atom_Coord[i][1] + Radius
        x2 = Atom_Coord[i][0] - Radius
        y2 = Atom_Coord[i][1] - Radius
        canevas.coords(ATOM[i], x1, y1, x2, y2)
    Iterations = Iterations + 1

    # print to terminal window
    if nout == 20:
        nout = 0
        print(
            "step: %d Time: %8.3f Etot: %6.1f Ekin: %6.1f Epot: %6.1f Elj: %6.1f Ecoul: %6.1f Temp: %6.1f"
            % (
                Iterations,
                float(Iterations) * timestep,
                Ene + Kin,
                Kin,
                Ene,
                EneLJ,
                EneCoul,
                temperature,
            )
        )

    canevas.after(speed, Go)


###########################################################################
###########################################################################
### energy functions  #####################################################
###########################################################################
###########################################################################


# classical LJ
def LJ(r, epsilon, rmin):
    return epsilon * ((rmin / r) ** 12 - (rmin / r) ** 6)


# calculate LJ from the squared distance
def LJ2(distsquare, epsilon, rmin_exp6):
    Z = (1 / distsquare) ** 3 * rmin_exp6
    return epsilon * Z * (Z - 1)


# classical Coulomb
def Coulomb(r, dielec, qa, qb):
    return qa * qb / (dielec * r)


# classical Coulomb from the squared distance
def Coulomb2(r, dielec, qa, qb):
    return qa * qb / (dielec * sqrt(r))


# Calculate energy Evdw + Ecoulomb (used squared distance)
# version with boundary conditions
def Calc_Ene2(coord, epsilon, rmin, dielec, cutoffsquare, boxdim, elec=1):
    Ene = 0.0
    distsquare = 0
    ELJ = 0.0
    ECoul = 0.0
    rmin_exp6 = rmin**6
    # doubly nested loop over all particule pairs
    for i in range(len(coord) - 1):
        for j in range(i + 1, len(coord)):
            # calculate the squared atomic distance
            distsquare = 0
            for k in range(2):
                tmp = coord[j][k] - coord[i][k]
                # chooses the nearest image
                halfbox = boxdim[k] / 2
                tmp = (
                    tmp - SignR(halfbox, tmp - halfbox) - SignR(halfbox, tmp + halfbox)
                )
                distsquare += tmp**2
            # compute vdw and Coulomb energy
            if distsquare < cutoffsquare:
                distsquare = max(distsquare, SoftCoreSquare)  # soft core
                qa = coord[i][2]
                qb = coord[j][2]
                vdw = LJ2(distsquare, epsilon, rmin_exp6)
                Ene += vdw
                ELJ += vdw
                if elec:
                    CC = Coulomb2(distsquare, dielec, qa, qb)
                    Ene += CC
                    ECoul += CC
    return Ene, ELJ, ECoul


# Calculate kinetic energy and temperature
def calc_temp(vel, nat, k, mass):
    v2 = 0.0
    for i in range(len(vel)):
        vx = vel[i][0]  # velocity
        vy = vel[i][1]
        v2 = v2 + vx**2 + vy**2
    nkt = 0.5 * mass * v2  # kinetic energy equals 0.5*m*v**2
    kin = v2 * 0.5 * mass
    temp = nkt / (nat * k)  # N*k*T=Kinetic Energy
    return kin, temp


###########################################################################
###########################################################################
### force functions  ######################################################
###########################################################################
###########################################################################


# force LJ (use squared distance)
def ForceLJ2(distsquare, epsilon, rmin_exp6, xi):
    rij = sqrt(distsquare)
    Z = (1 / distsquare) ** 3 * rmin_exp6
    dedz = epsilon * (2 * Z - 1)
    dzdr = rmin_exp6 * (-6.0 / rij ** (7.0))
    drdx = xi / rij
    return dedz * dzdr * drdx


# Force Coulomb (use squared distance)
def ForceCoulomb(distsquare, dielec, qa, qb, xi):
    rij = sqrt(distsquare)
    dedr = -1.0 * (qa * qb / dielec) * (1 / distsquare)
    drdx = xi / rij
    return dedr * drdx


# Calculate force from Evdw + Ecoulomb (uses squared distance)
def Calc_Force(coord, epsilon, rmin, dielec, cutoffsquare, boxdim):
    Force = []
    distsquare = 0
    rmin_exp6 = rmin**6
    # doubly nested loop over all particle pairs
    for i in range(len(coord)):
        tmpforce = [0.0, 0.0]
        for j in range(len(coord)):
            if not (i == j):
                # calculate the squared atomic distance
                distsquare = 0
                for k in range(2):
                    tmp = coord[j][k] - coord[i][k]
                    # chooses the nearest image
                    halfbox = boxdim[k] / 2
                    tmp = (
                        tmp
                        - SignR(halfbox, tmp - halfbox)
                        - SignR(halfbox, tmp + halfbox)
                    )
                    distsquare += tmp**2
                # compute vdw force
                if distsquare < cutoffsquare:
                    distsquare = max(distsquare, SoftCoreSquare)  # soft core
                    qa = coord[i][2]
                    qb = coord[j][2]
                    fflist = []
                    for k in range(2):
                        tmp = coord[j][k] - coord[i][k]
                        ff = ForceLJ2(distsquare, epsilon, rmin_exp6, tmp)
                        ff += ForceCoulomb(distsquare, dielec, qa, qb, tmp)
                        fflist.append(ff)
                    for k in range(2):
                        tmpforce[k] = tmpforce[k] + fflist[k]
        Force.append(tmpforce)
    return Force


###########################################################################
###########################################################################
### other functions  ######################################################
###########################################################################
###########################################################################


### distance ###
def dist(A, B):
    return sqrt((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2)


### squared distance ###
def dist2(A, B):
    return (A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2


### change sign ###
def SignR(a, b):
    if b > 0:
        return a
    else:
        return -a


### color particules based on charge ###
def charge_color(charge, qat):
    tmp = "#111111"
    if charge == qat:
        tmp = "#FFFFFF"
    else:
        tmp = "#333333"
    return tmp


### stop
def die(event=0):
    import sys

    sys.exit()


###########################################################################
###########################################################################
### initialization  #######################################################
###########################################################################
###########################################################################


### generates random coordinates ###
def InitConf(n, dim, radius, qat, frac_neg):
    seed(Seed)
    print("Initializing box, please wait...")
    # generate a list of random positions
    tmp_coord = []
    i = 0
    nneg = 0
    ntrial = 0
    # fix first atom
    x = random() * (dim[0] - 2 * radius) + radius  # dim[0]
    y = random() * (dim[1] - 2 * radius) + radius  # dim[1]
    nneg = int(float(n) * frac_neg)
    npos = n - nneg
    charge = -qat
    if npos == n:
        charge = qat
    i += 1
    if n == 2:
        tmp_coord.append([175, 300, charge])
    else:
        tmp_coord.append([x, y, charge])
    while i < nneg:
        x = random() * (dim[0] - 2 * radius) + radius  # dim[0]
        y = random() * (dim[1] - 2 * radius) + radius  # dim[1]
        # check wether the new particule ovelap an existing one
        OVERLAP = 1
        for j in range(i):
            if dist(tmp_coord[j], [x, y]) < (1 - OverlapFr) * 2 * radius:
                OVERLAP = 0
        if OVERLAP:
            charge = -qat
            if n == 2:
                tmp_coord.append([325, 300, charge])
            else:
                tmp_coord.append([x, y, charge])
            i += 1
        ntrial = ntrial + 1
        if ntrial > 100000:
            print("initialisation failed")
            print("==> reduce radius or number of atoms")
            sys.exit()
    while i < n:
        x = random() * (dim[0] - 2 * radius) + radius  # dim[0]
        y = random() * (dim[1] - 2 * radius) + radius  # dim[1]
        # check wether the new particule overlap an existing one
        OVERLAP = 1
        for j in range(i):
            if dist(tmp_coord[j], [x, y]) < (1 - OverlapFr) * 2 * radius:
                OVERLAP = 0
        if OVERLAP:
            charge = qat
            if n == 2:
                tmp_coord.append([325, 300, charge])
            else:
                tmp_coord.append([x, y, charge])
            i += 1
        ntrial = ntrial + 1
        if ntrial > 10**10:
            print("initialisation failed")
            print("==> reduce radius or number of atoms")
            sys.exit()
    return tmp_coord


# generates random charges
def InitCharge(n, dim, qat, frac_neg):
    global Atom_Coord
    print("Initializing charges, please wait...")
    i = 0
    nneg = 0
    nneg = int(float(n) * frac_neg)
    npos = n - nneg
    charge = -qat
    if npos == n:
        charge = qat
    Atom_Coord[i][2] = charge
    i += 1
    while i < nneg:
        charge = -qat
        Atom_Coord[i][2] = charge
        i += 1
    while i < n:
        charge = qat
        Atom_Coord[i][2] = charge
        i += 1


# generates initial velocities according to Maxwell distribution
def InitVel(n, temperature, cstboltz, mass):
    stdev = sqrt(cstboltz * temperature / mass)
    print("initializing velocities, please wait..")
    tmp_vel = []
    # for testing purposes:
    if n == 2:
        vel = [[0, 0], [0, 0]]
    else:
        i = 0
        while i < n:
            r1 = random()
            r2 = random()
            # generate random numbers according to Gaussian:
            x1 = sqrt(-2.0 * log(r1)) * cos(r2)
            x2 = sqrt(-2.0 * log(r1)) * sin(0.5 * r2)
            x1 = x1 * stdev
            x2 = x2 * stdev
            tmp_vel.append([x1, x2])
            i += 1
    # remove overall motion
    i = 0
    vxt = 0.0
    vyt = 0.0
    while i < n:
        vxt += tmp_vel[i][0]
        vyt += tmp_vel[i][1]
        i += 1
    i = 0
    while i < n:
        tmp_vel[i][0] = tmp_vel[i][0] - vxt / float(n)
        tmp_vel[i][1] = tmp_vel[i][1] - vyt / float(n)
        i += 1
    # scaling factor is used to get temperature exactly equal to desired temperature
    kin, tt = calc_temp(tmp_vel, n, cstboltz, mass)
    scaling = sqrt(temperature / tt)
    vel = []
    i = 0
    while i < n:
        vx, vy = tmp_vel[i][0], tmp_vel[i][1]
        vx = vx * scaling
        vy = vy * scaling
        vel.append([vx, vy])
        i += 1
    return vel


###########################################################################
###########################################################################
### various functions for input + layout ##################################
###########################################################################
###########################################################################


def stop():
    global Atom_Coord, ATOM, Color, canevas
    global r, nAtoms, size, Radius, drinit
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color = setupall(Atom_Coord, ATOM, Color, 0)
    update_ene()
    canevas.pack()
    update_ene()


### setup system ###
def setupall(Atom_Coord, ATOM, Color, Repack=1):
    global Iterations, Velocity, Temperature, Mass, cstboltz
    ATOM = []  # liste contenant les widgets
    if Repack == 1:
        Atom_Coord = InitConf(nAtoms, BoxDim, Radius, qat, frac_neg)
        Color = []
        for i in range(nAtoms):
            Color.append(charge_color(Atom_Coord[i][2], qat))
        Velocity = InitVel(nAtoms, Temperature, cstboltz, Mass)
    if Repack == 2:
        InitCharge(nAtoms, BoxDim, qat, frac_neg)
        Color = []
        for i in range(nAtoms):
            Color.append(charge_color(Atom_Coord[i][2], qat))
    for i in range(len(Atom_Coord)):
        x1 = Atom_Coord[i][0] + Radius
        y1 = Atom_Coord[i][1] + Radius
        x2 = Atom_Coord[i][0] - Radius
        y2 = Atom_Coord[i][1] - Radius
        ATOM.append(canevas.create_oval(x1, y1, x2, y2, fill=Color[i]))
        update_ene()
    return Atom_Coord, ATOM, Color


### set number of particules ###
def set_r(event):
    global Atom_Coord, ATOM, Color, canevas
    global r, nAtoms, Iterations
    nAtoms = int(r.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color = setupall(Atom_Coord, ATOM, Color)
    update_ene()
    canevas.pack()
    Iterations = 0


### set atom Radius ###
def set_size(event):
    global Atom_Coord, ATOM, Color, canevas, Radius, Rmin
    global r, nAtoms, size
    Radius = int(size.get())
    Rmin = 2 * Radius
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color = setupall(Atom_Coord, ATOM, Color)
    update_ene()
    canevas.pack()


### set epsilon for Lennard-Jones ###
def set_vdw1(event):
    global Atom_Coord, ATOM, Color, canevas, Radius
    global r, nAtoms, vdw1, Epsilon
    Epsilon = int(vdw1.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color = setupall(Atom_Coord, ATOM, Color, 0)
    update_ene()
    canevas.pack()


### set sigma for Lennard-Jones ###
def set_vdw2(event):
    global Atom_Coord, ATOM, Color, canevas, Radius
    global r, nAtoms, vdw2, Rmin
    Rmin = int(vdw2.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color = setupall(Atom_Coord, ATOM, Color, 0)
    update_ene()
    canevas.pack()


### set charge fraction ###
def set_frac(event):
    global Atom_Coord, ATOM, Color, canevas, frac_neg
    global r, nAtoms
    frac_neg = float(frac.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color = setupall(Atom_Coord, ATOM, Color, 2)
    update_ene()
    canevas.pack()


### set particule charge ###
def set_q(event):
    global Atom_Coord, ATOM, Color, canevas, qat
    global r, nAtoms
    qat = float(q.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color = setupall(Atom_Coord, ATOM, Color, 2)
    update_ene()
    canevas.pack()


### set dielectric constant ###
def set_diel(event):
    global Atom_Coord, ATOM, Color, canevas, Dielec
    global r, nAtoms
    Dielec = float(diel.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color = setupall(Atom_Coord, ATOM, Color, 0)
    update_ene()
    canevas.pack()


### set Temperature ###
def set_temp(event):
    global Atom_Coord, ATOM, Color, canevas, Radius
    global r, nAtoms, temp, Temperature, Velocity, Iterations, cstboltz
    Temperature = float(temp.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    canevas.bind("<Button-1>", Go)
    Velocity = InitVel(nAtoms, Temperature, cstboltz, Mass)
    Atom_Coord, ATOM, Color = setupall(Atom_Coord, ATOM, Color, 0)
    update_ene()
    canevas.pack()
    Iterations = 0


### set thermostat coupling time ###
def set_taut(event):
    global Atom_Coord, ATOM, Color, canevas
    global r, nAtoms, taut, tauT
    tauT = float(taut.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color = setupall(Atom_Coord, ATOM, Color, 0)
    update_ene()
    canevas.pack()


def set_spd(event):
    global Atom_Coord, ATOM, Color, canevas, Radius
    global r, nAtoms, spd, speed
    speed = max(1, int(spd.get()))
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color = setupall(Atom_Coord, ATOM, Color, 0)
    canevas.pack()


def set_tstep(event):
    global \
        Atom_Coord, \
        ATOM, \
        Color, \
        canevas, \
        Radius, \
        Iterations, \
        timestep, \
        Velocity, \
        cstboltz
    global r, nAtoms, tstep, Temperature
    timestep = float(tstep.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color = setupall(Atom_Coord, ATOM, Color, 0)
    Velocity = InitVel(nAtoms, Temperature, cstboltz, Mass)
    update_ene()
    canevas.pack()
    Iterations = 0


### update energy ###
def update_ene():
    global \
        Atom_Coord, \
        Radius, \
        BoxDim, \
        Epsilon, \
        Rmin, \
        CutOff, \
        CutOffSquare, \
        Iterations, \
        Ene
    global hwtextene1, hwtextene2
    global Color, sttext0, ptext2, paccept, Dielec, root, canevas, speed
    global Velocity, Mass, cstboltz
    Ene, EneLJ, EneCoul = Calc_Ene2(
        Atom_Coord, Epsilon, Rmin, Dielec, CutOffSquare, BoxDim
    )
    Kin, temperature = calc_temp(Velocity, nAtoms, cstboltz, Mass)
    hwtextene1.destroy()
    hwtextene2.destroy()

    mynewtext = "Etot: %6.1f Ekin: %6.1f Epot: %6.1f" % (Ene + Kin, Kin, Ene)
    hwtextene1 = Label(top2, text=mynewtext)
    hwtextene1.pack(side="left")

    mynewtext = "Elj: %6.1f Ecoul: %6.1f Temp: %6.1f" % (EneLJ, EneCoul, temperature)
    hwtextene2 = Label(top3, text=mynewtext)
    hwtextene2.pack(side="left")


###########################################################################
###########################################################################

################
# MAIN PROGRAM #
################

from tkinter import *
# from Canvas import *

root = Tk()  # root (main) window

canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
canevas.bind("<Button-1>", Go)
root.bind("<Escape>", die)
top = Frame(root)
top.pack(side="top")
top1 = Frame(root)
top1.pack(side="top")
top2 = Frame(root)
top2.pack(side="top")
top3 = Frame(root)
top3.pack(side="top")
low = Frame(root)
low.pack(side="bottom")
low2 = Frame(root)
low2.pack(side="bottom")
low3 = Frame(root)
low3.pack(side="bottom")
mynewtext = "Charged particles MD - velocity-Verlet + Berendsen thermostat"
hwtext = Label(top, text=mynewtext, foreground="red", font="times 18 bold")
hwtext.pack(side="left")
mynewtext = "nAtoms = "
hwtext = Label(low, text=mynewtext)
hwtext.pack(side="left")

r = DoubleVar()
size = DoubleVar()
vdw1 = DoubleVar()
vdw2 = DoubleVar()
frac = DoubleVar()
diel = DoubleVar()
q = DoubleVar()
temp = DoubleVar()
spd = DoubleVar()
tstep = DoubleVar()
taut = DoubleVar()
swp = DoubleVar()

Color = []
ATOM = []  # liste contenant les widgets
Atom_Coord = InitConf(nAtoms, BoxDim, Radius, qat, frac_neg)
Velocity = InitVel(nAtoms, Temperature, cstboltz, Mass)

r.set(nAtoms)
r_entry = Entry(low, width=6, textvariable=r)
r_entry.pack(side="left")
r_entry.bind("<Return>", set_r)

mynewtext2 = "VDW param: Radius = "
hwtext2 = Label(low3, text=mynewtext2)
hwtext2.pack(side="left")
size.set(Radius)
size_entry = Entry(low3, width=6, textvariable=size)
size_entry.pack(side="left")
size_entry.bind("<Return>", set_size)

mynewtext3 = "Epsilon = "
hwtext3 = Label(low3, text=mynewtext3)
hwtext3.pack(side="left")
vdw1.set(Epsilon)
vdw1_entry = Entry(low3, width=6, textvariable=vdw1)
vdw1_entry.pack(side="left")
vdw1_entry.bind("<Return>", set_vdw1)

mynewtext2 = "Coulomb param: frac neg charge = "
hwtext2 = Label(low2, text=mynewtext2)
hwtext2.pack(side="left")
frac.set(frac_neg)
frac_entry = Entry(low2, width=6, textvariable=frac)
frac_entry.pack(side="left")
frac_entry.bind("<Return>", set_frac)

mynewtext3 = "abs charge = "
hwtext3 = Label(low2, text=mynewtext3)
hwtext3.pack(side="left")
q.set(qat)
q_entry = Entry(low2, width=6, textvariable=q)
q_entry.pack(side="left")
q_entry.bind("<Return>", set_q)

mynewtext3 = "Dielec = "
hwtext3 = Label(low2, text=mynewtext3)
hwtext3.pack(side="left")
diel.set(Dielec)
diel_entry = Entry(low2, width=6, textvariable=diel)
diel_entry.pack(side="left")
diel_entry.bind("<Return>", set_diel)

mynewtext3 = "T[K] = "
hwtext3 = Label(low, text=mynewtext3)
hwtext3.pack(side="left")
temp.set(Temperature)
temp_entry = Entry(low, width=6, textvariable=temp)
temp_entry.pack(side="left")
temp_entry.bind("<Return>", set_temp)

mynewtext3 = "tauT = "
hwtext3 = Label(low, text=mynewtext3)
hwtext3.pack(side="left")
taut.set(tauT)
taut_entry = Entry(low, width=6, textvariable=taut)
taut_entry.pack(side="left")
taut_entry.bind("<Return>", set_taut)

mynewtext3 = "Timestep = "
hwtext3 = Label(low, text=mynewtext3)
hwtext3.pack(side="left")
tstep.set(timestep)
tstep_entry = Entry(low, width=6, textvariable=tstep)
tstep_entry.pack(side="left")
tstep_entry.bind("<Return>", set_tstep)

startbutton = Button(
    top3,
    text=" start ",
    command=Go,
    background="yellow",
    foreground="blue",
    relief="flat",
)
startbutton.pack(side="left", fill="x")
stopbutton = Button(
    top3,
    text=" stop ",
    command=stop,
    background="yellow",
    foreground="red",
    relief="flat",
)
stopbutton.pack(side="left")

mynewtext3 = "Delay = "
hwtext3 = Label(low, text=mynewtext3)
hwtext3.pack(side="left")
spd.set(speed)
spd_entry = Entry(low, width=6, textvariable=spd)
spd_entry.pack(side="left")
spd_entry.bind("<Return>", set_spd)

# set up a global var with the nb of iterations
Iterations = 0

# calculate first energy
Ene, EneLJ, EneCoul = Calc_Ene2(Atom_Coord, Epsilon, Rmin, Dielec, CutOffSquare, BoxDim)
Kin, temperature = calc_temp(Velocity, nAtoms, cstboltz, Mass)

mynewtext = "Etot: %6.1f Ekin: %6.1f Epot: %6.1f" % (Ene + Kin, Kin, Ene)
hwtextene1 = Label(top1, text=mynewtext)
hwtextene1.pack(side="left")

mynewtext = "Elj: %6.1f Ecoul: %6.1f Temp: %6.1f" % (EneLJ, EneCoul, temperature)
hwtextene2 = Label(top2, text=mynewtext)
hwtextene2.pack(side="left")

mynewtext = "step: %d Time: %8.3f" % (Iterations, float(Iterations) * timestep)
sttext0 = Label(top3, text=mynewtext)
sttext0.pack(side="left")

Atom_Coord, ATOM, Color = setupall(Atom_Coord, ATOM, Color)

canevas.pack()

print("Click on mouse button within the box to go ahead !!!")
print("<ESC> to quit")

# conserver la fenetre ouverte (inutile en interactif)
root.mainloop()
