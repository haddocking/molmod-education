#!/usr/local/bin/python
#
###########################################################################
###########################################################################
# Simplex (downhill / Nelder-Mead) EM of Lennard-Jones charged particles
# Alexandre Bonvin, Utrecht University
#
# adapted from a script from Patrick Fuchs, Uni. Paris VI
#
# Python 3 port, kept consistent with LJ-ELEC_EM-simplex.ipynb
# (same system + simplex parameters, same energy function with the
#  divide-by-zero floor, seeded initial configuration).
###########################################################################
###########################################################################

##################
# import modules #
##################
import sys
from math import sqrt, exp
from random import random, randint, seed

####################
# define constants #
####################

nDim = 2  # nr of spatial dimensions
nAtoms = 20  # number of atoms
Radius = 25.0  # beware that Radius must be in a good range (according to nAtoms)
# in order to be able to place all atoms
Rmin = 2.24 * Radius  # distance at which the LJ energy is minimal
BoxDim = [500, 500]  # box dimension
Atom_Coord = []  # list of the form : [NMAX][2]
Epsilon = 25.0  # well depth
Dielec = 1.0  # dielectric constant
qat = Radius  # Atom absolute charge
frac_neg = 0.5  # Fraction negative charges
OverlapFr = 0.0  # fraction of overlap allowed
CutOff = 250  # non-bonded cutoff
CutOffSquare = CutOff**2
speed = 20  # canvas update speed
Seed = 100  # random number seed (reproducibility)

# simplex parameters (tuned for the 20-atom system, see the notebook):
FracShrimp1 = 0.5  # contraction of the highest vertex toward the centroid
FracShrimp2 = (
    0.4  # shrink of all vertices toward the lowest (gentle -> avoids early collapse)
)
FracExpend = 2.0  # expansion factor after a successful reflection
deltaE = 0.001  # convergence criterion on the energy improvement
n2conv = 600  # number of sub-deltaE steps (cumulative) to declare convergence
Simplex_step = (
    300.0  # step used to initialize the simplex (large -> more exploration in 40-D)
)
simplex_canvas_factor = (
    1.0  # scale factor for the simplex canvas (only affects its appearance)
)

##################
# some functions #
##################


def dist(A, B):
    return sqrt((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2)


# change sign
def SignR(a, b):
    if b > 0:
        return a
    else:
        return -a


# generate a random rgb color like  #xxxxxx (xx should be an hexadecimal nb)
def random_color():
    tmp = "#"
    for i in range(6):
        rdm = randint(0, 15)
        tmp += hex(rdm)[-1]
    return tmp


# generate a rgb color based on charge like  #xxxxxx (xx should be an hexadecimal nb)
def charge_color(charge, qat):
    if charge == qat:
        tmp = "#FFFFFF"
    else:
        tmp = "#333333"
    return tmp


#########################
# initialize parameters #
#########################


# generates random coordinates
def InitConf(n, dim, radius, qat, frac_neg):
    seed(Seed)
    print("Initializing box, please wait...")
    # generate a list of random positions
    tmp_coord = []
    i = 0
    nneg = 0
    ntrial = 0
    # fix first atom
    x = random() * (dim[0] - radius) + radius
    y = random() * (dim[1] - radius) + radius
    nneg = int(float(n) * frac_neg)
    npos = n - nneg
    charge = -qat
    if npos == n:
        charge = qat
    i += 1
    tmp_coord.append([x, y, charge])
    while i < nneg:
        x = random() * (dim[0] - radius) + radius
        y = random() * (dim[1] - radius) + radius
        # check whether the new particle overlaps an existing one
        OVERLAP = 1
        for j in range(i):
            if dist(tmp_coord[j], [x, y]) < (1 - OverlapFr) * 2 * radius:
                OVERLAP = 0
        if OVERLAP:
            charge = -qat
            tmp_coord.append([x, y, charge])
            i += 1
        ntrial = ntrial + 1
        if ntrial > 100000:
            print("initialisation failed")
            print("==> reduce radius or number of atoms")
            sys.exit()
    while i < n:
        x = random() * (dim[0] - radius) + radius
        y = random() * (dim[1] - radius) + radius
        # check whether the new particle overlaps an existing one
        OVERLAP = 1
        for j in range(i):
            if dist(tmp_coord[j], [x, y]) < (1 - OverlapFr) * 2 * radius:
                OVERLAP = 0
        if OVERLAP:
            charge = qat
            tmp_coord.append([x, y, charge])
            i += 1
        ntrial = ntrial + 1
        if ntrial > 100000:
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


class Simplex:
    def __init__(self, ac, step):
        # one vertex per coordinate direction: the base config displaced by `step`
        tmpss = []
        for pp in range(len(ac)):
            for qq in range(len(ac[pp]) - 1):
                ss = make_simplex_coor(ac, pp, qq, step)
                tmpss.append(ss)
        self.points = tmpss
        self.simplex_energy, self.highest_nr, self.lowest_nr = self.energy()

    def energy(self):
        simplex_count = 0
        simplex_energy = []
        temp_coord = self.points[simplex_count]
        simplex_ene = Calc_Ene(temp_coord, Epsilon, Rmin, Dielec, CutOffSquare, BoxDim)
        highest_ene = simplex_ene
        highest_nr = 0
        lowest_ene = simplex_ene
        lowest_nr = 0
        simplex_count = 1
        simplex_energy.append(simplex_ene)
        for pp in self.points[1:]:
            temp_coord = self.points[simplex_count]
            simplex_ene = Calc_Ene(
                temp_coord, Epsilon, Rmin, Dielec, CutOffSquare, BoxDim
            )
            simplex_energy.append(simplex_ene)
            if simplex_ene > highest_ene:
                highest_ene = simplex_ene
                highest_nr = simplex_count
            if simplex_ene < lowest_ene:
                lowest_ene = simplex_ene
                lowest_nr = simplex_count
            simplex_count = simplex_count + 1
        return simplex_energy, highest_nr, lowest_nr

    def Update_Simplex(self, new, nr):
        count = 0
        news = []
        for i in self.points:
            if count == nr:
                tmplist = new
            else:
                tmplist = i
            news.append(tmplist)
            count += 1
        self.points = news

    def shrimp(self, lowest):
        global nAtoms
        count = 0
        totalnr = nDim * nAtoms
        flat = []
        lowestcoor = self.points[lowest]
        for i in self.points:
            new = []
            if count != lowest:
                for x1 in range(len(i)):
                    for jj in range(len(i[x1][:-1])):
                        tmp = i[x1][jj] + FracShrimp2 * (lowestcoor[x1][jj] - i[x1][jj])
                        new.append(tmp)
            else:
                for x1 in range(len(i)):
                    for jj in range(len(i[x1][:-1])):
                        tmp = lowestcoor[x1][jj]
                        new.append(tmp)
            flat.append(new)
            count += 1
        new = []
        for h2av in flat:
            coord = [
                [h2av[x], h2av[x + 1], self.points[0][x // 2][2]]
                for x in range(0, len(h2av), 2)
            ]
            new.append(coord)
        self.points = new

    def boundary(self, n, nr=0):
        # apply boundary conditions
        if n == 1:
            count = 0
            for i in self.points:
                if count == nr:
                    coord = i
                    for x1 in range(len(coord)):
                        for ii in range(
                            2
                        ):  # ii=0 -> x coordinate ; ii=1 -> y coordinate
                            halfbox = BoxDim[ii] / 2
                            Z = coord[x1][ii]
                            coord[x1][ii] = (
                                Z - SignR(halfbox, Z) - SignR(halfbox, Z - BoxDim[ii])
                            )
                    self.Update_Simplex(coord, count)
                count += 1
        if n == 2:
            count = 0
            for i in self.points:
                coord = i
                for x1 in range(len(coord)):
                    for ii in range(2):  # ii=0 -> x coordinate ; ii=1 -> y coordinate
                        halfbox = BoxDim[ii] / 2
                        Z = coord[x1][ii]
                        coord[x1][ii] = (
                            Z - SignR(halfbox, Z) - SignR(halfbox, Z - BoxDim[ii])
                        )
                self.Update_Simplex(coord, count)
                count += 1

    def updateCanevas(self, can, cg):
        count = -1
        count2 = -1
        for i in self.points:
            count2 = count2 + 1
            tmpcount = -1
            for ii in i:
                count = count + 1
                tmpcount = tmpcount + 1
                can.delete(SIMPLEXATOM[count])
                if count2 == cg:
                    x1 = ii[0] + Radius * 1.2
                    y1 = ii[1] + Radius * 1.2
                    x2 = ii[0] - Radius * 1.2
                    y2 = ii[1] - Radius * 1.2
                else:
                    x1 = ii[0] + Radius
                    y1 = ii[1] + Radius
                    x2 = ii[0] - Radius
                    y2 = ii[1] - Radius
                if count2 == cg:
                    tmp = "#FFF000"
                else:
                    tmp = Color[tmpcount]
                SIMPLEXATOM[count] = can.create_oval(
                    simplex_canvas_factor * x1,
                    simplex_canvas_factor * y1,
                    simplex_canvas_factor * x2,
                    simplex_canvas_factor * y2,
                    fill=tmp,
                )


def Simplex_Reflection(simplex, highest, reflnr=2.0):
    global nAtoms
    count = 0
    totalnr = nDim * nAtoms
    av = [0.0] * totalnr  # centroid of the non-highest vertices
    h2av = [0.0] * totalnr  # reflected point
    h = [0.0] * totalnr  # the highest vertex
    alt2av = [0.0] * totalnr  # contracted point
    # getting average of the non-highest points:
    for i in simplex.points:
        if count != highest:
            tmp = 0
            for xx in i:
                for jj in xx[:-1]:
                    av[tmp] = av[tmp] + jj
                    tmp = tmp + 1
        else:
            tmp = 0
            for xx in i:
                for jj in xx[:-1]:
                    h[tmp] = h[tmp] + jj
                    tmp = tmp + 1
        count += 1
    for i in range(len(av)):
        av[i] = av[i] / float(count - 1)
    for i in range(len(av)):
        # reflection / expansion:
        h2av[i] = h[i] + reflnr * (av[i] - h[i])
        # contraction:
        alt2av[i] = h[i] + FracShrimp1 * (av[i] - h[i])

    coord = [
        [h2av[x], h2av[x + 1], simplex.points[0][x // 2][2]]
        for x in range(0, len(h2av), 2)
    ]
    coord2 = [
        [alt2av[x], alt2av[x + 1], simplex.points[0][x // 2][2]]
        for x in range(0, len(alt2av), 2)
    ]

    return coord, coord2


####################
# calculate energy #
####################


# calculate LJ from the squared distance
def LJ2(distsquare, epsilon, rmin_exp6):
    Z = (1 / distsquare) ** 3 * rmin_exp6
    return epsilon * Z * (Z - 1)


# classical Coulomb from the squared distance
def Coulomb2(r, dielec, qa, qb):
    return qa * qb / (dielec * sqrt(r))


# Calculate energy Evdw + Ecoulomb (uses squared distance)
# version with boundary conditions
def Calc_Ene(coord, epsilon, rmin, dielec, cutoffsquare, boxdim, elec=1):
    Ene = 0.0
    rmin_exp6 = rmin**6
    # doubly nested loop over all particle pairs
    for i in range(len(coord) - 1):
        for j in range(i + 1, len(coord)):
            # calculate the squared atomic distance (nearest image)
            distsquare = 0
            for k in range(2):
                tmp = coord[j][k] - coord[i][k]
                # chooses the nearest image
                halfbox = boxdim[k] / 2
                tmp = (
                    tmp - SignR(halfbox, tmp - halfbox) - SignR(halfbox, tmp + halfbox)
                )
                distsquare += tmp**2
            distsquare = max(
                distsquare, 1e-6
            )  # avoid division by zero (simplex moves can overlap atoms)
            # compute vdw and Coulomb energy
            if distsquare < cutoffsquare:
                qa = coord[i][2]
                qb = coord[j][2]
                Ene += LJ2(distsquare, epsilon, rmin_exp6)
                if elec:
                    Ene += Coulomb2(distsquare, dielec, qa, qb)
    return Ene


def make_simplex_coor(coord, nr, xory, step):
    listnew = []
    for i in range(len(coord)):
        tmplist = []
        for j in range(len(coord[i]) - 1):
            tmp = float(coord[i][j])
            if (i == nr) and (j == xory):
                tmp = tmp + step
            tmplist.append(tmp)
        tmplist.append(coord[i][-1])
        listnew.append(tmplist)
    return listnew


################################
# Simplex energy minimization  #
################################
def Go(*args):
    global \
        Atom_Coord, \
        Radius, \
        FracShrimp1, \
        FracShrimp2, \
        BoxDim, \
        Epsilon, \
        Rmin, \
        CutOffSquare, \
        Iterations, \
        Ene, \
        Accepted, \
        deltaE, \
        Simplex_step
    global \
        Color, \
        sttext0, \
        entext1, \
        ptext2, \
        tmove, \
        paccept, \
        Dielec, \
        Simp, \
        nconv, \
        FracExpend
    ConvCrit1 = deltaE
    if Iterations > 0:
        sttext0.destroy()
        entext1.destroy()
        ptext2.destroy()
        tmove.destroy()

    simplex_energy, highest_nr, lowest_nr = Simp.energy()
    HighEne = simplex_energy[highest_nr]
    LowEne = simplex_energy[lowest_nr]

    mynewtext = "step %d" % Iterations
    sttext0 = Label(top, text=mynewtext)
    sttext0.pack(side="left")
    mynewtext = "E=%.3f" % LowEne
    entext1 = Label(top, text=mynewtext)
    entext1.pack(side="left")
    mynewtext = "Paccept= %.2f" % paccept
    ptext2 = Label(top, text=mynewtext)
    ptext2.pack(side="left")
    mynewtext = "   "
    tmove = Label(top, text=mynewtext)
    tmove.pack(side="left")

    New_Atom_Coord, Alt_Coord = Simplex_Reflection(Simp, highest_nr)
    New_Ene = Calc_Ene(New_Atom_Coord, Epsilon, Rmin, Dielec, CutOffSquare, BoxDim)
    if (HighEne - New_Ene) > ConvCrit1:
        # reflection was successful
        simplexmovetype = "Refl"
        canvasgreen = highest_nr
        Simp.Update_Simplex(New_Atom_Coord, highest_nr)
        Simp.boundary(1, highest_nr)
        # reflection was successful, so try to do an additional expansion
        New_Atom_Coord2, Alt_Coord2 = Simplex_Reflection(Simp, highest_nr, -FracExpend)
        New_Ene2 = Calc_Ene(
            New_Atom_Coord2, Epsilon, Rmin, Dielec, CutOffSquare, BoxDim
        )
        if (New_Ene - New_Ene2) > ConvCrit1:
            simplexmovetype = "Re-E"
            Simp.Update_Simplex(New_Atom_Coord2, highest_nr)
            Simp.boundary(1, highest_nr)
    else:
        New_Ene = Calc_Ene(Alt_Coord, Epsilon, Rmin, Dielec, CutOffSquare, BoxDim)
        if (HighEne - New_Ene) > ConvCrit1:
            # contraction of the highest vertex succeeded
            simplexmovetype = "High"
            canvasgreen = highest_nr
            Simp.Update_Simplex(Alt_Coord, highest_nr)
            Simp.boundary(1, highest_nr)
        else:
            # nothing helped -> shrink the whole simplex toward the lowest vertex
            Simp.shrimp(lowest_nr)
            canvasgreen = lowest_nr
            Simp.boundary(2)
            simplexmovetype = "Lowe"

    simplex_energy, highest_nr, lowest_nr = Simp.energy()
    # get the coordinates corresponding to the current lowest energy configuration
    # in the simplex, to be plotted on the screen
    Atom_Coord = Simp.points[lowest_nr]
    NewLowEne = simplex_energy[lowest_nr]
    Iterations = Iterations + 1
    if Iterations > 0:
        sttext0.destroy()
        entext1.destroy()
        ptext2.destroy()
        tmove.destroy()
    mynewtext = "step %d" % Iterations
    sttext0 = Label(top, text=mynewtext)
    sttext0.pack(side="left")
    mynewtext = "E= %.3f %s" % (NewLowEne, simplexmovetype)
    entext1 = Label(top, text=mynewtext)
    entext1.pack(side="left")
    atoms_to_do = range(nAtoms)

    for ii in atoms_to_do:
        canevas.delete(ATOM[ii])
        # draw new
        x1 = Atom_Coord[ii][0] + Radius
        y1 = Atom_Coord[ii][1] + Radius
        x2 = Atom_Coord[ii][0] - Radius
        y2 = Atom_Coord[ii][1] - Radius
        ATOM[ii] = canevas.create_oval(x1, y1, x2, y2, fill=Color[ii])
    Ene_diff = LowEne - NewLowEne
    Simp.updateCanevas(simplexcanevas, canvasgreen)
    if abs(Ene_diff) < ConvCrit1:
        nconv = nconv + 1
        if nconv > n2conv:
            print("convergence")
            stop()
        else:
            canevas.after(speed, Go)
    else:
        canevas.after(speed, Go)


########################
# some other functions #
########################
def die(event=0):
    import sys

    sys.exit()


###########################################################################
###########################################################################


# some functions added for extra buttons etc.
def stop():
    global Atom_Coord, ATOM, Color, canevas, Simp, simplexcanevas
    global r, nAtoms, size, Radius
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color, 0)
    canevas.pack()
    nconv = 0


def setupall(Atom_Coord, ATOM, Color, Repack=1):
    global Iterations, Simp, nconv
    nconv = 0
    ATOM = []  # list holding the widgets
    SIMPLEXATOM = []  # list for the simplex

    if Repack == 1:
        Atom_Coord = InitConf(nAtoms, BoxDim, Radius, qat, frac_neg)
        Color = []
        for i in range(nAtoms):
            Color.append(charge_color(Atom_Coord[i][2], qat))
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
    if Repack == 1:
        Simp = Simplex(Atom_Coord, Simplex_step)
    for i in Simp.points:
        count = -1
        for p in i:
            count = count + 1
            x1 = p[0] + Radius
            y1 = p[1] + Radius
            x2 = p[0] - Radius
            y2 = p[1] - Radius
            tmp = Color[count]
            SIMPLEXATOM.append(
                simplexcanevas.create_oval(
                    simplex_canvas_factor * x1,
                    simplex_canvas_factor * y1,
                    simplex_canvas_factor * x2,
                    simplex_canvas_factor * y2,
                    fill=tmp,
                )
            )

    return Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM


def set_r(event):
    global Atom_Coord, ATOM, Color, canevas, Simp, simplexcanevas, SIMPLEXATOM
    global r, nAtoms
    nAtoms = int(r.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    canevas.bind("<Button-1>", Go)
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color)
    canevas.pack()


def set_size(event):
    global Atom_Coord, ATOM, Color, canevas, Radius, Simp, simplexcanevas, SIMPLEXATOM
    global r, nAtoms, size
    Radius = int(size.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color)
    canevas.pack()


def set_vdw1(event):
    global Atom_Coord, ATOM, Color, canevas, Radius, Simp, simplexcanevas, SIMPLEXATOM
    global r, nAtoms, vdw1, Epsilon
    Epsilon = int(vdw1.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color, 0)
    canevas.pack()


def set_vdw2(event):
    global Atom_Coord, ATOM, Color, canevas, Radius, Simp, simplexcanevas, SIMPLEXATOM
    global r, nAtoms, vdw2, Rmin
    Rmin = int(vdw2.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color, 0)
    canevas.pack()


def set_frac(event):
    global Atom_Coord, ATOM, Color, canevas, frac_neg, Simp, simplexcanevas, SIMPLEXATOM
    global r, nAtoms, temp
    frac_neg = float(frac.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color, 2)
    canevas.pack()


def set_q(event):
    global Atom_Coord, ATOM, Color, canevas, qat, Simp, simplexcanevas, SIMPLEXATOM
    global r, nAtoms, temp
    qat = float(q.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color, 2)
    canevas.pack()


def set_diel(event):
    global Atom_Coord, ATOM, Color, canevas, Dielec, Simp, simplexcanevas, SIMPLEXATOM
    global r, nAtoms, temp
    Dielec = float(diel.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color)
    canevas.pack()


def set_dmv(event):
    global Atom_Coord, ATOM, Color, canevas, Radius, Simp, simplexcanevas, SIMPLEXATOM
    global r, nAtoms, dmv, FracShrimp1, FracShrimp2, Simplex_step
    FracShrimp1 = float(dmv.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color, 0)
    canevas.pack()


def set_dmv2(event):
    global Atom_Coord, ATOM, Color, canevas, Radius, Simp, simplexcanevas, SIMPLEXATOM
    global r, nAtoms, dmv, dmv2, FracShrimp1, FracShrimp2, Simplex_step
    FracShrimp2 = float(dmv2.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color, 0)
    canevas.pack()


def set_ss(event):
    global Atom_Coord, ATOM, Color, canevas, Radius, Simp, simplexcanevas, SIMPLEXATOM
    global r, nAtoms, ss, FracShrimp1, FracShrimp2, Simplex_step
    Simplex_step = float(ss.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color)
    canevas.pack()


def set_swp(event):
    global Atom_Coord, ATOM, Color, canevas, deltaE, Simp, simplexcanevas, SIMPLEXATOM
    global r, nAtoms, dmv, dmv2, FracShrimp1, FracShrimp2, Simplex_step
    deltaE = float(swp.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color, 0)
    canevas.pack()


def set_n2c(event):
    global Atom_Coord, ATOM, Color, canevas, n2conv, Simp, simplexcanevas, SIMPLEXATOM
    global r, nAtoms, dmv, dmv2, FracShrimp1, FracShrimp2, Simplex_step
    n2conv = float(n2c.get())
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color, 0)
    canevas.pack()


def set_spd(event):
    global Atom_Coord, ATOM, Color, canevas, Radius, Simp, simplexcanevas, SIMPLEXATOM
    global r, nAtoms, spd, speed
    speed = max(1, int(spd.get()))
    canevas.destroy()
    canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
    simplexcanevas.destroy()
    simplexcanevas = Canvas(
        root,
        width=BoxDim[0] * simplex_canvas_factor,
        height=BoxDim[1] * simplex_canvas_factor,
        bg="#cdddff",
    )
    simplexcanevas.pack(side="right")
    canevas.bind("<Button-1>", Go)
    Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color, 0)
    canevas.pack()


###########################################################################
###########################################################################

################
# MAIN PROGRAM #
################

from tkinter import *

root = Tk()  # root (main) window

root.bind("<Escape>", die)
top2 = Frame(root)
top2.pack(side="top")
top = Frame(root)
top.pack(side="top")
low = Frame(root)
low.pack(side="top")
low2 = Frame(root)
low2.pack(side="top")
low3 = Frame(root)
low3.pack(side="top")
low4 = Frame(root)
low4.pack(side="top")
mynewtext = "Charged particles Simplex"
hwtext = Label(top2, text=mynewtext, foreground="red", font="times 18 bold")
hwtext.pack(side="left")
mynewtext = "nAtoms = "
hwtext = Label(low, text=mynewtext)
hwtext.pack(side="left")

canevas = Canvas(root, width=BoxDim[0], height=BoxDim[1], bg="#ccddff")
simplexcanevas = Canvas(
    root,
    width=BoxDim[0] * simplex_canvas_factor,
    height=BoxDim[1] * simplex_canvas_factor,
    bg="#cdddff",
)
simplexcanevas.pack(side="right")
canevas.bind("<Button-1>", Go)

r = DoubleVar()
size = DoubleVar()
vdw1 = DoubleVar()
vdw2 = DoubleVar()
frac = DoubleVar()
diel = DoubleVar()
q = DoubleVar()
spd = DoubleVar()
dmv = DoubleVar()
dmv2 = DoubleVar()
ss = DoubleVar()
swp = DoubleVar()
n2c = DoubleVar()

Color = []
ATOM = []  # list holding the widgets
Atom_Coord = InitConf(nAtoms, BoxDim, Radius, qat, frac_neg)

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

mynewtext3 = "Rmin = "
hwtext3 = Label(low3, text=mynewtext3)
hwtext3.pack(side="left")
vdw2.set(Rmin)
vdw2_entry = Entry(low3, width=6, textvariable=vdw2)
vdw2_entry.pack(side="left")
vdw2_entry.bind("<Return>", set_vdw2)


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

mynewtext3 = "Shrimp = "
hwtext3 = Label(low4, text=mynewtext3)
hwtext3.pack(side="left")
dmv.set(FracShrimp1)
dmv_entry = Entry(low4, width=4, textvariable=dmv)
dmv_entry.pack(side="left")
dmv_entry.bind("<Return>", set_dmv)
dmv2.set(FracShrimp2)
dmv2_entry = Entry(low4, width=4, textvariable=dmv2)
dmv2_entry.pack(side="left")
dmv2_entry.bind("<Return>", set_dmv2)

mynewtext3 = "StepSize = "
hwtext3 = Label(low4, text=mynewtext3)
hwtext3.pack(side="left")
ss.set(Simplex_step)
ss_entry = Entry(low4, width=6, textvariable=ss)
ss_entry.pack(side="left")
ss_entry.bind("<Return>", set_ss)

mynewtext3 = "Conv = "
hwtext3 = Label(low4, text=mynewtext3)
hwtext3.pack(side="left")
swp.set(deltaE)
swp_entry = Entry(low4, width=6, textvariable=swp)
swp_entry.pack(side="left")
swp_entry.bind("<Return>", set_swp)

mynewtext3 = "N2conv = "
hwtext3 = Label(low4, text=mynewtext3)
hwtext3.pack(side="left")
n2c.set(n2conv)
n2c_entry = Entry(low4, width=6, textvariable=n2c)
n2c_entry.pack(side="left")
n2c_entry.bind("<Return>", set_n2c)

startbutton = Button(
    top,
    text=" start ",
    command=Go,
    background="yellow",
    foreground="blue",
    relief="flat",
)
startbutton.pack(side="left", fill="x")
stopbutton = Button(
    top,
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

Atom_Coord, ATOM, Color, Simp, SIMPLEXATOM = setupall(Atom_Coord, ATOM, Color)

# set up a global var with the nb of iterations and nb of accepted conf
Iterations = 0
Accepted = 0
paccept = 1
# set up Ene as a global variable
Ene = 0.0
counter = 0
nconv = 0
canevas.pack()
simplexcanevas.pack(side="right")

print("Click on mouse button within the box to go ahead !!!")
print("<ESC> to quit")

# keep the window open
root.mainloop()
