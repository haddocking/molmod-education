"""
PyMolTris 2.0 - Tetris inside PyMOL.  (Python 3 / numpy / PyMOL 3.x port of pymoltris2.py)

Run it from the command line:

    pymol pymoltris3.py

The game is driven entirely through PyMOL - keyboard keys plus commands typed in the
PyMOL command line (no separate window, so it works under the default Qt GUI).

  keys : Left/Right = move    Up/Down = rotate    PgDn = drop
         (PyMOL only sends keys to the 3-D view once you CLICK on it - do that first)

  commands (type in the PyMOL command line):
    start                start / pause the game
    pause                pause                     drop        drop the block
    savegame [file]      loadgame [file]           quitgame    quit PyMOL
    periodic             toggle periodic field     shownext    toggle next-block preview
    moreblocks           toggle extra blocks
    blockcolor <name>    blockstyle <n | Default | Small | ...>
    jokers: atwa  collapse  turnaround  destroyc  destroyr
    tetris_help          print this help again

Original (Python 2) author: MK <mkrzemin@nmr.chem.uu.nl>.  Ported to Python 3 / numpy /
PyMOL 3.x, and the Tkinter control panel replaced with native PyMOL keys/commands
(a Tk window cannot coexist with PyMOL's Qt GUI on macOS - it hard-crashes the app).
"""
from random import randint
from threading import Thread, Event
import pickle as c, copy
import socket
import os, sys, math, time, random
import numpy as np

Name    = 'PyMolTris'
Author  = 'MK'
Contact = 'mkrzemin@nmr.chem.uu.nl'
Game    = 'TETRIS'
Version = 2.0
PyMolWebsite = 'http://pymol.sourceforge.net/'


"""
TODO = 	Rotation in the periodic field

Jokers	1. Allthewayaround : the field is returned and every ball is falling down
	2. Shake : All ball are falling down
	3. Turn around : ONLY the field rotates 180 degres

"""

class Control:

	police={'button':("Times", "12", "bold"), 'cbutton':("Times", "12"), 'scale':("Times", "10", "bold"), 'menu':("Arial", "10", "bold")}

	def __init__(self):

		self._PAUSE 	= 1
		self._ANIM  	= 0
		self._AS 	= 0
		self.show_next 	= 1
		self.scani      = 1
		self._destroyc  = 0
		self._destroyr  = 0


		# --- controls via PyMOL (no Tkinter window; works under the default Qt GUI) ---
		self.per = 0		# Periodic field
		self.shn = 0		# Hide next block
		self.mbl = 0		# More blocks
		self._start_action = self.PoP	# what the 'start' command does

		cmd.set_key("left",  self.left)		# move the block left
		cmd.set_key("right", self.right)	# move the block right
		# rotation on several keys - on some macOS setups the Up/Down arrows never reach
		# the viewer, so PgUp/Home/CTRL-R (and the 'rotate' command) rotate as well
		for _k in ("up", "pgup", "home", "CTRL-R"):
			try:    cmd.set_key(_k, self.rotatel)	# rotate one way
			except Exception: pass
		for _k in ("down", "end", "CTRL-E"):
			try:    cmd.set_key(_k, self.rotater)	# rotate the other way
			except Exception: pass
		for _k in ("pgdn", "insert", "CTRL-D"):
			try:    cmd.set_key(_k, self.fall)	# drop the block
			except Exception: pass

		# PyMOL calls extended commands with a _self= keyword, so every callback swallows **k
		cmd.extend("start",      lambda *a, **k: self._start_action())
		cmd.extend("pause",      lambda *a, **k: self.pause())
		cmd.extend("drop",       lambda *a, **k: self.fall())
		cmd.extend("rotate",     lambda *a, **k: self.rotatel())
		cmd.extend("rotback",    lambda *a, **k: self.rotater())
		cmd.extend("savegame",   lambda fname="pymoltris.pmt", **k: self.savegame(fname))
		cmd.extend("loadgame",   lambda fname="pymoltris.pmt", **k: self.loadgame(fname))
		cmd.extend("quitgame",   lambda *a, **k: self.quit())
		cmd.extend("periodic",   lambda *a, **k: self.ponoff(1))
		cmd.extend("shownext",   lambda *a, **k: self.togglenext())
		cmd.extend("moreblocks", lambda *a, **k: self.togglemore())
		cmd.extend("blockcolor", lambda color, **k: self.chgcol(color))
		cmd.extend("blockstyle", lambda style, **k: self.quality(style))
		for _j in ("atwa", "collapse", "turnaround", "destroyc", "destroyr"):
			cmd.extend(_j, (lambda j: (lambda *a, **k: getattr(self, j)()))(_j))
		cmd.extend("tetris_help", lambda *a, **k: self.helptext())

		self.helptext()

	def _status(self, txt):
		print("[PyMolTris] " + str(txt))

	def helptext(self, *a):
		print("")
		print("===================== PyMolTris =====================")
		print(" keys  : Left/Right = move   Up or PgUp = rotate   PgDn = drop")
		print("         >>> first CLICK on the 3-D display so the keys reach the game <<<")
		print("         (if Up/Down do nothing on macOS, use PgUp to rotate, or type 'rotate')")
		print(" start        start / pause the game")
		print(" pause        pause the game            drop        drop the block")
		print(" rotate       rotate the block          rotback     rotate the other way")
		print(" savegame [file]   loadgame [file]      quitgame    quit PyMOL")
		print(" periodic     toggle periodic field     shownext    toggle next preview")
		print(" moreblocks   toggle extra blocks")
		print(" blockcolor <name>     blockstyle <n|Default|Small|...>")
		print(" jokers: atwa  collapse  turnaround  destroyc  destroyr")
		print("=====================================================")
		print("")

	def togglenext(self, *a):
		self.shn = 0 if self.shn else 1
		self.shownext()
		self._status("Next block " + ("hidden" if self.shn else "shown"))

	def togglemore(self, *a):
		self.mbl = 0 if self.mbl else 1
		self._status("Extra blocks " + ("ON" if self.mbl else "OFF"))


	def quality(self, t):
		try:
			cmd.set('sphere_quality', int(t))
		except:
			if (t == "Small"):
				cmd.set("sphere_scale", 1.95)
			elif (t == "Default"):
				cmd.set("sphere_scale", 2.9)
				cmd.set('sphere_transparency', 0.00)
				cmd.set('sphere_quality', 2)
			else:
				cmd.set('sphere_transparency', 0.75)
	def startgame(self, event=None):
		if not self._AS: Mick.start()
		self._AS = 1
		self._PAUSE = 0
		self._status("running")

	def chgcol(self, color):
		self.ANIM = 1
		global _COLF

		init = cmd.get_color_tuple(_COLF)
		end  = cmd.get_color_tuple(color)


		nb_trans = 20
		trans = [(float(end[i])-init[i])/nb_trans for i in range(len(init))]

		for i in range(nb_trans):
			cmd.set_color("new_color", init)
			cmd.color("new_color", "_Field")
			init = [init[i]+trans[i] for i in range(3)]
			time.sleep(0.1)
		_COLF = color

		self._ANIM=0

	def restart(self):
		block.gameover._stag.set()
		cmd.delete("GameOver")
		globals()['grid'] = np.zeros(area, dtype=int)
		init_setup()
		self._status("running")
		control._start_action = self.PoP
		self._PAUSE = 0
			
	def pause(self, event=None):
		if self._PAUSE:
			self._PAUSE = 0
			self._status("running")
		else:
			self._PAUSE = 1
			self._status("paused")

	def ponoff(self, toggle=0):
		if toggle:
			self.per = 0 if self.per else 1
		if not self.per:
			posmin = min([a[1] for a in block.cur_pos])
			posmax = max([a[1] for a in block.cur_pos])
			if (posmax-posmin > 3):
				self.per = 1
				self._status("Cannot turn the periodic field off while the block is split.")
				return
		self._status("Periodic field " + ("ON" if self.per else "OFF"))

	def shownext(self):
		if self.shn:
			cmd.disable("_Next")
		else:
			cmd.enable("_Next")

	def collapse(self):
		if self._ANIM or not self._AS: return
		global grid
		self._ANIM = 1
		self.scani = 0
		self.joker.entryconfig(1, state="disabled")

		if (min([a[0] for a in block.cur_pos])<4):
			m = min([a[0] for a in block.cur_pos])
			tmp_cur_pos = [[a[0]+4-m, a[1]] for a in block.cur_pos]
			for i in tmp_cur_pos:
				if (grid[tuple(i)]):
					self._status("You can not use this option now!")
					return
				else:
					for i in range(len(tmp_cur_pos)):
						block.cur_pos[i] = tmp_cur_pos[i]
		for i in block.cur_pos: grid[tuple(i)] = 1

		block.build_field()
		cmd.delete("_block")
		cmd.origin("_Fake")
		for i in range(90):
			cmd.rotate("x", 2., "_Field", 1, 1, "_Field")
			time.sleep(0.01)

		ntb  = int(np.sum(grid != 0))
		grid = np.zeros(area, dtype=int)	# Reinitialization of the grid

		tmps = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
		for i in range (1, ntb+1, 1):
			pos = cmd.get_atom_coords("_Field and resi "+str(i))
			index_x  = int(round((pos[0] + 45) / _SIZE))
			init_pos = 200 - pos[1]
			fina_pos = -5. + _SIZE * (np.sum(grid, 0)[index_x] + 1)
			diff = (fina_pos - init_pos) / 10.

			if (not fina_pos+5 < _eps):
				for j in tmps:
					cmd.translate([0., -diff, 0.], "resi "+str(i))
					time.sleep(0.01/j)
			grid[area[0]-np.sum(grid, 0)[index_x]-1, index_x] = 1

		block.build_field()
		block.newblock()
		cmd.enable("_block")
		block.checkline()
		cmd.enable("_block")

		self.scani = 1
		self._ANIM = 0

	def turnaround(self):
		if self._ANIM or not self._AS: return
		self.joker.entryconfig(2, state="disabled")
		pass

	def atwa(self):
		if self._ANIM or not self._AS: return
		global grid
		self.joker.entryconfig(0, state="disabled")
		if grid.any():
			self._ANIM = 1

			cmd.delete("_block")

			for k in range(4, area[0], 1):
				if grid[k].any():
					beg = k
					break

			if (area[0]-k)%2:	# odd
				for k in range(area[1]):
					if grid[(area[0]+beg-1)//2,k]:
						orig = str(grid[(area[0]+beg-1)//2,k])
						break
			else:			# even
				orig = ""
				for k in range(area[1]):
					if grid[(area[0]+beg-2)//2,k]:
						orig += str(grid[(area[0]+beg-2)//2,k])+"-"
						break
				for k in range(area[1]):
					if grid[(area[0]+beg)//2,k]:
						orig += str(grid[(area[0]+beg)//2,k])
						break

			cmd.origin("_Field and resi "+orig)

			for i in range(90):
				cmd.rotate("x", 2., "_Field", 1, 1, "_Field")
				time.sleep(0.01)

			tmp_grid = grid.copy()

			for i in range(area[1]):
				for j in range(beg, area[0], 1):
					grid[j,i] = tmp_grid[area[0]-j-1+beg,i]

			block.build_field()

			block.newblock()
			self._ANIM = 0

	def destroyc(self):
		if self._ANIM or not self._AS: return
		global grid
		self._ANIM     = 1
		self._destroyc = 1
		self.joker.entryconfig(3, state="disabled")

		# First, we check the field in
		self.allow = []
		for i in range(area[1]-1, -1, -1):
			if np.take(grid,(i,), 1).any():
				self.maxix = i
				self.allow.append(i)
		if self.allow:
			self.allow.reverse()
		else:
			self._ANIM     = 0
			self._destroyc = 0
			return

		# Then, the text
		
		axes = [[15.0,0.0,0.0],[0.0,15.0,0.0],[0.0,0.0,15.0]]
		des = []
		des_txt= []
		cyl_text(des_txt, plain, [-59.0, 220.0, 25.0], "Select a column", 1.2, color=[.2, .2, .2], axes=axes)
		des.extend(des_txt)
		cmd.load_cgo(des, '_destroy', 1)

		# Finally, the animation

		for maxiy in range(4, area[0]):
			if grid[maxiy].any():
				startx = -45.+self.maxix*_SIZE
				starty = 235.-(maxiy-2)*_SIZE

				arrow = [COLOR,1.0,0.1,0.1,BEGIN, TRIANGLES,\
					 VERTEX, startx-5., starty, 15.,\
					 VERTEX, startx, starty-_SIZE, 15.,\
					 VERTEX, startx+5, starty, 15.,\
					 END, BEGIN, LINE_STRIP,\
					 VERTEX, startx, starty, 15.,\
					 VERTEX, startx, starty+2*_SIZE, 15.,\
					 END]

				cmd.load_cgo(arrow, '_Arrow', 1)
				break



	def destroyr(self):
		if self._ANIM or not self._AS: return
		global grid
		self._ANIM     = 1
		self._destroyr = 1
		self.joker.entryconfig(4, state="disabled")

		# First, we check the field in
		self.allow = []
		for i in range(4, area[0]):
			if grid[i].any():
				self.allow.append(i)

		if not self.allow:
			self._ANIM     = 0
			self._destroyr = 0
			return


		# Then, the text
		
		axes = [[15.0,0.0,0.0],[0.0,15.0,0.0],[0.0,0.0,15.0]]
		des = []
		des_txt= []
		cyl_text(des_txt, plain, [-49.0, 220.0, 25.0], "Select a raw", 1.2, color=[.2, .2, .2], axes=axes)
		des.extend(des_txt)
		cmd.load_cgo(des, '_destroy', 1)

		# Finally, the animation
		

		for self.maxiy in range(4, area[0]):
			if grid[self.maxiy].any():
				self.allowr = list(range(self.maxiy, area[0]))
				startx = -45. - 2 * _SIZE
				starty = 235.- self.maxiy * _SIZE

				arrow = [COLOR,1.0,0.1,0.1,BEGIN, TRIANGLES,\
					 VERTEX, startx, starty-5., 15.,\
					 VERTEX, startx+_SIZE, starty, 15.,\
					 VERTEX, startx, starty+5., 15.,\
					 END, BEGIN, LINE_STRIP,\
					 VERTEX, startx, starty, 15.,\
					 VERTEX, startx-2*_SIZE, starty, 15.,\
					 END]

				cmd.load_cgo(arrow, '_Arrow', 1)
				break
		else:
			self._ANIM = 0


	def checkit(self, np):
		for i in np:
			try:
				if (grid[i[0]][i[1]] or i[1] < 0): return 0
			except:
				return 0
		return 1

	def savegame(self, fname="pymoltris.pmt"):
		if not self._PAUSE: self.pause()
		try:
			data = {"field":grid, "block":block.cur_pos, "next":block.next, "score":score.score, "color":block.rand_bl}
			with open(fname, "wb") as f: c.dump(data, f)
			self._status("Game saved to " + fname + " - type 'start' to go on playing.")
		except Exception as e:
			self._status("Could not save the game: " + str(e))
				
	def loadgame(self, fname="pymoltris.pmt"):
		if not self._PAUSE: self.pause()
		try:
			with open(fname, "rb") as f: data = c.load(f)
			globals()['grid'] = data['field']
			block.cur_pos = data['block']
			block.next = data['next']
			score.score = data['score']
			block.rand_bl = data["color"]
			block.build_field()
			block.newblock(block.cur_pos)
			score.updatescore(score.score)
			self._status("Game loaded from " + fname + " - type 'start' to go on playing.")
		except Exception as e:
			self._status("Could not load the game: " + str(e))

	def PoP(self):
		if not (self._PAUSE or self._ANIM):
			self.pause(None)
		else:
			self._PAUSE = 0
			self.startgame(None)
	

	def left(self, event=None):

		if not (self._ANIM):
			if self.per:
				for i in block.cur_pos:
					if (i[1] < 1):
						if (grid[i[0],area[1]-1]):
							return
					elif grid[i[0],i[1]-1]:
						return
			else:
				for i in block.cur_pos:
					if (i[1] < 1):
						return
					elif (grid[i[0],i[1]-1]):
						return

			if (self.per and min([a[1] for a in block.cur_pos]) < 1):
				for i in range(len(block.blocktype[block.rand_bl])):
					if (block.cur_pos[i][1] < 1):
						cmd.translate([(area[1]-1)*_SIZE, 0, 0], "_block and resi "+str(i+300))
						block.cur_pos[i][1] = area[1]-1
					else:
						cmd.translate([-_SIZE, 0, 0], "_block and resi "+str(i+300))
						block.cur_pos[i][1] -= 1
				return
			else:
				cmd.translate([-_SIZE, 0, 0], "_block")
				for i in block.cur_pos: i[1] -= 1
				return
		elif(self._destroyc):
			try:
				if (self.allow.index(self.maxix) != 0):
					cmd.translate([_SIZE*(self.allow[self.allow.index(self.maxix)-1]-self.maxix), 0, 0], object="_Arrow")
					self.maxix = self.allow[self.allow.index(self.maxix)-1]
				else:
					cmd.translate([_SIZE*(self.allow[-1]-self.maxix), 0, 0], object="_Arrow")
					self.maxix = self.allow[-1]
			except:
				pass

	def right(self, event=None):
		if not (self._PAUSE or self._ANIM):
			if self.per:
				for i in block.cur_pos:
					if (i[1] > area[1]-2):
						if (grid[i[0],0]):
							return
					elif (grid[i[0],i[1]+1]):
						return
			else:
				for i in block.cur_pos:
					if (i[1] > area[1]-2):
						return
					elif (grid[i[0],i[1]+1]):
						return

			if (self.per and max([a[1] for a in block.cur_pos]) > area[1]-2):
				for i in range(len(block.blocktype[block.rand_bl])):
					if (block.cur_pos[i][1] > area[1]-2):
						cmd.translate([-(area[1]-1)*_SIZE, 0, 0], "_block and resi "+str(i+300))
						block.cur_pos[i][1] = 0
					else:
						cmd.translate([_SIZE, 0, 0], "_block and resi "+str(i+300))
						block.cur_pos[i][1] += 1
			else:
				cmd.translate([_SIZE, 0, 0], "_block")
				for i in block.cur_pos: i[1] += 1
		elif (self._destroyc):
			try:				
				if (self.allow.index(self.maxix) != len(self.allow)-1):
					cmd.translate([_SIZE*(self.allow[self.allow.index(self.maxix)+1]-self.maxix), 0, 0], object="_Arrow")
					self.maxix = self.allow[self.allow.index(self.maxix)+1]
				else:
					cmd.translate([_SIZE*(self.allow[0]-self.maxix), 0, 0], object="_Arrow")
					self.maxix = self.allow[0]
			except:
				pass

	def rotatel(self, event=None):	# Up
		if not (self._PAUSE or self._ANIM):
			if not self.per:		#(max([a[1] for a in block.cur_pos]) - min([a[1] for a in block.cur_pos]) > 3)):
				if block.rand_bl not in block.norot:
					tmp_pos = [[block.cur_pos[0][0],block.cur_pos[0][1]],]
			
					for i in block.cur_pos[1:]:
						tmp_pos.append([tmp_pos[0][0]+tmp_pos[0][1]-i[1], tmp_pos[0][1]-tmp_pos[0][0]+i[0]])
						if (tmp_pos[-1][1]>9 or tmp_pos[-1][1]<0 or tmp_pos[-1][0]<0):		# Testing the position about the edges
							return
						if (grid[tmp_pos[-1][0], tmp_pos[-1][1]]):				# Testing the position about the other blocks
							return

					cmd.origin("_block and resi 300")
					cmd.rotate("z", 90, "_block")			
					block.cur_pos = []
					block.upload_cur_pos(tmp_pos)
					#for i in tmp_pos: block.cur_pos.append([i[0], i[1]])
			else:
				# Rotation of the splitted block
				if block.rand_bl not in block.norot:
					tmp_pos = [[block.cur_pos[0][0],block.cur_pos[0][1]],]

					if (max([a[1] for a in block.cur_pos]) - min([a[1] for a in block.cur_pos]) > 3):
						# Re-assembling the block
						for i in range(1, len(block.cur_pos)):
							if (block.cur_pos[0][1]-block.cur_pos[i][1] > 3):
								block.cur_pos[i][1] += area[1]
								cmd.translate([area[1]*_SIZE,0,0], "resi "+str(300+i))
							elif (block.cur_pos[0][1]-block.cur_pos[i][1] < -3):
								block.cur_pos[i][1] -= area[1]
								cmd.translate([-area[1]*_SIZE,0,0], "resi "+str(300+i))

					for i in block.cur_pos[1:]:
						x = tmp_pos[0][1] - tmp_pos[0][0]+i[0]
						if (x > area[1]-1):
							x -= area[1]
						elif (x < 0):
							x += area[1]

						tmp_pos.append([tmp_pos[0][0]+tmp_pos[0][1]-i[1], x])
						if (grid[tmp_pos[-1][0], tmp_pos[-1][1]]): return			# Testing the position about the other blocks

					
					for i in range(1, len(block.cur_pos)):
						xstep = (tmp_pos[i][1] - block.cur_pos[i][1]) * _SIZE
						ystep = (tmp_pos[i][0] - block.cur_pos[i][0]) * _SIZE
						cmd.translate([xstep, -ystep, 0], "resi "+str(300+i))
					block.upload_cur_pos(tmp_pos)
		elif (self._destroyr):
			try:
				if (self.allow.index(self.maxiy) != 0):
					cmd.translate([0, _SIZE*(self.maxiy-self.allow[self.allow.index(self.maxiy)-1]), 0], object="_Arrow")
					self.maxiy = self.allow[self.allow.index(self.maxiy)-1]
				else:
					cmd.translate([0, _SIZE*(self.maxiy-self.allow[len(self.allow)-1]), 0], object="_Arrow")
					self.maxiy = self.allow[-1]
			except:
				pass

	def rotater(self, event=None):	# Down
		if not (self._PAUSE or self._ANIM):
			if not self.per:
				if block.rand_bl not in block.norot:
					tmp_pos = [[block.cur_pos[0][0],block.cur_pos[0][1]],]
			
					for i in block.cur_pos[1:]:
						tmp_pos.append([tmp_pos[0][0]-tmp_pos[0][1]+i[1], tmp_pos[0][0]+tmp_pos[0][1]-i[0]])
						if (tmp_pos[-1][1]>9 or tmp_pos[-1][1]<0 or tmp_pos[-1][0]<0):		# Testing the position about the edges
							return
						if (grid[tmp_pos[-1][0], tmp_pos[-1][1]]):				#Testing the position about the other blocks
							return

					cmd.origin("_block and resi 300")
					cmd.rotate("z", -90, "_block")
					block.cur_pos = []
					block.upload_cur_pos(tmp_pos)
					#for i in tmp_pos: block.cur_pos.append([i[0], i[1]])
			else:
				# Rotation of the splitted block
				if block.rand_bl not in block.norot:
					tmp_pos = [[block.cur_pos[0][0],block.cur_pos[0][1]],]

					if (max([a[1] for a in block.cur_pos]) - min([a[1] for a in block.cur_pos]) > 3):
						# Re-assembling the block
						for i in range(1, len(block.cur_pos)):
							if (block.cur_pos[0][1]-block.cur_pos[i][1] > 3):
								block.cur_pos[i][1] += area[1]
								cmd.translate([area[1]*_SIZE,0,0], "resi "+str(300+i))
							elif (block.cur_pos[0][1]-block.cur_pos[i][1] < -3):
								block.cur_pos[i][1] -= area[1]
								cmd.translate([-area[1]*_SIZE,0,0], "resi "+str(300+i))

					for i in block.cur_pos[1:]:
						x = tmp_pos[0][0]+tmp_pos[0][1]-i[0]
						if (x > area[1]-1):
							x -= area[1]
						elif (x < 0):
							x += area[1]

						tmp_pos.append([tmp_pos[0][0]-tmp_pos[0][1]+i[1], x])
						if (grid[tmp_pos[-1][0], tmp_pos[-1][1]]): return			# Testing the position about the other blocks

					
					for i in range(1, len(block.cur_pos)):
						xstep = (tmp_pos[i][1] - block.cur_pos[i][1]) * _SIZE
						ystep = (tmp_pos[i][0] - block.cur_pos[i][0]) * _SIZE
						cmd.translate([xstep, -ystep, 0], "resi "+str(300+i))
					block.upload_cur_pos(tmp_pos)
		elif (self._destroyr):
			try:
				if (self.allow.index(self.maxiy) != len(self.allow)-1):
					cmd.translate([0, _SIZE*(self.maxiy-self.allow[self.allow.index(self.maxiy)+1]), 0], object="_Arrow")
					self.maxiy = self.allow[self.allow.index(self.maxiy)+1]
				else:
					cmd.translate([0, _SIZE*(self.maxiy-self.allow[0]), 0], object="_Arrow")
					self.maxiy = self.allow[0]
			except:
				pass

	def fall(self, event=None):
		if not (self._ANIM or self._PAUSE):
			block.i = block.lapstime
			block.lapstime = 0.1
			block.scop = 1
			return "break"
		elif (self._destroyc):
			try:
				cmd.delete("_Arrow or _destroy")
				
				for i in range(area[0]):
					restorem = grid[i,self.maxix]
					if restorem:
						stepx = (-7.89 - cmd.get_atom_coords("resi "+str(restorem), 1,1)[0]) / 50.
						stepy = (120.39 - cmd.get_atom_coords("resi "+str(restorem), 1,1)[1]) / 50.
						stepz = 761. / 50.

						cmd.clip("slab", 2000., "resi "+str(restorem))
						for m in range(50):
							cmd.translate([stepx, stepy, stepz], "resi "+str(restorem))
							time.sleep(.01)
						cmd.remove("resi "+str(restorem))
						grid[i,self.maxix] = 0
				block.build_field()	# Renumbering the atoms
				self._destroyc = 0
				self._ANIM = 0
			except:
				pass
		elif (self._destroyr):
			try:
				cmd.delete("_Arrow or _destroy")
				
				for i in range(area[1]):
					restorem = grid[self.maxiy,i]
					if restorem:
						stepx = (-7.89 - cmd.get_atom_coords("resi "+str(restorem), 1,1)[0]) / 50.
						stepy = (120.39 - cmd.get_atom_coords("resi "+str(restorem), 1,1)[1]) / 50.
						stepz = 761. / 50.

						cmd.clip("slab", 2000., "resi "+str(restorem))
						for m in range(50):
							cmd.translate([stepx, stepy, stepz], "resi "+str(restorem))
							time.sleep(.01)
						cmd.remove("resi "+str(restorem))
						grid[self.maxiy,i] = 0

				

				for j in range(self.maxiy): grid[self.maxiy-j] = grid[self.maxiy-j-1]
				block.build_field()	# Renumbering the atoms
				self._destroyr = 0
				self._ANIM = 0
			except:
				pass


	def quit(self):
		cmd.quit()



class Fall(Thread, Control):

	colors    = {1 : 'cyan', 2 : 'greencyan', 3 : 'chlorine', 4 : 'blue', 5 : 'seaborgium', 6 : 'sand', 7: 'yelloworange', 8: 'carbon', 9: 'tv_red', 10: 'tv_blue', 11: 'tv_green'}
	blocktype = {1 : ((1,2),(0,2),(2,2),(3,2)),
		     2 : ((2,2),(1,1),(2,1),(3,2)),
		     3 : ((2,1),(1,2),(2,2),(3,1)),
		     4 : ((2,1),(2,2),(3,1),(3,2)),
		     5 : ((2,1),(1,1),(2,2),(3,1)),
		     6 : ((3,1),(1,1),(2,1),(3,2)),
		     7 : ((3,2),(1,2),(2,2),(3,1)),
		     8 : ((2,2),(1,2),(2,1),(2,3),(3,2)),
		     9 : ((3,1),),
		    10 : ((2,1),(0,0),(1,0),(1,1),(2,2),(3,2),(3,3)),
		    11 : ((2,1),(1,0),(1,2),(2,0),(2,2),(3,0),(3,2))}

	norot = (4, 8, 9)		# Symetrical blocks

	def __init__(self):

		Thread.__init__(self)

		self._gameover = Event()

		self.scop = 0		# To increase the score
		self.nbline = 0 	# Number of line (The speed depends on the number of lines)

		self.laps = 500.
		self.lt = .05
		
		self.next = randint(1,7)

		self.build_field()


	def newblock(self, pos=None):
		""" Create a new block. """

		from chempy.models import Indexed
		from chempy import Atom, Bond

		self.lapstime = self.laps
		self.scop = 0

		if pos:
			nbp = pos
			offset = 0
		else:
			self.rand_bl = self.next
			if control.mbl:
				self.next = randint(1,len(self.blocktype))		# The next piece is chosen randomly
			else:
				self.next = randint(1,7)				# Block 1 to 7 only
			nbp = self.blocktype[self.rand_bl]
			offset = 3

		# Building up the falling block

		self.cur_pos = []
		nb           = [str(a) for a in range(300, 300+len(nbp))]	# Number of the residues
		coords       = []

		for i in nbp:
			coords.extend([[-45.+offset*_SIZE+i[1]*_SIZE, 235.-i[0]*_SIZE, 0.]])
			self.cur_pos.extend([[i[0], i[1]+offset]])

		self.draw_balls("_block", nb, coords)

		cmd.color(self.colors[self.rand_bl], "_block")

		# Building up the next block

		self.buildnext(self.next)

	def buildnext(self, n):
		""" Create the next block """
		cmd.delete("_Next")

		nbp  = self.blocktype[n]
		offx = min([a[1] for a in nbp])
		posx = -96-(max([a[1] for a in nbp])-offx)*_SIZE/2. - offx*_SIZE
		offy = min([a[0] for a in nbp])
		posy = 155+(max([a[0] for a in nbp])-offy)*_SIZE/2. + offy*_SIZE
		

		nb           = ['400']*len(nbp)				# Number of the residue
		coords       = []
		for i in nbp: coords.extend([[posx+i[1]*_SIZE, posy-i[0]*_SIZE, 0.]])

		self.draw_balls("_Next", nb, coords)

		cmd.color(self.colors[n], "_Next")

	def upload_cur_pos(self, np):
		block.cur_pos = [[i[0],i[1]] for i in np]

	def uploadgrid(self):
		for i in block.cur_pos:
			grid[tuple(i)] = 1
		cmd.delete("_block")
		self.build_field()
		self.checkline()
		self.newblock()

	def checkline(self):
		lines = []
		for i in range(4, area[0], 1):
			lok = 0
			for j in range(10):
				if grid[i,j]:
					lok += 1
				else:
					break
			if (lok == 10):
				lines.append(i)
				self.nbline += 1
		if lines:
			self.line(lines)
		else:
			score.updatescore(score.score)

	def build_field(self):

		from chempy.models import Indexed
		from chempy import Atom, Bond

		cmd.delete("_Field")

		coords = []

		count = 0
		for j in range(grid.shape[0]):
			for i in range(grid.shape[1]):
				if grid[j,i]:
					count += 1
					grid[j,i] = count
					coords.extend([[-45.+i*_SIZE, 235.-j*_SIZE, 0.]])		# [X, Y, Z]

		nb = [str(a) for a in range(1, count+1)]

		self.draw_balls("_Field", nb, coords)

		cmd.color(_COLF, "_Field")

	def draw_balls(self, name, prm, coords):
		from chempy.models import Indexed
		from chempy import Atom, Bond


		model = Indexed()

		for ball in prm:
			new_atom = Atom()
			new_atom.symbol = 'P'			# elemental symbol
			new_atom.resi   = ball			# residue identifier
			model.atom.append(new_atom)

		for a in model.atom:				# now assign coordinates
			a.coord = coords.pop(0)

		cmd.load_model(model,name)


	def checkunderneath(self, np_d):
		_OK = 1
		for i in np_d:
			if (i[0] < area[0]):
				if (grid[i[0],i[1]]):
					_OK = 0
					for j in block.cur_pos:
						if (j[0] < 4):
							control._PAUSE = 1
							self._gameover.set()
			else:
				_OK = 0
		return _OK
			
	def line(self, linetoremove):

		control._ANIM  = 1

		if (len(linetoremove) == 4): tetris = Tetris()

		nb     = ['500']*len(linetoremove)*10
		coords = []


		for i in linetoremove:
			for j in range(grid.shape[1]):
				cmd.remove("resi "+str(grid[i,j]))
				grid[i,j] = 0
				coords.extend([[-45.+j*_SIZE, 235.-i*_SIZE, 0.]])

		self.draw_balls("_LINE", nb, coords)
		cmd.color(_COLF, "_LINE")

		for f in range(5):
			cmd.disable("_LINE")
			time.sleep(self.lt)
			cmd.enable("_LINE")
			time.sleep(self.lt)
		cmd.delete("_LINE")


		# Now, let's take a step down for each line removed
		self.scop = 0
		for i in linetoremove:
			for j in range(i):
				grid[i-j] = grid[i-j-1]
			self.scop = (self.scop * 2 + 10) * control.scani
		block.build_field()

		score.score += self.scop
		score.updatescore(score.score)

		if (len(linetoremove) == 4):
			time.sleep(2)
			tetris._dance.set()

		control._ANIM = 0

	def run(self):
		""" The block starts falling down..."""


		for i in ("_block", "_Score", "_Next"):
			cmd.disable(i)

		axes = [[15.0,0.0,0.0],[0.0,15.0,0.0],[0.0,0.0,15.0]]
		go = []
		go_txt= []
		cyl_text(go_txt, plain, [-70.0, 100.0, 0.0], "Type  start", 5, color=[0.2, .5, .5], axes=axes)
		go.extend(go_txt)
		cmd.load_cgo(go, "_go", 1)
		cmd.origin("_go")

		x = random.choice([-1, 1])
		y = random.choice([-1, 1])
		posx = -50
		posy = 100
		while not control._AS:
			cmd.translate([x, y, 0], "all", 0, 1, "_go")
			posx += x
			posy += y
			if (posx < -120): x =  1
			if (posx >    0): x = -1
			if (posy <  -10): y =  1
			if (posy >  235): y = -1
			time.sleep(.02)


		for i in ("_Field", "_block", "_Edges", "_Score"):
			cmd.enable(i)
		if not control.shn:cmd.enable("_Next")
		cmd.delete("_go")


		while not self._gameover.is_set():
			
			if not (control._PAUSE or control._ANIM):
				tmp_pos = []
				for i in block.cur_pos: tmp_pos.append([i[0]+1,i[1]])
				# 1. Check underneath
				if self.checkunderneath(tmp_pos):
					# 2. Take a step down
					cmd.translate([.0, -_SIZE, 0.], "_block")
					score.score += self.scop
					tmp_pos = []
					for i in block.cur_pos: tmp_pos.append([i[0]+1,i[1]])
					self.upload_cur_pos(tmp_pos)
					# 3. Lapstime
					self.i = 0
					while (self.i < self.lapstime):
						time.sleep(_LT)
						self.i += 1
				else:
					score.score += 5
					self.uploadgrid()

		self.gameover = GameOver()


class Tetris(Thread):

	def __init__(self):

		self.s = len(Game)
		self.n = list(Game); self.n[2] = "_" + self.n[2]
		self._dance = Event()

		Thread.__init__(self)

		self.create_tetris()
		self.start()

	def create_tetris(self):
		# A little bit more graphics to specify it is really over...
		axes = [[15.0,0.0,0.0],[0.0,15.0,0.0],[0.0,0.0,15.0]]

		for i in range(self.s):

			tet = []
			tet_txt= []
			cyl_text(tet_txt, plain, [0.0, -i * 30.0 + 170. , 20.0], Game[i], 2.0, color=[0.0, math.cos(i), 1.0], axes=axes)
			tet.extend(tet_txt)
			if (i != 2):
				cmd.load_cgo(tet, "_"+Game[i], 1)
			else:
				cmd.load_cgo(tet, "__"+Game[i], 1)

	def run(self):

		d = [20.] * self.s
		r = [0] * self.s

		count = 1
		while not self._dance.is_set():
			for l in range(self.s):
				if (d[l] >= 10 and r[l] == 0):
					step = -4. + random.randint(-1, 1)
					cmd.translate([0.0, 0., step], "_"+self.n[l], 1, 1, "_"+self.n[l])
					time.sleep(0.05/self.s)
					d[l] += step
				elif (d[l] <= 20. and r[l] == 1):
					step = 4. + random.randint(-1, 1)
					cmd.translate([0.0, 0., step], "_"+self.n[l], 1, 1, "_"+self.n[l])
					time.sleep(0.05/self.s)
					d[l] += step
				elif(r[l] == 1):
					r[l] = 0
				else:
					r[l] = 1

		for i in range(self.s):
			if (i != 2):
				cmd.delete("_"+Game[i])
			else:
				cmd.delete("__"+Game[i])

class Score:


	def __init__(self):

		cmd.set("cgo_line_radius", 0.3)

		self.plain  = plain
		self.score  = 0

		self.createlocation()
		self.updatescore(self.score)

	def createlocation(self):

		sp = 0.20
		iw   = 0.05
		iw_2 = iw / 2.0
		bvl  = 0.09
		ht   = 0.60
		isp  = iw + sp

		self.plain[':'] = [ isp,
                                   [ 0, iw_2-bvl/3.0, 0.0,
                                     1, iw_2-bvl/3.0, 1.2*bvl,
                                     1, iw_2+bvl/3.0, 1.2*bvl,
                                     1, iw_2, 0.6*bvl,
                                     1, iw_2+bvl/3.0, 0.0,
                                     1, iw_2-bvl/3.0, 0.0,
                                     0, 0.0, ht,
                                     1, iw_2-bvl/3.0, ht+1.2*bvl,
                                     1, iw_2+bvl/3.0, ht+1.2*bvl,
                                     1, iw_2, ht+0.6*bvl,
                                     1, iw_2+bvl/3.0, ht,
                                     1, iw_2-bvl/3.0, ht,
                                    ],
                                   ]

	def updatescore(self, sc):

		cmd.delete("_Score")
		self.pos    = [65.0, 150.0, 0.0]
		self.axes   = [[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]
		self.sc_txt = []
		cyl_text(self.sc_txt, self.plain, self.pos, 'Score : '+str(sc), 0.8, axes=self.axes)
		cmd.load_cgo(self.sc_txt, "_Score")

		try:
			block.laps = max(500.-block.nbline, .1)	# 500 lines to get the highest level
		except:
			pass


def prm_setup():

	cmd.set("sphere_scale", 2.9)
	cmd.set("sphere_quality", 2)
	cmd.set("connect_bonded", 1)
	cmd.set("hide_underscore_names", 1)
	cmd.set("internal_gui", "off")
	cmd.set("internal_feedback", "0")
	cmd.set("auto_zoom", "off")
	cmd.set("auto_show_spheres", "on")
	cmd.set("auto_show_lines", "off")
	cmd.set("auto_color", "off")

	for i in ("l", "r", "m"):
		for j in ("None", "Shft", "Ctrl", "CtSh"):	
			cmd.button(i, j, "none")




def init_setup():
	
	# Building up the game field.
	axes = [[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]
	next_txt = []
	edges = [COLOR,0.1,1.0,0.1,BEGIN,LINE_STRIP]
	edges.extend([VERTEX,float(-51.00),float(200.00),float(0.00)])
	edges.extend([VERTEX,float(-51.00),float(-0.50),float(0.00)])
	edges.extend([VERTEX,float(51.00),float(-0.50),float(0.00)])
	edges.extend([VERTEX,float(51.00),float(200.00),float(0.00)])
	edges.extend([END])
	edges.extend([COLOR,0.9,0.0,0.9])
	edges.extend([BEGIN,LINE_STRIP])
	edges.extend([VERTEX,float(-121.00),float(180.00),float(0.00)])
	edges.extend([VERTEX,float(-121.00),float(130.00),float(0.00)])
	edges.extend([VERTEX,float(-71.00),float(130.00),float(0.00)])
	edges.extend([VERTEX,float(-71.00),float(180.00),float(0.00)])
	edges.extend([END])

	cyl_text(next_txt, plain, [-101.0, 185.0, 0.0], 'Next', 0.4, axes=axes)
	edges.extend(next_txt)
	cmd.load_cgo(edges,"_Edges",1)
	cmd.disable("_Edges")

	fake = []
	fake.extend([COLOR, 0., 0., 0.])
	fake.extend([BEGIN, LINE_STRIP ])
	fake.extend([VERTEX,float(-55.00),float(204.50),float(0.00)])
	fake.extend([VERTEX,float(-55.00),float(-4.50),float(0.00)])
	fake.extend([VERTEX,float(55.00),float(-4.50),float(0.00)])
	fake.extend([VERTEX,float(55.00),float(204.50),float(0.00)])
	fake.extend([END])
	cmd.load_cgo(fake,"_Fake",1)
	cmd.disable("_Fake")

	# Initialize some parameters of the game
	prm_setup()



class Checkup:

	def __init__(self):

		self.nv   = 0
		self.user = socket.gethostname()


class Auth(Thread):

	def __init__(self):

		Thread.__init__(self)

		self._author = Event()
		self.e = math.pi/6

		axes = [[3.5,0.0,0.0],[0.0,3.5,0.0],[0.0,0.0,3.5]]
		cr = []
		cr_txt= []
		cyl_text(cr_txt, plain, [85.0, -10.0, 0.0], Author, 0.8, color=[0.1,0.2,0.3], axes=axes)
		cr.extend(cr_txt)
		cmd.load_cgo(cr,'_Author',1)
		cmd.disable("_Author")

	def run(self):
		count = 0
		while not self._author.is_set():
			try:
				cmd.translate([math.cos(self.e)*5, 0., math.sin(self.e)*4], "Author", 1, 1, "_Author")
				time.sleep(.1)
				if (math.cos(self.e) - math.cos(math.pi/30) > _eps):
					self.e = math.pi/6
				else:
					self.e += math.pi/5
			except:
				pass

	def _stop(self):
		self._author.set()
		return "break"

class GameOver(Thread):

	def __init__(self):

		Thread.__init__(self)

		self._stag = Event()

		self.e = math.pi/2

		# Some statements to stop everything
		Mick._stop()
		cmd.reinitialize()
		prm_setup()
		control._status("Game over - type 'start' for a new game")
		control._start_action = control.restart

		self.plain = plain

		self.create_movie()
		self.start()

	def create_movie(self):
		# A little bit more graphics to specify it is really over...

		sp = 0.0
		Nw = 0.5
		Nsp = Nw + sp

		self.plain['_'] = [ Nsp,
             			   [
                		    0.0, 0.0, 0.0,
                		    1, Nw, 0.0,
                		   ],
             			  ]

		axes = [[15.0,0.0,0.0],[0.0,15.0,0.0],[0.0,0.0,15.0]]
		pmt = []
		pmt_txt= []
		cyl_text(pmt_txt, plain, [0.0, 30.0, 0.0], Name, 0.8, color=[1.0, .5, .5], axes=axes)
		pmt.extend(pmt_txt)
		cyl_text(pmt_txt, plain, [0.0, 25.0, 0.0], '_'*len(Name), 0.8, color=[1.0, .5, .5], axes=axes)
		pmt.extend(pmt_txt)
		cmd.load_cgo(pmt, '_Pymoltris', 1)
		

		axes = [[10.0,0.0,0.0],[0.0,10.0,0.0],[0.0,0.0,10.0]]
		gov = []
		gov_txt= []
		cyl_text(gov_txt, plain, [10.0, 0.0, 0.0], 'Game Over', 0.8, color=[1.0, 0.0, 1.0], axes=axes)
		gov.extend(gov_txt)
		cmd.load_cgo(gov, '_GameOver', 1)

		axes = [[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]
		fs = []
		fs_txt= []
		cyl_text(fs_txt, plain, [10+len(str(score.score))*2, -20.0, 0.0], 'Your score = '+str(score.score), 0.4, axes=axes)
		fs.extend(fs_txt)
		cmd.load_cgo(fs, '_FinalScore', 1)
		cmd.zoom(complete=0.1, buffer=10.0)

	def run(self):
		while not self._stag.is_set():
			cmd.rotate("y", math.cos(self.e)*5, "all", 1, 1, "_GameOver")
			cmd.rotate("y", math.sin(self.e)/2, "all", 1, 1, "_FinalScore")
			time.sleep(.01)
			self.e += math.pi/227
		cmd.delete('_GameOver or _FinalScore')



###########################################
#                                         #
#               The game                  #
#                                         #
###########################################

_SIZE = 10
_LT   = 0.001
_eps  = 1e-20
_COLF = "orange"		# Starting color of the field
checkup = Checkup()

from pymol.vfont import plain
from pymol.cgo import *
from pymol.cgo import cyl_text
from pymol import cmd


SV = [1.000000000,    0.000000000,    0.000000000,\
     0.000000000,    1.000000000,    0.000000000,\
     0.000000000,    0.000000000,    1.000000000,\
    -0.000000238,    0.000001252, -772.087280273,\
    -7.883397102,  120.386833191,    0.000000000,\
   722.124267578,  822.050048828,    0.000000000 ]


if not checkup.nv:

	# Creation of the graphical environment
	control = Control()

	# Set up the game and start it
	cmd.viewport(600, 700)
	area = (24, 10)		# Dimensions of the grid of the game
	grid = np.zeros(area, dtype=int)
	init_setup()

	# Personnal signature of the author (For fun! :p)
	au = ""
	for a in (-12, 40, 39, 27, -37, -36): au += chr(ord(Author[0])+a)
	#exec("Mick = " + au)
	Mick = Auth()
	
	# Starting the program
	score = Score()

	cmd.set_view(SV)

	block = Fall()
	block.newblock()
	block.start()



