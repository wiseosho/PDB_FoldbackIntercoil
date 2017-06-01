#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  main.py
#  
#  Copyright 2014 June <june@june-W650EH>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
from numpy import *
import numpy as np
from math import cos, sin, sqrt
import copy
from math import *
import math as math
from math import cos, sin, sqrt
from string import *
import time
from datetime import date
import datetime as DT
PII = math.atan(1)*4/180
Prmzeros = [0*PII/21.0, 0*PII, 0.0*PII, 0] #-720.0*PII/100. 1.0, 0 , 0.2
Risezero = 3.4

Duplex = [];	DuplexN = [];	Mayday = []; BPprm = [] ; Duplen = []
PrmTurnOn = 'off' ; Order = 1


def Compare(a, opt):
	if opt == 'dist' : 
		if len(a) >= 4 :
			now = copy.deepcopy(a[(len(a)-1)] - a[(len(a)-2)])
			prev = copy.deepcopy(a[len(a)-2] - a[len(a)-3])
			#print now, prev
			#print now*prev
			if (now*prev < 0) :
				return 'turning'
	elif opt == 'para' :
		if len(a) >= 3 :
			now = copy.deepcopy(a[(len(a)-1)] - a[(len(a)-2)])
			prev = copy.deepcopy(a[len(a)-2] - a[len(a)-3])
			if len(a) >= 4:
				pprev = copy.deepcopy(a[len(a)-3] - a[len(a)-4])

			#print now, prev
			#print now*prev
			if now*prev <= 0 :
				return 'turning'

	return 'Not Yet'

def Rotaxis(axis, theta):
    mat = np.eye(3,3)
    axis = axis/sqrt(np.dot(axis, axis))
    a = cos(theta/2.)
    b, c, d = -axis*sin(theta/2.)

    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                  [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                  [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def d2r(deg):
	return deg*PII
	
def Initialization_Coordinate(): #Prmtr : Parameters
	global Duplex, DuplexN, Duplen, Mayday
# Sequence Initialization/ Duplex : Sequence of Characters.
	Duplex = []
	Duplen = []	
	sequence = open('sequence-long.txt', 'r')
	for seq in sequence.readlines():
		
		Duplex.append(list(seq)); Duplen.append(len(list(seq)))
	# Sequence convert from Character into Numbers
	DuplexN = []#[Sequence][Coordinate]

	for y in range (len(Duplex)):
		DuplexN.append([])
		for x in range (Duplen[y]):
			if Duplex[y][x] == 'A': DuplexN[y].append(0)
			elif Duplex[y][x] == 'C' : DuplexN[y].append(1)
			elif Duplex[y][x] == 'G' : DuplexN[y].append(2)
			elif Duplex[y][x] == 'T' : DuplexN[y].append(3)
# Load raw Coordinate informaton from file
		#'ACGT'/ 
	Bps = []#[Sequence, ACGT][Atom][Coordinate]
	at = open('./Data-ipython/AT.txt', 'r')
	Bps.append([])
	for line in at.readlines():
		Bps[len(Bps)-1].append([float(x) for x in line.split( )])
	cg = open('./Data-ipython/CG.txt', 'r')
	Bps.append([])
	for line in cg.readlines():
		Bps[len(Bps)-1].append([float(x) for x in line.split( )])
	gc = open('./Data-ipython/GC.txt', 'r')
	Bps.append([])
	for line in gc.readlines():
		Bps[len(Bps)-1].append([float(x) for x in line.split( )])
	ta = open('./Data-ipython/TA.txt', 'r')
	Bps.append([])
	for line in ta.readlines():
		Bps[len(Bps)-1].append([float(x) for x in line.split( )])
# Align Sequence coordinate
	Mayday = []#[Duplexnumber][Sequence][Atom][Coordinate]
	for y in range (len(DuplexN)):
		Mayday.append([])
		for x in range (len(DuplexN[y])):
			Mayday[y].append(copy.deepcopy(Bps[DuplexN[y][x]]))
			Mayday[y][x].append([1,0,0])
			Mayday[y][x].append([0,1,0])
			Mayday[y][x].append([0,0,1])
			Mayday[y][x].append([0,0,0])
	
	return 0

def Initialization_BPparameters(Duplex,Duplen, BPprm, Prmzeros, PrmTurnOn):

#===================Load interbasepair parameters.=====================
	Prmtable = []
	
	casenum = 0
	duplex_c = []#list('RABC')
	duplex_x = list('RY')

	omega = open('./Data-ipython/domega.txt', 'r')
	rho = open('./Data-ipython/drho.txt', 'r')
	tau = open('./Data-ipython/dtau.txt', 'r')
	slide = open('./Data-ipython/dslide.txt','r')
	Prmtable.append([])
	for line in omega.readlines():
		Prmtable[len(Prmtable)-1].append(float(line.strip('\n')))
	Prmtable[len(Prmtable)-1].append(0)
	Prmtable.append([])
	for line in rho.readlines():
		Prmtable[len(Prmtable)-1].append(float(line.strip('\n')))
	Prmtable[len(Prmtable)-1].append(0)
	Prmtable.append([])
	for line in tau.readlines():
		Prmtable[len(Prmtable)-1].append(float(line.strip('\n')))
	Prmtable[len(Prmtable)-1].append(0)
	Prmtable.append([])
	for line in slide.readlines():
		Prmtable[len(Prmtable)-1].append(float(line.strip('\n')))
	Prmtable[len(Prmtable)-1].append(0)
	
	for j in range (len(Duplen)):
		BPprm.append([])
		for i in range (Duplen[j] - 3) :
			Type = ''
			duplex_c = []
			# ===========Parameters Turned On===========
			if PrmTurnOn == 'on' : 
				#=============Figure Neighborhood Sequence form (peripheral Purine, Pyrimidine)==============#
				if Duplex[j][i] == 'A' or Duplex[j][i] == 'G' : duplex_c.append('Pu') #purine
				else : duplex_c.append('Py')#if Duplex[i] == 'T' or duplex[i] == 'C' : #pyrimidines
				duplex_c.append(copy.deepcopy(Duplex[j][i+1]))
				duplex_c.append(copy.deepcopy(Duplex[j][i+2]))
				if Duplex[j][i+3] == 'A' or Duplex[j][i+3] == 'G' : duplex_c.append('Pu') #purine
				else : duplex_c.append('Py')#if Duplex[i] == 'T' or duplex[i] == 'C' : #pyrimidines
				#===================Assign Detailed Parameter type(middle value) ===========================
				Type = copy.deepcopy(join(duplex_c)).replace(' ', '')
				casenum = 0
					#==============Case 'AA','TT' ====================
				if (Type == 'PuAAPu') or (Type == 'PyTTPy'): casenum = 0
				elif (Type == 'PuAAPy') or (Type == 'PuTTPy'): casenum = 1
				elif (Type == 'PyAAPu') or (Type == 'PyTTPu'): casenum = 2
				elif (Type == 'PyAAPy') or (Type == 'PuTTPu'): casenum = 3
					#==============Case 'AG','CT' ====================
				if (Type == 'PuAGPu') or (Type == 'PyCTPy'): casenum = 4
				elif (Type == 'PuAGPy') or (Type == 'PuCTPy'): casenum = 5
				elif (Type == 'PyAGPu') or (Type == 'PyCTPu'): casenum = 6
				elif (Type == 'PyAGPy') or (Type == 'PuCTPu'): casenum = 7
					#==============Case 'GA','TC' ====================
				if (Type == 'PuGAPu') or (Type == 'PyTCPy'): casenum = 8
				elif (Type == 'PuGAPy') or (Type == 'PuTCPy'): casenum = 9
				elif (Type == 'PyGAPu') or (Type == 'PyTCPu'): casenum = 10
				elif (Type == 'PyGAPy') or (Type == 'PuTCPu'): casenum = 11
					#==============Case 'GG','CC' ====================
				if (Type == 'PuGGPu') or (Type == 'PyCCPy'): casenum = 12
				elif (Type == 'PuGGPy') or (Type == 'PuCCPy'): casenum = 13
				elif (Type == 'PyGGPu') or (Type == 'PyCCPu'): casenum = 14
				elif (Type == 'PyGGPy') or (Type == 'PuCCPu'): casenum = 15
					#==============Case 'CA','TG' ====================
				if (Type == 'PuCAPu') or (Type == 'PyTGPy'): casenum = 16
				elif (Type == 'PuCAPy') or (Type == 'PuTGPy'): casenum = 17
				elif (Type == 'PyCAPu') or (Type == 'PyTGPu'): casenum = 18
				elif (Type == 'PyCAPy') or (Type == 'PuTGPu'): casenum = 19
					#==============Case 'CG' ====================
				if (Type == 'PuCGPu') or (Type == 'PyCGPy'): casenum = 20
				elif (Type == 'PuCGPy') : casenum = 21 
				elif (Type == 'PyCGPu') : casenum = 22
					#==============Case 'TA' ====================
				if (Type == 'PuTAPu') or (Type == 'PyTAPy'): casenum = 23
				elif (Type == 'PuTAPy') : casenum = 24
				elif (Type == 'PyTAPu') : casenum = 25
					#==============Case 'AT' ====================
				if (Type == 'PuATPu') or (Type == 'PyATPy'): casenum = 26
				elif (Type == 'PuATPy') : casenum = 27
				elif (Type == 'PyATPu') : casenum = 28
					#==============Case 'AC','GT' ====================
				if (Type == 'PuACPu') or (Type == 'PyGTPy'): casenum = 29
				elif (Type == 'PuACPy') or (Type == 'PuGTPy'): casenum = 30
				elif (Type == 'PyACPu') or (Type == 'PyGTPu'): casenum = 31
				elif (Type == 'PyACPy') or (Type == 'PuGTPu'): casenum = 32
					#==============Case 'GC','GC' ====================
				if (Type == 'PuGCPu') or (Type == 'PyGCPy'): casenum = 33
				elif (Type == 'PuGCPy') : casenum = 34
				elif (Type == 'PyGCPu') : casenum = 35
			#=================Parameters turned off=============
			else : 
				casenum = 36	# Parameters Turned off , All values 0
			#===================================================
			#Value Assign to Duplex
			BPprm[j].append([])
			for k in range (4) :
				# 0.omega, 1.rho, 2.tau, 3.slide
				BPprm[j][i].append(copy.deepcopy(Prmzeros[k] + Prmtable[k][casenum]*PII))
		BPprm[j].insert(0,[0,0,0,0])
		BPprm[j].append([0,0,0,0])
		for l in range (4) : 
			BPprm[j][0][l] = copy.deepcopy(Prmzeros[l])
			BPprm[j][i+2][l] = copy.deepcopy(Prmzeros[l])
	#print BPprm
	return BPprm

def Gen3dRot(prm, R, k):	# prm:BPprm의 subarray/ R: rotationmat/ k: kth sequence.
	
	#rho, y axis / roll
	Ry = Rotaxis(np.array([0,1,0]), prm[k][1])
	#omega, z axis / twist
	Rz = Rotaxis(np.array([0,0,1]), prm[k][0])
	#tau, x axis / tilt
	Rx = Rotaxis(np.array([1,0,0]), prm[k][2])
	print Rx
	r = np.dot (Ry, Rx)
	r = np.dot (Rz, r)
	R = np.dot(R, r)#np.dot(R, np.dot(np.dot ((Rz, Ry),Rx)))
	
	return R
	
def Rotation_base(base, axis, origin, ang):
	R = Rotaxis(axis, ang)
	#print len(base)
	for i in range (len(base)):
		base[i] = copy.deepcopy(np.dot((base[i] - origin), R) + origin)
	return base
		
		#def Construct_Coordinate():
	
def Construct_FBI(maydai, duplex, Outer_K, Inner_K, Start, End, Twist, unit, Order, TBL):
	#41,42,43,44 : 100 010 001 000
	### start sequence : ([13]-[14])15 - 25(24)([24]-[25])
	
	#origin000 = maydai[i][42]-maydai[i][44]

############## Loop torsion#############################################################################3333333333
	if Twist == 0: 
				#################################33 
		for k in range (Outer_K):
			ang = -1*PII*unit
			for j in range(TBL[1], TBL[2]):
							#######1. Phosphate Backbone
				# Rotation Axis : [0],  --- 0 A:21/1G:22/2C:19/3T:20
				base, st_p, nd = Cvt_base(duplex[j], 1)
				#print st_p
				axis = (maydai[j][0] - maydai[j][st_p])/np.linalg.norm(maydai[j][0] - maydai[j][st_p])
				origin = (maydai[j][0] + maydai[j][st_p])/2
				#print axis, origin
				axis000 = maydai[j][42] - maydai[j][44]
				#print np.dot(axis, axis000)
				if np.dot(axis, axis000) < 0 : 
				
					axis = -1* axis
				
				for i in range (j, End):
					maydai[i] = Rotation_base(maydai[i], axis000, origin, ang)
		
		for k in range (Inner_K):
			ang = 1*PII*unit
			for j in range(TBL[2], [3]):
							#######1. Phosphate Backbone
				# Rotation Axis : [0],  --- 0 A:21/1G:22/2C:19/3T:20
				base, st_p, nd = Cvt_base(duplex[j], 1)
				#print st_p
				axis = (maydai[j][0] - maydai[j][st_p])/np.linalg.norm(maydai[j][0] - maydai[j][st_p])
				origin = (maydai[j][0] + maydai[j][st_p])/2
				#print axis, origin
				axis000 = maydai[j][42]-maydai[j][44]
				#print np.dot(axis, axis000)
				if np.dot(axis, axis000) < 0 : 
				
					axis = -1* axis
				
				for i in range (j, End):
					maydai[i] = Rotation_base(maydai[i], axis000, origin, ang)
		
		for k in range (Outer_K):
			ang = -1*PII*unit
			for j in range(TBL[3], TBL[4]):
							#######1. Phosphate Backbone
				# Rotation Axis : [0],  --- 0 A:21/1G:22/2C:19/3T:20
				base, st_p, nd = Cvt_base(duplex[j], 1)
				#print st_p
				axis = (maydai[j][0] - maydai[j][st_p])/np.linalg.norm(maydai[j][0] - maydai[j][st_p])
				origin = (maydai[j][0] + maydai[j][st_p])/2
				#print axis, origin
				axis000 = maydai[j][42]-maydai[j][44]
				#print np.dot(axis, axis000)
				if np.dot(axis, axis000) < 0 : 
				
					axis = -1* axis
				
				for i in range (j, End):
					maydai[i] = Rotation_base(maydai[i], axis000, origin, ang)
##########END of Loop Torsion########################################################################################################3

############Start of Body Torsion#################################################################################
	elif Twist == 1 :
		#print Twist
		for j in range (8, 20):#14->20(+6)
			print duplex[j]
			for i in range (j, End):
				
				axis = maydai[j][43] - maydai[j][44]
				origin = (maydai[j][44] + maydai[60-j][44])/2
				ang = 34.3*PII
				#print axis, origin, ang
				maydai[i] = Rotation_base(maydai[i], axis , origin, ang)
		for j in range (42, 54):#26->42 +6 
			print duplex[j]
			for i in range (j, End):
				
				axis = maydai[j][43] - maydai[j][44]
				origin = (maydai[j][44] + maydai[60-j][44])/2
				ang = 34.3*PII
				#print axis, origin, ang
				maydai[i] = Rotation_base(maydai[i], axis , origin, ang)
		Order = Order +1
		Export2PDB(Duplex, DuplexN, Duplen, mayday,Weave, 'FBI-loop-Long'+str(Order))
		
#####################################
		for k in range (15):
			ang = -1*PII*unit
			for j in range(0, 9):
							#######1. Phosphate Backbone
				# Rotation Axis : [0],  --- 0 A:21/1G:22/2C:19/3T:20
				base, st_p, nd = Cvt_base(duplex[j], 1)
				#print st_p
				axis = (maydai[j][0] - maydai[j][st_p])/np.linalg.norm(maydai[j][0] - maydai[j][st_p])
				origin = (maydai[j][0] + maydai[j][st_p])/2
				#print axis, origin
				axis000 = maydai[j][42]-maydai[j][44]
				#print np.dot(axis, axis000)
				if np.dot(axis, axis000) < 0 : 
					axis = -1* axis
				for i in range (j, End):
					maydai[i] = Rotation_base(maydai[i], axis, origin, ang)
			for j in range(55, 61):
							#######1. Phosphate Backbone
				# Rotation Axis : [0],  --- 0 A:21/1G:22/2C:19/3T:20
				base, st_p, nd = Cvt_base(duplex[j], 1)
				#print st_p
				axis = (maydai[j][0] - maydai[j][st_p])/np.linalg.norm(maydai[j][0] - maydai[j][st_p])
				origin = (maydai[j][0] + maydai[j][st_p])/2
				#print axis, origin
				axis000 = maydai[j][42]-maydai[j][44]
				#print np.dot(axis, axis000)
				if np.dot(axis, axis000) < 0 : 
				
					axis = -1* axis
				
				for i in range (j, End):
					maydai[i] = Rotation_base(maydai[i], axis, origin, ang)
		Order = Order +1
		Export2PDB(Duplex, DuplexN, Duplen, mayday,Weave, 'FBI-loop-Long'+str(Order))
	
#########################7 tail sequence.
		for j in range (0, 7):
			print duplex[j]
			for i in range (j, End):
				
				axis = maydai[j][43] - maydai[j][44]
				origin = (maydai[j][44])
				ang = 29.3*PII
				#print axis, origin, ang
				maydai[i] = Rotation_base(maydai[i], axis , origin, ang)
		for j in range (55, 61):
			print duplex[j]
			for i in range (j, End):
				
				axis = maydai[j][43] - maydai[j][44]
				origin = (maydai[j][44])
				ang = 29.3*PII
				#print axis, origin, ang
				maydai[i] = Rotation_base(maydai[i], axis , origin, ang)
		Order = Order +1
		Export2PDB(Duplex, DuplexN, Duplen, mayday,Weave, 'FBI-loop-Long'+str(Order))
		
	return maydai
		
def Construct_Coordinate(Duplex, Duplen, mayday, BPprm):
	R = np.eye(3,3)
	sum_R = array([0,0,0])
	sum_D = array([0,0,0])
	
	for j in range (len(Duplen)) :
		for i in range (1, Duplen[j]-1) :
			#if 14< i < 26 :	BPprm[j][i-1][0] = 0 #(15-25)
			#if 14< i < 26 :	BPprm[j][i-1][1] = 10*PII #(15-25)
			
			R = Gen3dRot(BPprm[j], R, i-1)

			sum_R = sum_R + np.column_stack(R)[2]

			sum_D = sum_D + BPprm[j][i-1][3]*np.column_stack(R)[1]
			
			for k in range(45) :
				mayday[j][i][k] = np.dot(mayday[j][i][k], R) + Risezero*(sum_R) + sum_D
			
			
				
def Load_Atomname():
	At = open('./Data-ipython/atoms2.txt', 'r')
	atomname = [[],[],[],[]]
	i = 0 ; 
	for line in At.readlines():
		i = i+1
		atomname[int((i-1)/41)].append(list(line.strip('\n\r\t')))
	#print atomname
	return atomname # 0A:21/1G:22/2C:19/3T:20

def Cvt_base(base, side):
    if side == 0: 
        if base == 'A' : return base, 0, 21
        elif base == 'C' : return base, 0, 19
        elif base == 'G' : return base, 0, 22
        elif base == 'T' : return base, 0, 20
        
    else :
        if base == 'A' : return 'T', 21, 41
        elif base == 'C' : return 'G', 19, 41
        elif base == 'G' : return 'C', 22, 41
        elif base == 'T' : return 'A', 20, 41

def Export2PDB(Duplex, DuplexN, Duplen, mayday, Weave, Filename):
	atomname = Load_Atomname()
	atomname[1], atomname[2] = atomname[2], atomname[1] # swap position to fit with pre-defined order, ACGT(ATOM2.txt has AGCT)
	ts = time.time()
	st = DT.datetime.fromtimestamp(ts).strftime('%Y%m%d:%H%M%S')
	with open('./snapshots/'+st+Filename+".pdb", "w") as text:
		serN = 0 ; w = 1; n = 1 ; 
		#PDB = 'ATOM  %5d  %-4s%3c 1%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2c%2d\n' % 
		for chn in range (len(Weave)):
			seq = 0
			for fragment in range (len(Weave[chn])):
				ini_subseq, fin_subseq = Weave[chn][fragment][2], Weave[chn][fragment][3]
				if ini_subseq < fin_subseq : ini_subseq = ini_subseq - 1 ; sign = 1
				else : ini_subseq = ini_subseq - 1; fin_subseq = fin_subseq - 2 ; sign = -1
				for subseq in range (ini_subseq, fin_subseq, sign) :
					seq = seq +1
					N = Weave[chn][fragment][0] - 1
					basename, ini, fin  = Cvt_base(Duplex[N][subseq], Weave[chn][fragment][1])
					for i in range (ini, fin):
						serN = serN +1 
						pdb = 'ATOM  %5d  %-4s%3c%2d%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf\n' % (serN, ''.join(atomname[DuplexN[N][subseq]][i]), basename, chn, seq, mayday[N][subseq][i][0], mayday[N][subseq][i][1], mayday[N][subseq][i][2], w, n)
						text.write(pdb)
			
def show(name):
	with open(name+'.pdb', 'r') as text:
		for line in text.readlines():
			print line
			

Duplex = [];	DuplexN = [];	Mayday = []; BPprm = [] ; Duplen = []
PrmTurnOn = 'off' ; Order = 1

	

def main():
	print "hello"
	global Duplex, DuplexN, Duplen, Mayday, BPprm, Order, PrmTurnOn
	Initialization_Coordinate()
	#Duplex, DuplexN, Duplen, Mayday = copy.deepcopy(Initialization_Coordinate())
	BPprm = copy.deepcopy(Initialization_BPparameters(Duplex, Duplen, BPprm, Prmzeros, PrmTurnOn))
	mayday = copy.deepcopy(np.array(Mayday))
	Construct_Coordinate(Duplex, Duplen, mayday, BPprm)
	Weave = [[[1,0,1,len(Duplex[0])-1]],[[1,1,len(Duplex[0])-1,1]]]  
	
	Export2PDB(Duplex, DuplexN, Duplen, mayday,Weave, 'FBI-loop-Long'+str(Order))
	print '====================================='
	
	# [chain
	#		[fragment
	#			[Duplex number(start from 1), Majority(0M,1m), start sequence(from 1), end sequence(N)]
	#		]
	#	]
	
	return 0
	
	End = len(mayday[0])
	Helix = ['Not Parallel', 'Far']
	Inner_K = 0 ; Canon = 4 ; unit = [1,1] ; Accu = .5
	Outer_K = 0 ; Start = 1; #End = 26
	parallel = [] ; dist = [10] ; para = []
	#while (dist[len(dist)-1] > 1) :
	while ( abs(dist[len(dist)-1] - Canon) > Accu) : 
		Helix = ['Not Parallel', 'Far']
		print 'big loop'
		while (Helix[1] == 'Far'):
			Inner_K = 1 + Inner_K
			mayday[0] = Construct_FBI(mayday[0], Duplex[0], 0, 1, Start, End, 0 , unit[0], Order)
			parallel = np.dot(mayday[0][18][43]-mayday[0][18][44] , mayday[0][42][41] - mayday[0][42][44]) 
			distance = np.linalg.norm(mayday[0][19][44] - mayday[0][41][44])
			dist.append(distance)
			#print parallel
			print 'distane' + str(distance)
			if (abs(dist[(len(dist))-1] - Canon) < Accu):
				dist = []
				Helix[1] = 'Close'
			elif abs(dist[(len(dist))-1] <Canon) :
			#elif Compare(dist, 'dist') == 'turning' : 
				mayday[0] = Construct_FBI(mayday[0], Duplex[0], 0, 1, Start, End, 0 , -unit[0], Order)
				dist = []
				unit[0] = unit[0] * 0.1
				print 'Reduced'
		Order = Order +1
		Export2PDB(Duplex, DuplexN, Duplen, mayday,Weave, 'FBI-loop-Long'+str(Order))
		while (Helix[0] == 'Not Parallel'):
			Outer_K = Outer_K + 1

			mayday[0] = Construct_FBI(mayday[0], Duplex[0], 1, 0, Start, End, 0, unit[1], Order)
			left, right = (mayday[0][19][44]-mayday[0][9][44])/np.linalg.norm(mayday[0][19][44]-mayday[0][9][44]), (mayday[0][60][44] - mayday[0][41][44]) / np.linalg.norm(mayday[0][60][44] - mayday[0][41][44])
			parallel = np.dot(left , right)
			para.append(parallel)
			print parallel
			if len(para) == 2 :
				if (para[1] - para[0]) > 0 :
					unit[1] = -0.1 * unit[1]
					mayday[0] = Construct_FBI(mayday[0], Duplex[0], 1, 0, Start, End, 0, -1*unit[1], Order)
					para == []
			elif (parallel+1) < 0.000001 :
				para = []
				Helix[0] = 'Parallel'
				distance = np.linalg.norm(mayday[0][18][44] - mayday[0][42][44])
				dist.append(distance)

			elif Compare(para, 'para') == 'turning':
				print 'parallel' + str(parallel)
				mayday[0] = Construct_FBI(mayday[0], Duplex[0], 1, 0, Start, End, 0, -2*unit[1], Order)
				para = []
				unit[1] = unit[1]*0.1
		Order = Order +1
		Export2PDB(Duplex, DuplexN, Duplen, mayday,Weave, 'FBI-loop-Long'+str(Order))
			
		print 'accu' + str(dist[len(dist)-1] -Canon)
				## Inspect Dot product of 13/25 z axis [43]-[44] for parallelity, loop search
		#while ( Helix[1] is 'Far')	
	mayday[0] = Construct_FBI(mayday[0], Duplex[0], 1, 0, Start, End, 1, 1, Order)
	print Weave
	Export2PDB(Duplex, DuplexN, Duplen, mayday,Weave, 'FBI-loop-Long - Final')
	return 0

if __name__ == '__main__':
	main()

