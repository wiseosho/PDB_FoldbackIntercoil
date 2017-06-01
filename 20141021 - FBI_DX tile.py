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
Prmzeros = [720*PII/21.0, 0*PII, 0.0*PII, 0] #-720.0*PII/100. 1.0, 0 , 0.2
Risezero = 3.4

Duplex = [];	DuplexN = [];	Mayday = []; BPprm = [] ; Duplen = []
PrmTurnOn = 'off' ; Order = 0


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
	Dict = {}
	Footnote = 0
	sequence = open('sequence.txt', 'r')
	for seq in sequence.readlines():
		if Footnote == 0 :
			if seq[0:4] == 'Name' :
				name = seq.split('.')[1]
				Dict[name] = []
			elif seq[0:2] == '//' :
				#print seq.split(' ')
				Dict[name].append(seq.split(' ')[1].strip('\n')) ## 4 spaces.
			elif seq[0:3] == 'END':
				break
		elif seq[0:8] == 'REMARK//':
			Footnote == 1
		elif seq[0:8] == '//REMARK':
			Footnote == 0
		
				
	for val in Dict.values():
		Duplex.extend(list(val))
	for i in range(len(Duplex)):
		Duplex[i] = list(Duplex[i])
		Duplen.append(len(Duplex[i]))
	
	#Duplex.reverse()
	#Duplen.reverse()
	


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

def Initialization_BPparameters():
	global Duplex, Duplen, BPprm, Prmzeros, PrmTurnOn, Prmzeros
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
				if j == 0 : Prmzeros[0] = 0	# j==0 : 1st Duplex, FBI
				BPprm[j][i].append(copy.deepcopy(Prmzeros[k] + Prmtable[k][casenum]*PII))
				Prmzeros[0] = 720/21*PII
				
		BPprm[j].insert(0,[0,0,0,0])
		BPprm[j].append([0,0,0,0])
		for l in range (4) : 
			BPprm[j][0][l] = copy.deepcopy(Prmzeros[l])
			BPprm[j][i+2][l] = copy.deepcopy(Prmzeros[l])
	#print BPprm
	return 0

def Gen3dRot(prm, R, k):	# prm:BPprm의 subarray/ R: rotationmat/ k: kth sequence.
	
	#rho, y axis / roll
	Ry = Rotaxis(np.array([0,1,0]), prm[k][1])
	#omega, z axis / twist
	Rz = Rotaxis(np.array([0,0,1]), prm[k][0])
	#tau, x axis / tilt
	Rx = Rotaxis(np.array([1,0,0]), prm[k][2])
	#print Rx
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
	
def Construct_FBI(maydai, duplex, Operation, K, unit,TBL):
	#print type (maydai[0][43])
	global snshot,  Order
	#41,42,43,44 : 100 010 001 000
	### start sequence : ([13]-[14])15 - 25(24)([24]-[25])
	#origin000 = maydai[i][42]-maydai[i][44]

	End = len(duplex)
	#print End
	Rng = [0,0]
	
############## Loop torsion#############################################################################3333333333	
	if list(Operation)[0] == 'L' :
		Rng[0] = range(TBL[2]+TBL[6]+1, TBL[3]-1)
		Rng[1] = range(TBL[2]+1, TBL[2] + TBL[6]+1) + range(TBL[3]-1, TBL[3]  + TBL[6] -1)
		if Operation == 'Loop-In': m = 0; ang = 1*PII*unit
		elif Operation == 'Loop-Out' : m = 1 ; ang = -1*PII*unit
				#################################33 
		for k in range (K[m]):
			for j in Rng[m]:
							#######1. Phosphate Backbone
				# Rotation Axis : [0],  --- 0 A:21/1G:22/2C:19/3T:20
				base, st_p, nd = Cvt_base(duplex[j], 1)
				axis = (maydai[j][0] - maydai[j][st_p])/np.linalg.norm(maydai[j][0] - maydai[j][st_p])
				origin = (maydai[j][0] + maydai[j][st_p])/2
				axis000 = maydai[j][42] - maydai[j][44]
				#print np.dot(axis, axis000)
				if np.dot(axis, axis000) < 0 : 
				
					axis = -1* axis
				
				for i in range (j, End):
					maydai[i] = Rotation_base(maydai[i], axis000, origin, ang)
		return maydai
##########END of Loop Torsion########################################################################################################3
	elif Operation == 'Body':
		print 'Body'
		rng = range(TBL[1]+2, TBL[2]+1) + range(TBL[3]+2, TBL[4]+1)
		print rng, End
		#return maydai
############Start of Body Torsion#################################################################################
		for k in range(343):
			for j in rng :#14->20(+6)
				for i in range (j, End):
					axis = maydai[j][43] - maydai[j][44]
					origin = (maydai[j][44] + maydai[End-j-1][44])/2
					ang = 0.1*PII
					#print axis, origin, ang
					maydai[i] = Rotation_base(maydai[i], axis , origin, ang)
			if (snshot == 'Y') : Export2PDB('FBI-loop-short')
			print k
		print 'end for body'			
		return maydai
		
#####################################
#######################################Tail Torsion###################################
	elif Operation == 'Tail':
		
		#########################7 tail sequence.
		Rng = range (TBL[0], TBL[1]+1) + range (TBL[4]+1, TBL[5]+1)
		for k in range(293):
			for j in Rng :
				for i in range (j, End):
					axis = maydai[j][43] - maydai[j][44]
					#print j, type(maydai[j][43])
					#print type(maydai[j][44])
					origin = (maydai[j][44])
					ang = 0.1*PII
					maydai[i] = Rotation_base(maydai[i], axis , origin, ang)
			if (snshot == 'Y') : Export2PDB('FBI-loop-short')
		#########################################################################################################
		Rng = range(TBL[0]+1, TBL[1]+2) + range(TBL[4]+1, TBL[5]+1)
		for k in range (15):
			ang = -1*PII*unit
			for j in Rng:
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
			if (snshot == 'Y') : Export2PDB('FBI-loop-short')
		#########################################################################################################
		print 'End for Tail'
	return maydai
		
def Construct_Coordinate():
	global Duplex, Duplen, mayday, BPprm
	
	for j in range (len(Duplen)) :
		R = np.eye(3,3)
		sum_R = array([0,0,0])
		sum_D = array([0,0,0])
		for i in range (1, Duplen[j]) :
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

def Export2PDB(Filename):
	import os
	global Duplex, DuplexN, Duplen, mayday, Weave, imel, Order
	
	Order = Order +1
	atomname = Load_Atomname()
	atomname[1], atomname[2] = atomname[2], atomname[1] # swap position to fit with pre-defined order, ACGT(ATOM2.txt has AGCT)
	ts = time.time()
	daytime = DT.datetime.fromtimestamp(ts).strftime('%Y%m%d - %H%M%S')
	day = DT.datetime.fromtimestamp(ts).strftime('%Y%m%d')
	ime = DT.datetime.fromtimestamp(ts).strftime('%H%M%S')
	
	imel.append(ime)
	if not os.path.exists(day):
		os.makedirs(day)
	if not os.path.exists('./'+day+'/'+imel[0]):
		os.makedirs('./'+day+'/'+imel[0])
	print './'+day+'/'+imel[0]+'/'+Filename+str(Order)+".pdb"
	with open('./'+day+'/'+imel[0]+'/'+Filename+str(Order)+".pdb", "w") as text:
		serN = 0 ; w = 1; n = 1 ; 
		#PDB = 'ATOM  %5d  %-4s%3c 1%4d    %8.3lf%8.3lf%8.3lf%6.2lf%6.2lf          %2c%2d\n' % 
		for chn in range (len(Weave)):
			#print chn
			seq = 0
			for fragment in range (len(Weave[chn])):
				ini_subseq, fin_subseq = Weave[chn][fragment][2], Weave[chn][fragment][3]
				if ini_subseq < fin_subseq : ini_subseq = ini_subseq - 1 ; sign = 1
				else : ini_subseq = ini_subseq - 1; fin_subseq = fin_subseq - 2 ; sign = -1
				#print fragment
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
			
def Construct_FBI_Loop1(TBL):
	global snshot,  mayday, Duplex, Weave, Order, snshot
	
	End = len(mayday[0])
	#print type(mayday[0][0][43]), End
	Helix = ['Not Parallel', 'Far']
	K = [1,1] # [inner, outer]
	unit = [1,1] # [inner, outer] 
	Accu = .5 ; Canon = 4
	Operation = '' # 'Loop-In'/'Loop-Out'/'Body'
	# Unit : Radius Factor for [Inner K, Outer K]
	parallel = [] ; dist = [10] ; para = []
	#while (dist[len(dist)-1] > 1) :
	
	#TailBodyLoop = [0,7,13,24,30,38, 2]
	#TailBodyLoop = []
	# [7,18,22,24,37,
	# 0,1,2,3,4,5,6,7 | 8,9,10,11,12,13 | 14,15,16,17,18 , 19,20,21,22,23,24 | 25,26,27,28,29,30 | 31,32,33,34,35,36,37,38
	# 0,1,2,3,4,5,6,7 | 8,9,10,11,12,13 , 14,15,16,17,18 | 19,20,21,22,23,24   25,26,27,28,29,30   31,32,33,34,35,36,37,38,39,40 | 41,42,43,44,45,46,47,48,49,50,51 | 52,53,54,55,56,57,58,59
	
	####################Rotate Inner Circled Loop till Close.
	while ( abs(dist[len(dist)-1] - Canon) > Accu) : 
		Helix = ['Not Parallel', 'Far']
		while (Helix[1] == 'Far'):
			#Inner_K = 1 + Inner_K
			mayday[0] = Construct_FBI(mayday[0], Duplex[0], 'Loop-In', K, unit[0], TBL)
			if (snshot == 'Y') : Export2PDB('FBI-loop-short')
			parallel = np.dot(mayday[0][TBL[2]][43]-mayday[0][TBL[2]][44] , mayday[0][TBL[3]+1][43] - mayday[0][TBL[3]+1][44]) 
			distance = np.linalg.norm(mayday[0][TBL[2]][44] - mayday[0][TBL[3]+1][44])
			#raw_input(str(distance))
		
			dist.append(distance)
			if (abs(dist[(len(dist))-1] - Canon) < Accu):
				dist = [10]
				Helix[1] = 'Close'
				continue
			elif abs(dist[(len(dist))-1] < Canon) or (Compare(dist, 'dist') == 'turning'):	#Too close or Past
			#elif Compare(dist, 'dist') == 'turning' : 
				mayday[0] = Construct_FBI(mayday[0], Duplex[0], 'Loop-In', K, -unit[0], TBL)
				if (snshot == 'Y') : Export2PDB('FBI-loop-short')
				dist = []
				unit[0] = unit[0] * 0.1
				print 'Reduced'
				continue
			#raw_input()


		if (snshot == 'Y') : Export2PDB('FBI-loop-short')
		
		print 'good'
		####################Rotate Outer Circled Loop till parallel
		while (Helix[0] == 'Not Parallel'):
			mayday[0] = Construct_FBI(mayday[0], Duplex[0], 'Loop-Out', K, unit[1], TBL)
			if (snshot == 'Y') : Export2PDB('FBI-loop-short')
			left, right = (mayday[0][TBL[2]][44]-mayday[0][TBL[0]+2][44])/np.linalg.norm(mayday[0][TBL[2]][44]-mayday[0][TBL[0]+2][44]), (mayday[0][TBL[5]-2][44] - mayday[0][TBL[3]+1][44]) / np.linalg.norm(mayday[0][TBL[5]-2][44] - mayday[0][TBL[3]+1][44])
			parallel = np.dot(left , right)
			para.append(parallel)
			#print parallel
			if (len(para) >= 2) and (para[1] - para[0]) > 0 :
				mayday[0] = Construct_FBI(mayday[0], Duplex[0], 'Loop-Out', K, -unit[1], TBL)
				if (snshot == 'Y') : Export2PDB('FBI-loop-short')
				para == []
				unit[1] = -1*unit[1] # If direction is wrong, change rotation direction
				#raw_input('(len(para) >= 2) and (para[1] - para[0]) > 0 ')
				continue
			elif abs(parallel +1) < 0.000001 :
				para = []
				Helix[0] = 'Parallel'
				distance = np.linalg.norm(mayday[0][TBL[2]][44] - mayday[0][TBL[3]+1][44])
				dist.append(distance)
				
		#		raw_input('abs(parallel) < 1.000001')
				continue

			elif Compare(para, 'para') == 'turning':
				#print 'parallel' + str(parallel)
				mayday[0] = Construct_FBI(mayday[0], Duplex[0], 'Loop-Out', K, -unit[1], TBL)
				if (snshot == 'Y') : Export2PDB('FBI-loop-short')
				para = []
				unit[1] = unit[1]*0.1
		#		raw_input('Compare(para, ''parapara'') == turning')
		#raw_input('Parallel!')
			#Export2PDB('FBI-loop-Long')			
		print 'accu' + str(dist[len(dist)-1] -Canon)
				## Inspect Dot product of 13/25 z axis [43]-[44] for parallelity, loop search
		#while ( Helix[1] is 'Far')	
	
def Rotate_Dup(N):
	global Duplex, Duplen, mayday
	End = len(Duplex[N])
	for j in END :
		for i in range (j, End):
			axis = mayday[N][j][43] - mayday[N][j][44]
			origin = mayday[N][j][44]
			ang = (720/21)*PII
			#print axis, origin, ang
			mayday[N][i] = Rotation_base(mayday[N][i], axis , origin, ang)
	return 0

def Move_Duplex(Dup, unit_dir, quant):
	global mayday
	for seq in range(len(mayday[Dup])) : 
		for atm in range(len(mayday[Dup][seq])):
			mayday[Dup][seq][atm] = mayday[Dup][seq][atm] + unit_dir*quant
	return 0
	
		

def DX_Move_Duplex(jct):
	global Duplex, DuplexN, Mayday, BPprm, Duplen, mayday
	rad = 22
	#Cvt_base => (, residue), residue: 1 : complementary
	#Duplex[N][seq]
	for i in range (2) :
		print Cvt_base(Duplex[1][jct[i]],1)[1]
		Mv_dir = mayday[1][jct[i]][Cvt_base(Duplex[1][jct[i]],0)[1]+8] - mayday[2][jct[i]][Cvt_base(Duplex[2][jct[i]],1)[1]+8] + mayday[1][jct[3-i]][Cvt_base(Duplex[1][jct[3-i]],1)[1]+8] - mayday[2][jct[3-i]][Cvt_base(Duplex[2][jct[3-i]],0)[1]+8]
		Mv_dir = Mv_dir/np.linalg.norm(Mv_dir)
	
	for dup in [1,2] :
		if dup == 1: Dir = -1
		else : Dir = 1
		Move_Duplex(dup, Mv_dir, rad/2*Dir)
				
	return 0
def FBI_rotate2normal():
	global Duplex, Duplen, mayday, BPprm
	#chain4, T20/A21 Duplex[3][19,20][residue] :: chain3, G25/G26 Duplex[3][17,16][complement]
	Dir = []
	Dir.append(np.array(mayday[2][19][Cvt_base(Duplex[2][19],0)[0]] - mayday[2][19][44]))
	Dir.append(np.array(mayday[2][20][Cvt_base(Duplex[2][20],0)[0]] - mayday[2][20][44]))
	Dir.append(np.array(mayday[2][17][Cvt_base(Duplex[2][17],1)[0]] - mayday[2][17][44]))
	Dir.append(np.array(mayday[2][16][Cvt_base(Duplex[2][16],1)[0]] - mayday[2][16][44]))
	print Dir
	#[oxygen point]-mayday[#][residue][44]
	return 0
	
	pass
	
def FBI_Translate():
	pass

def FBI_Fine_adjust():
	pass
	

Duplex = [];	DuplexN = [];	Mayday = []; BPprm = [] ; Duplen = [] ; mayday = []
PrmTurnOn = 'off' ; Order = 1
imel = []
Weave = []
snshot = 'N'

def main():
	print "hello"
	global Duplex, DuplexN, Duplen, Mayday, BPprm, Order, PrmTurnOn, mayday, Weave, snshot
	Initialization_Coordinate()
	Initialization_BPparameters()
	mayday = []
	for i in range(len(Mayday)):
		mayday.append(np.array(Mayday[i]))
	
	Weave = [ [[1,0,1,len(Duplex[0])]],[[1,1,len(Duplex[0]),1]] ]
	# , [[2, 0, 1, jct[0]+1], [3,1,jct[0]+1, 1]], [[3,1, len(Duplex[2]), jct[1]+1], [2, 0, jct[1]+1, len(Duplex[2])]], [[3, 0, 1, jct[2]+1], [2, 1, jct[2]+1, 1]], [[2, 1, len(Duplex[1]), jct[3]+1], [3, 0, jct[3]+1, len(Duplex[1])]] ]#,[[2,1,len(Duplex[1]), 1]], [[3,0,1,len(Duplex[2])]],[[3,1,len(Duplex[2]), 1]]]
	Construct_Coordinate()
	if (snshot == 'Y') : Export2PDB('FBI-loop-short')
	jct = [12, 13, 28, 29]
	
	
	print '====================================='
	
	# [chain
	#		[fragment
	#			[Duplex number(start from 1), Majority(0M,1m), start sequence(from 1), end sequence(N)]
	#		]
	#	]
	# Export2PDB('FBI-loop-Long')
	# return 0
	
	#============================Make DX tile=================================
	#DX_Move_Duplex(jct)
	#Move_Duplex(0, np.array([1,1,0]), 30)
	
	#============================Construct FBI================================
	TailBodyLoop = [0,7,13,24,30,38, 2]
	# [0,0,20,24,37,
	# 0,1,2,3,4,5,6,7 | 8,9,10,11,12,13 | 14,15,16,17,18 , 19,20,21,22,23,24 | 25,26,27,28,29,30 | 31,32,33,34,35,36,37,38
	# 0,1,2,3,4,5,6,7 | 8,9,10,11,12,13 , 14,15,16,17,18 | 19,20,21,22,23,24
	Construct_FBI_Loop1(TailBodyLoop)
	
	mayday[0] = Construct_FBI(mayday[0], Duplex[0], 'Body', [0,0], 0, TailBodyLoop)
	mayday[0] = Construct_FBI(mayday[0], Duplex[0], 'Tail', [0,0], 1, TailBodyLoop)
	#============================Construct FBI================================
	#================================Join FBI=================================
	#FBI_rotate2normal()
	#FBI_Translate()
	#FBI_Fine_adjust()
	
	
	Export2PDB('FBI-loop-short')
	return 0

if __name__ == '__main__':
	main()

