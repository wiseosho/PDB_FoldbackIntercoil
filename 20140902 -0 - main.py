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
PII = math.atan(1)*4/180
Prmzeros = [720.0*PII/21.0, 1.0*PII, 0.0*PII, 0.2]
Risezero = 3.4


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
	
def Gen3dRot(prm, R, k):	# prm:BPprm의 subarray/ R: rotationmat/ k: kth sequence.
	
	#rho, y axis
	Ry = Rotaxis(np.array([0,1,0]), prm[k][1])
	#omega, z axis
	Rz = Rotaxis(np.array([0,0,1]), prm[k][0])
	#tau, x axis
	Rx = Rotaxis(np.array([1,0,0]), prm[k][2])
	R = np.dot(R, np.dot (Rz, Ry, Rx))
	
	return R
                  
def Initialization_Coordinate(Duplex, DuplexN, Duplen, Mayday): #Prmtr : Parameters
# Sequence Initialization/ Duplex : Sequence of Characters.
	Duplex = []
	Duplen = []
	Duplex.append([])
	Duplex[0] = list('ACTG') ; Duplen.append(len(Duplex[0]))
	Duplex.append([])
	Duplex[1] = list('GTCA') ; Duplen.append(len(Duplex[1]))
	print Duplen
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
	
	return Duplex, DuplexN, Duplen, Mayday

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
	print Prmtable
	for j in range (len(Duplen)):
		BPprm.append([])
		for i in range (Duplen[0] - 3) :
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
				Type = join(duplex_c)
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

	return BPprm

def Construct_Coordinate(Duplex, Duplen, mayday, BPprm):
	R = np.eye(3,3)
	sum_R = array([0,0,0])
	sum_D = array([0,0,0])
	
	for j in range (len(Duplen)) :
		for i in range (1, Duplen[j]) :
			Gen3dRot(BPprm[j], R, i-1)
			
			sum_R = sum_R + np.column_stack(R)[2]
			
			sum_D = sum_D + BPprm[j][i-1][3]*np.column_stack(R)[1]
			
			for k in range(45) :
				mayday[j][i][k] = np.dot(R, mayday[j][i][k]) + Risezero*(sum_R) + sum_D
				
def Load_Atomname():
	At = open('./Data-ipython/atoms2.txt', 'r')
	atomname = [[],[],[],[]]
	i = 0 ; 
	for line in At.readlines():
		i = i+1
		atomname[int((i-1)/41)].append([list(line.strip('\n\r\t'))])
	return atomname

def Export2PDB(Duplex, DuplexN, Duplen, mayday, Filename):
	initatom = Load_Atomname()
	
	print initatom
	

def main():
	print "hello"
	Duplex = [];	DuplexN = [];	Mayday = []; BPprm = [] ; Duplen = []
	PrmTurnOn = 'off'
	Duplex, DuplexN, Duplen, Mayday = copy.deepcopy(Initialization_Coordinate(Duplex, DuplexN, Duplen, Mayday))
	BPprm = copy.deepcopy(Initialization_BPparameters(Duplex, Duplen, BPprm, Prmzeros, PrmTurnOn))
	mayday = copy.deepcopy(np.array(Mayday))
	print mayday
	Construct_Coordinate(Duplex, Duplen, mayday, BPprm)
	print '====================================='
	print mayday
	Export2PDB(Duplex, DuplexN, Duplen, mayday, 'test')
	
	return 0

if __name__ == '__main__':
	main()

