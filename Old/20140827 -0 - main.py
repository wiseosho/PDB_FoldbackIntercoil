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
import copy

def Initialization_Coordinate(Duplex, DuplexN, Mayday): #Prmtr : Parameters
# Sequence Initialization/ Duplex : Sequence of Characters.
	Duplex = []
	Duplex.append([])
	Duplex[0] = list('ACTG')
	Duplex.append([])
	Duplex[1] = list('GTCA')
	# Sequence convert from Character into Numbers
	DuplexN = []#[Sequence][Coordinate]

	for y in range (len(Duplex)):
		DuplexN.append([])
		for x in range (len(Duplex[y])):
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

#===================Load interbasepair parameters.=====================
	Prmtable = []
	casenum = 0
	duplex_c = list('RABC')
	duplex_x = list('RY')

	omega = open('./Data-ipython/domega.txt', 'r')
	rho = open('./Data-ipython/drho.txt', 'r')
	tau = open('./Data-ipython/dtau.txt', 'r')
	slide = open('./Data-ipython/dslide.txt','r')
	Prmtable.append([])
	for line in omega.readlines():
		Prmtable[len(Prmtable)-1].append(float(line.strip('\n')))
	Prmtable.append([])
	for line in rho.readlines():
		Prmtable[len(Prmtable)-1].append(float(line.strip('\n')))
	Prmtable.append([])
	for line in tau.readlines():
		Prmtable[len(Prmtable)-1].append(float(line.strip('\n')))
	Prmtable.append([])
	for line in slide.readlines():
		Prmtable[len(Prmtable)-1].append(float(line.strip('\n')))
	#print Prmtable
	
	return Duplex, DuplexN, Mayday



def main():
	print "hello"
	Duplex = [];	DuplexN = [];	Mayday = []; Parameters = []
	Duplex, DuplexN, Mayday = copy.deepcopy(Initialization_Coordinate(Duplex, DuplexN, Mayday))
	print Parameters

# Test
	return 0

if __name__ == '__main__':
	main()

