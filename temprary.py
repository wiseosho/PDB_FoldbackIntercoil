#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  temprary.py
#  
#  Copyright 2014 Junyoung Son <june@june-VirtualBox>
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

TailBodyLoop = [7,13,15,24,30]
# [0,0,20,24,37,
# 0,1,2,3,4,5,6,7 | 8,9,10,11,12,13 | 14,15,16,17,18 , 19,20,21,22,23,24 | 25,26,27,28,29,30 | 31,32,33,34,35,36,37,38
# 0,1,2,3,4,5,6,7 | 8,9,10,11,12,13 , 14,15,16,17,18 | 19,20,21,22,23,24

def Construct_FBI_Loop1(TailBodyLoop):
	global mayday, Duplex
	
	End = len(mayday[0])
	Helix = ['Not Parallel', 'Far']
	Inner_K = 0 ; Canon = 4 ; unit = [1,1] ; Accu = .5 
	# Unit : Radius Factor for [Inner K, Outer K]
	Outer_K = 0 ; Start = 1; #End = 26
	parallel = [] ; dist = [10] ; para = []
	#while (dist[len(dist)-1] > 1) :
	
	####################Rotate Inner Circled Loop till Close.
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
			print 'distance' + str(distance)
			if (abs(dist[(len(dist))-1] - Canon) < Accu):
				dist = []
				Helix[1] = 'Close'
			elif abs(dist[(len(dist))-1] < Canon) :
			#elif Compare(dist, 'dist') == 'turning' : 
				mayday[0] = Construct_FBI(mayday[0], Duplex[0], 0, 1, Start, End, 0 , -unit[0], Order)
				dist = []
				unit[0] = unit[0] * 0.1
				print 'Reduced'
		Order = Order +1
		Export2PDB(Duplex, DuplexN, Duplen, mayday,Weave, 'FBI-loop-Long'+str(Order))
		
		####################Rotate Outer Circled Loop till parallel
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
def Construct_FBI_Loop2():
	

def main():
	
	return 0

if __name__ == '__main__':
	main()

