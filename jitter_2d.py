import os
import random
from random import uniform
import math
from numpy import *
from numpy.random import standard_normal
from numpy.random import *


def xyz_morph(xyzFile,size):
	# perturbed the geometry based on size
	# this function calls 2 other functions:
	# 1. xyz_perturber(x,y,z,size)
	# 2. too_close_checker(morphed_xyzFile)

	inFile = open(xyzFile,'r')
	temp_lines = inFile.readlines()
	
	size = float(size)

	# look to see if there is a file freeze.xyz telling which coords not to perturb
	freeze = os.path.isfile("freeze.xyz")
	
	# passes indicates whether the perturbed geometry is accepted or if we should continue looping to try a new set of coords
	passes = 1
	while passes == 1:
		outstring = xyzFile[:-4] + '_morphed.xyz'
		outFile = open(outstring,'w')
		new_lines = []
		for indx,line in enumerate(temp_lines):
			new_line = line.lstrip()
			new_lines = new_line.split()
			if len(new_lines) > 0:
				#calculate the perturbation vector, convert it to Cartesians, add it to the old coordinates
				# if input file has F, don't move, if it hase M, move (F for freeze, M for move)
				move = new_lines[0]
				if new_lines[0] == "F":
					print "coordinate %d is frozen" % indx
					new_x = float(new_lines[2])
					new_y = float(new_lines[3])
					new_z = float(new_lines[4])
				
				elif new_lines[0] == "S":
					print "coordinate %d is to be moved" % indx
					theta = vonmises(math.pi,0)
					phi = (vonmises(math.pi,0)/2)
					# r is divided by 4 to move these coordinates less, S is for small move
				#	r = standard_normal()*size/4
					r = size/2
					new_x = float(new_lines[2]) + r*math.cos(phi)*math.sin(theta)
					new_y = float(new_lines[3]) + r*math.sin(theta)*math.sin(phi)
					new_z = float(new_lines[4]) + r*math.cos(theta)/4
					#new_z = float(new_lines[4])
				else:
					print "coordinate %d is to be moved" % indx
					theta = vonmises(math.pi,0)
					phi = (vonmises(math.pi,0)/2)
					#r = standard_normal()*size
					r = size
					new_x = float(new_lines[2]) + r*math.cos(phi)*math.sin(theta)
					new_y = float(new_lines[3]) + r*math.sin(theta)*math.sin(phi)
					new_z = float(new_lines[4]) + r*math.cos(theta)/2
			#		new_z = float(new_lines[4])
					#write the new coords to the xyzFIle_morphed.xyz
				outFile.write(new_lines[1] + "   " + str(new_x) + "   " + str(new_y) + "   " + str(new_z) + "\n")
		outFile.close()
			#check to see if any atoms are too close, the function will update passes accordingly and we will continue or exit the while loop based on passes
		passes = too_close_checker(outstring)
	inFile.close()
	#return the new set of coords to the main program loop
	return outstring


def too_close_checker(xyz_file):
	too_close = 1.2
	inFile = open(xyz_file,'r')
	temp_lines = inFile.readlines()

	its2close_number = 0
	for line in temp_lines:
		new_line = line.lstrip()
		new_line = new_line.split()
		if len(new_line) > 0:
			x1 = new_line[1]
			y1 = new_line[2]
			z1 = new_line[3]
		for line2 in temp_lines:
			new_line2 = line2.lstrip()
			new_line2 = new_line2.split()
			if len(new_line2) > 0:
				x2 = new_line2[1]
				y2 = new_line2[2]
				z2 = new_line2[3]
				distance = math.sqrt(math.pow(float(x2)-float(x1),2) + math.pow(float(y2)-float(y1),2) + math.pow(float(z2)-float(z1),2))
				#check if it's too close
				if distance < too_close:
					#print "distance is " + str(distance)
					its2close_number += 1
		if its2close_number > 1:
			print "bad geom"
			its2close_number = 0
			bad_struct = 1
			break
		else:
			its2close_number = 0
			bad_struct = 0
		

	inFile.close()
	return bad_struct





# ask user which .xyz file to import
xyz_file = raw_input("Enter coordinate file (Ex. geom.xyz): ")


#ask user how many jitter sizes to use, how many jitters per size, then ask each of the jitter size
N_jitter_sizes = int(raw_input("Enter number of jitter sizes: "))
jitters_per_size = int(raw_input("Enter number of jitters per size: "))
jitter_sizes = []
i = 0
while i < N_jitter_sizes:
	j_size = float(raw_input("Enter jitter size %d: " % i))
	jitter_sizes.append(j_size)
	i += 1
	

for x in range(len(jitter_sizes)):
	print "x is %d" % x
	jitter_size = jitter_sizes[x]
	for y in range(jitters_per_size):
	#	jittered_file = xyz_morph(xyzfile,jitter_size)
		new_xyzFile = xyz_morph(xyz_file,jitter_size)
		print new_xyzFile
		# rename the new_xyzFile accordingly
		os.rename(new_xyzFile,xyz_file[:-4] + "_%d_%d.xyz"%(x,y))
				


