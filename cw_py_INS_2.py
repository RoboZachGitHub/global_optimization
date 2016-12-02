from random import randint
from math import sqrt
import sys

class Atom:
	def __init__(self):
		self.x_index = 0
		self.y_index = 0
		self.z_index = 0
		self.x_cart = 0.0
		self.y_cart = 0.0
		self.z_cart = 0.0
		self.element = "B"

        def set_element(self, elem):
                self.element = elem

        def set_cartesians(self, grid_params):  
                self.x_cart = float(self.x_index*grid_params.major_mesh/grid_params.minor_mesh)
                self.y_cart = float(self.y_index*grid_params.major_mesh/grid_params.minor_mesh)
                self.z_cart = float(self.z_index*grid_params.major_mesh/grid_params.minor_mesh)

	


class Grid:
	def __init__(self):
		# grid_x etc. is the number of major grid points in arbitrary units
		# max_step_x etc is the number of minor grid points that can be jumped 
		self.grid_x = 6
		self.grid_y = 6
		self.grid_z = 1
		self.max_step_x = 3 
		self.max_step_y = 3
		self.max_step_z = 0
		# major_mesh is the distance between two major grid points
		# major_mesh is in units of Angstroms
		# imagine graph paper with bolded lines and lighter lines, major_mesh is the distance between bold lines
		# minor_mesh is the number of small blocks between each major line
		# minor_mesh has arbitrary units
		# (i.e. if minor_mesh is 2, there is one extra point between the major lines
		self.major_mesh = 1.6
		self.minor_mesh = 2
		self.num_x_points = (self.grid_x-1)*(self.minor_mesh) + 1
		self.num_y_points = (self.grid_y-1)*(self.minor_mesh) + 1
		self.num_z_points = (self.grid_z-1)*(self.minor_mesh) + 1
		# too_close is the constraint that no two atoms can bee too close to each other
		# too_close is in units of angstroms
		self.too_close = 1.0
		self.num_atoms = 3
		self.formula = ["B", "B", "C"]



def dist(a, b, grid):
	# a and b are Atom objects
	# returns the cartesian difference in Angstroms between two points on the grid (i.e. two atoms)
	# used for checking if two atoms are too close to each other
	dx = a.x_cart - b.x_cart
	dy = a.y_cart - b.y_cart
	dz = a.z_cart - b.z_cart
	d2 = dx*dx + dy*dy + dz*dz
	#print "distance is: %f" % sqrt(d2)
	#if sqrt(d2) <= grid.too_close:
		#print "fail. sqrt(d2)/distance is: %f" % sqrt(d2)
		#print "    grid.too_close is: %f" % grid.too_close
	return sqrt(d2)

def bound_check(grid, a):
	# determines if an atom to be placed is out of the grid boundaries
	# returns True if within boundaries or False if outside boudaties
	check_x = a.x_index >= 0 and a.x_index < grid.num_x_points
	check_y = a.y_index >= 0 and a.y_index < grid.num_y_points
	check_z = a.z_index >= 0 and a.z_index < grid.num_z_points
	return check_x and check_y and check_z	

def rand_atom(grid):
	# places and atom randomly on the grid and returns the atom with proper cartesians set
	atom = Atom()
	atom.x_index = randint(0, grid.num_x_points-1) 
	atom.y_index = randint(0, grid.num_y_points-1)
	atom.z_index = randint(0, grid.num_z_points-1)
	atom.set_cartesians(grid)
	# randomly select an atom from the formula list, then delete it from the formula list
	#print "formula is: " + str(grid.formula)
	#rand_atom = grid.formula.pop(randint(0,len(grid.formula)))
	#atom.element = rand_atom
	#print "rand_atom is: " + str(rand_atom)
	#print "formula is now: " + str(grid.formula)
	gum_machine(grid,atom)
	
	return atom

def take_step(grid, start):
	# takes a random step constrained by grid.max_step_x etc. from current atom position
	# returns new Atom object
	atom = Atom()

	# set constraints based on max_steps and bounds of grid
	start_x = start.x_index
	start_y = start.y_index
	start_z = start.z_index
	#print "start_z: " + str(start_z)
	temp_max_x = grid.max_step_x
	temp_max_y = grid.max_step_y
	temp_max_z = grid.max_step_z
	#print "temp_max_z: " + str(temp_max_z)
	lower_x = start_x - temp_max_x
	upper_x = start_x + temp_max_x
	lower_y = start_y - temp_max_y
	upper_y = start_y + temp_max_y
	lower_z = start_z - temp_max_z
	upper_z = start_z + temp_max_z
	#print " num_xyz_points:  %d  %d  %d " % (grid.num_x_points,grid.num_y_points,grid.num_z_points)
	if lower_x < 0:
		lower_x = 0
	if upper_x > grid.num_x_points - 1:
		upper_x = grid.num_x_points - 1
	if lower_y < 0:
		lower_y = 0
	if upper_y > grid.num_y_points - 1:
		upper_y = grid.num_y_points - 1
	if lower_z < 0:
		lower_z = 0
	if upper_z > grid.num_z_points - 1:
		upper_z = grid.num_z_points - 1

	# calculate the new positions
	new_x = randint(lower_x, upper_x)
	new_y = randint(lower_y, upper_y)
	#print "(lower_x, upper_x): (%d, %d)" % (lower_x, upper_x)
	#print "(lower_y, upper_y): (%d, %d)" % (lower_y, upper_y)
	new_z = randint(lower_z, upper_z)

	atom.x_index = new_x
	atom.y_index = new_y
	atom.z_index = new_z
	atom.set_cartesians(grid)

#	print "(lower_x, upper_x):  (%d, %d)" % (lower_x, upper_x)
#	print "(lower_y, upper_y):  (%d, %d)" % (lower_y, upper_y)
#	print "(lower_z, upper_z):  (%d, %d)" % (lower_z, upper_z)

	return atom

def gum_machine(grid, a):
	# randomly select an atom from the formula list, then delete it from the formula list
	#print "formula is: " + str(grid.formula)
	rand_atom = grid.formula.pop(randint(0,len(grid.formula)-1))
	a.element = rand_atom
	#print "rand_atom is: " + str(rand_atom)
	#print "formula is now: " + str(grid.formula)
	

def read_INS():
    inFile = open("INS_cw_py.txt",'r')
    temp_lines = inFile.readlines()
    l = 0
    for line in temp_lines:
        useful_line = temp_lines[l]
        new_lines = useful_line.split()
        if l == 1:
#           print "setting grid_x as: " + str(new_lines[2])
            grid_x = int(new_lines[2])
            l = l + 1
        elif l == 2:
            grid_y = int(new_lines[2])
#           print "setting grid_y as: " + str(new_lines[2])
            l = l + 1
        elif l == 3:
            grid_z = int(new_lines[2])
#           print "setting grid_z as: " + str(new_lines[2])
            l = l + 1
        elif l == 4:
            major_mesh = float(new_lines[2])
#           print "setting major_mesh  as: " + str(new_lines[2])
            l = l + 1
        elif l == 5:
#           print "setting mesh as: " + str(new_lines[1])
            mesh = float(new_lines[1])
            l = l + 1
        elif l == 6:
            max_step_x = int(new_lines[3])
            l = l + 1
        elif l == 7:
            max_step_y = int(new_lines[3])
            l = l + 1
        elif l == 8:
            max_step_z = int(new_lines[3])
            l = l + 1
        elif l == 9:
            too_close = float(new_lines[2])
            l = l + 1
	elif l == 10:
	    population = int(new_lines[1])
	    l = l + 1
	elif l == 13:
	    chem_formula = []
	    atom_amounts = []
	    n_atom_types = 0
	    n_atoms = 0
	    new_lines = temp_lines[l].split()
	    k = 1
	    while k < len(new_lines):
	        if k%2 != 1:
	            atom_amounts.append(int(new_lines[k]))
		    n_atoms += int(new_lines[k])
		if k%2 == 1:
		    n_atom_types += 1
		    for x in range(0,int(new_lines[k+1])):
		        chem_formula.append(new_lines[k])
		k = k + 1
	    l = l + 1
        else:
            l = l + 1

 
    parameters = [grid_x ,grid_y, grid_z, major_mesh, mesh, max_step_x, max_step_y, max_step_z, too_close, n_atoms, chem_formula, population]
    return parameters

def create_xyz_file(atom_list,tag):
	outFile = open(str(tag)+".xyz",'w')
	outString = ''
	for a in atom_list:
	    outString += "%s  %f  %f  %f\n" %(a.element,a.x_cart,a.y_cart,a.z_cart)  
	outFile.write(outString)
	outFile.close()

def main():
	#print "Started"
	# read info from INS_cw_py.txt, create grid, and set parameters
	grid = Grid()
	grid_params = read_INS()	
	population = grid_params[11]
	grid.grid_x = grid_params[0]
	grid.grid_y = grid_params[1]
	grid.grid_z = grid_params[2]
	grid.major_mesh = grid_params[3]
	grid.minor_mesh = grid_params[4]
	initial_max_x = grid_params[5]	
	initial_max_y = grid_params[6]
	initial_max_z = grid_params[7]
	grid.max_step_x = initial_max_x
	grid.max_step_y = initial_max_y
	grid.max_step_z = initial_max_z
	grid.too_close = grid_params[8]
	grid.num_atoms = grid_params[9]
	initial_formula = grid_params[10]
	#print "initial_formula is: " + str(initial_formula)
	# make a copy of a list, grid.formula = initial_formula just makes a new pointer
	grid.formula = list(initial_formula)
	grid.num_x_points = (grid.grid_x-1)*grid.minor_mesh + 1
	grid.num_y_points = (grid.grid_y-1)*grid.minor_mesh + 1
	grid.num_z_points = (grid.grid_z-1)*grid.minor_mesh + 1
	
	max_fails = 10000
	exit_fails = 1000000

	for n in range(0,population):
		#print "initial_formula remains: " + str(initial_formula)
		grid.formula = list(initial_formula)
		#print "grid.formula is: " + str(grid.formula)
		#start with a single random atom appended to list of atoms
		atoms = [rand_atom(grid)]
		#print "starting with %f %f %f" % (atoms[0].x_cart,atoms[0].y_cart,atoms[0].z_cart)	
	
		#take the walk
		fails = 0
		total_fails = 0
		while len(atoms) < grid.num_atoms:
			# deal with corner and surrounded atom cases
			#print "fails is: " + str(fails)
			if fails < max_fails and total_fails <= exit_fails:		
			    next_atom = take_step(grid, atoms[-1])
			elif fails < exit_fails and total_fails <= exit_fails:
			    # atom has walkied into corner or surrounded itself or there are no more legal moves
			    # increase step size dynamically and continue trying
			    print "fails is greater than max_fails less than exit_fails.   fails: %d   exit_fails: %d" % (fails, exit_fails)
			    print "\n\n total_fails: %d    exit_fails: %d  \n\n" %(total_fails, exit_fails)
			    if grid.max_step_x < grid.num_x_points:
				grid.max_step_x += 1
			    if grid.max_step_y < grid.num_y_points:
				grid.max_step_y += 1
			    if grid.max_step_z < grid.num_z_points:
				grid.max_step_z += 1
			    next_atom = take_step(grid, atoms[-1])
			    print "max steps are now: (max_step_x, max_step_y, max_step_z)  (%d, %d, %d)" %(grid.max_step_x, grid.max_step_y, grid.max_step_z)	
			   
			    fails = 0
			
			elif total_fails >= exit_fails:
			    # there have been many fails, should consider new grid options
			    print "error! grid parameters may be too restrictive for a large sample of legal walks!"
			    sys.exit()
	
			# only append the new atom and assign it an element if it passes some constraint checks
		
			# first check if out of bounds (should rarely if ever fail this check)
			if not bound_check(grid, next_atom):
#				print "    * Failed boundary check"
				fails += 1
				total_fails += 1
				continue

			# now check if new_atom is too close to previously placed atoms
			if any([dist(atom, next_atom, grid) < grid.too_close for atom in atoms]):
#				print "    * Failed closeness check"
				fails += 1
				total_fails += 1
				#print "total_fails is: " + str(total_fails)
				continue
					
#			print "(fails: %d   total_fails: %d)" % (fails, total_fails)
			print "  Adding atom %d" % (len(atoms) + 1)
			gum_machine(grid,next_atom)
			atoms.append(next_atom)
			# reset max step sizes to their initial values and reset fail count
			grid.max_step_x = initial_max_x
			grid.max_step_y = initial_max_y
			grid.max_step_z = initial_max_z
			fails = 0
			total_fails = 0
		# write trial structure to text file
		print "creating file"
		create_xyz_file(atoms,n)


main()






	
#grid = Grid()
#for x in range(0,100000):
#	a = rand_atom(grid)
#	b = take_step(grid,a)		



