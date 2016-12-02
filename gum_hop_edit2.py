#standard modules
import sys
import os
import os.path
import fnmatch
import time
import shutil
import subprocess
from random import seed
from random import uniform
from random import randrange
import math
##mpi modules
from mpi4py import MPI
##numpy modules, for random perturbations in spherical coords and such
from numpy import *
from numpy.random import standard_normal
from numpy.random import *



###################################
######## CLASS DEFINITION #########
###################################

class trial_geom:

        trial_generation = 0
        trial_id = 0
        trial_energy = 0.0
	trial_geom_type = 0
        trial_gjf = "initiate.gjf"
        trial_log = "initiate.log"

        def __call__(self):
                print "%s object has been called." % self.trial_gjf

        def set_genID(self,gen,ID):
                self.trial_generation = gen
                self.trial_id = ID

	def set_geom_type(self,geom):
		self.trial_geom_type = geom

        def set_energy(self,E):
                self.trial_energy = E

        def set_names(self,name):
                self.trial_gjf = name[:-4] + '.gjf'
                self.trial_log = name[:-4] + '.log'

        def print_info(self):
                print "[generation, proc ID, name, energy]: (%d,%d,%s,%f)" % (self.trial_generation,self.trial_id,self.trial_gjf[:-4],self.trial_energy)

	def info_return(self):
		info_string = "[generation, proc ID, name, energy]: (%d,%d,%s,%f)" % (self.trial_generation,self.trial_id,self.trial_log,self.trial_energy)
		return info_string

###################################
######  FUNCTIONS  ################
###################################  

def run_gjf(gjf_name):

        file2run = gjf_name

        #copy the gjf to gjf_log_copies with appropriate indexing
        x = 0
        gjf2copy = file2run
        gjf2copy_string = file2run[:-4] + '_%d.gjf' % x
        gjf2copy_path = copies_path + '/' + gjf2copy_string
        while os.path.isfile(gjf2copy_path):
                x = x + 1
                gjf2copy_string = gjf2copy_string[:gjf2copy_string.find('_')] + '_' + str(x) + '.gjf'
                gjf2copy_path = copies_path + '/' + gjf2copy_string
        shutil.copy2(gjf2copy,gjf2copy_path)

        #run the gjf
        print "running structure " + file2run
        commands = "g09 " + file2run
        subprocess.Popen(commands, shell=True).wait()
        print "structure " + file2run + " done."

        #copy the .log file
        log2copy = file2run[:-4] + '.log'
        log2copy_path = copies_path + '/' + log2copy
        shutil.copy2(log2copy,log2copy_path)

        #check the log file to see if calculation ran fine
        check = check_log(log2copy)

        if check == 0:
                #calculation ran fine
                pass
        elif check == 1:
                resub(log2copy)
        elif check == 2:
                try_again(log2copy)


        return log2copy


def check_log(logfile):
        ##checks a Gaussian *.log file and returns a different value if it completed fine, has an error, or simply didn't complete
        inFile = open(logfile,'r')

        okay = -1
        text = inFile.read()
        search = 'Normal termination of Gaussian'
        index = text.find(search)
        okay = index
        if okay != -1:
                return 0
        else:
                search = 'termination request processed by link 9999'
                index = text.find(search)
                okay = index
                if okay != -1:
                        return 1
                else:
                        print 'gaussian error during optimization'
                        return 2




def resub(log_file):
        print "in resub"
        #grabs current geometry from a log file with link9999 error and resubmits a new job to complete the optimization
        xyz_file = get_xyz(log_file)
        gjf_name = make_gjf(xyz_file,0)
        run_gjf(gjf_name)

def try_again(log_file):
        print "in try_again"
        #submits an entirely new random seed in replacement of a calculation with error
	box_size = int(num_atoms)/2
        tmp_name = random_seed(box_size,box_size,box_size)
        gum_machine(tmp_name)
        new_name = log_file[:-4] + ".xyz"
        if tmp_name != new_name:
                shutil.copy2(tmp_name,new_name)
        gjf_name = make_gjf(new_name,0)
        run_gjf(gjf_name)

def random_seed(box_x,box_y,box_z):
        ## need to replace the X's with the normal atom list order
        its2close = 2.3
        its2close_check = 0
        seed()
        outFileName = "gen%d_id%d.xyz" % (generation,my_id)
        outFile = open(outFileName,'w')
        lines = []
        writeit_check = 0
        i = 0
        while i < num_atoms:
                #print "i is: " + str(i)
                x_coord = uniform(0.0,box_x)
                y_coord = uniform(0.0,box_y)
                z_coord = uniform(0.0,box_z)
                #X is used as a placeholder, gum_machine() will replace X with an element type
                if len(lines) > 0:
                        for line in lines:
                               # "print line is " + str(line)
                               # "print line index is " + str(lines.index(line))
                                set_coord = line.split()
                                x1 = float(set_coord[1])
                                y1 = float(set_coord[2])
                                z1 = float(set_coord[3])
                                distance = math.sqrt(math.pow(float(x1)-float(x_coord),2) + math.pow(float(y1)-float(y_coord),2) + math.pow(float(z1)-float(z_coord),2))
                                #print "distance is : " + str(distance) + " too close is: " + str(its2close)
                                if distance < its2close:
                                #       print "distance was too close: " + str(distance)
                                        its2close_check = its2close_check + 1
                                else:
                                        continue
                                #       print "distance " + str(i) + "-" + str(lines.index(line)) +" was fine: " + str(distance)

                        if its2close_check > 0:
                        #       print "bad geometry, its2close_check is :" + str(its2close_check)
                                its2close_check = 0
                                writeit_check = 1

                        else:
                                writeit_check = 0
                                        #coord_string = "X    %f     %f     %f  \n" % (x_coord,y_coord,z_coord)
                                        #lines.append(coord_string)
                                        #i = i + 1

                #elif len(lines) == 0:
                #       print "len was 0"
                #       coord_string = "X    %f    %f    %f  \n" % (x_coord,y_coord,z_coord)
                #       lines.append(coord_string)
                #       i = i + 1       
                #else:
                #       print "there was an error"
                if writeit_check == 0:
                #       print "setting coordinate " + str(i)
                        coord_string = " X    %f    %f    %f \n" % (x_coord,y_coord,z_coord)
                        lines.append(coord_string)
                        i = i + 1
                elif writeit_check == 1:
                #       print "not setting coordinate, failed distance check"
                        i = i

        outFile.writelines(lines)

        return outFileName




def make_gjf(xyzfile,mod):
        #mod can be 0,1. mod 1 is for the altRoute
        inFile = open("INS_gum.txt",'r')
        temp_lines = inFile.readlines()
        l = 0
        for line in temp_lines:
                if l == 10:
                        useful_line = temp_lines[l]
                        Route = useful_line.lstrip('Route:')
                        Route = Route.lstrip()
                        l = l + 1
                if l == 11:
                        useful_line = temp_lines[l]
                        altRoute = useful_line.lstrip('altRoute:')
                        altRoute = altRoute.lstrip()
                        l = l + 1
                if l == 12:
                        useful_line = temp_lines[l]
                        new_lines = useful_line.split()
                        mem = new_lines[1]
                        l = l + 1
                if l == 13:
                        useful_line = temp_lines[l]
                        new_lines = useful_line.split()
                        nproc = new_lines[1]
                        l = l + 1
                else:
                        l = l + 1

        inFile.close()
        inFile = open(xyzfile,'r')
        gjf_name = xyzfile[:-4] + '.gjf'
        outFile = open(gjf_name,'w')

        lines = []
        lines.append('%mem=' + mem + '\n')
        lines.append('%nproc=' + nproc + '\n')
        if mod == 0:
                lines.append('# ' + Route + '\n')
        elif mod == 1:
                lines.append('# ' + altRoute + '\n')
        lines.append('structure ' + xyzfile + '\n\n')
        lines.append(charge + ' ' + mult + '\n')
        #while 1 means infinite loop until break
        #cycle through the coords in .xyz and add them to lines until reaching end of .xyz
        while 1:
                line = inFile.readline()
                if not line:
                        break
                lines.append('' + str(line))
        lines.append('\n')

        #write it, bitch!!!!
        outFile.writelines(lines)

        #close up them files, Son!!!
        inFile.close()
        outFile.close()

        return  gjf_name



def gum_machine(xyzfile):
        inFile = open(xyzfile,'r')
        tmp_out_name = xyzfile[:xyzfile.find('.xyz')] + '_tmp.xyz'
        outFile = open(tmp_out_name,'w')
        #create a list holding integers 1 through number of atoms
        list1 = []
        rand_list = []
        for i in range(0,num_atoms):
                list1.append(i)
        #now randomize them in a new list
        x = num_atoms
        while(x>0):
                my_rand = randrange(0,x)
                rand_list.append(list1[my_rand])
                list1.remove(list1[my_rand])
                x = x - 1
        #now use the globally known list of atom types with the random integer list to selects the elements from the element list in random order
        i = 0
        lines = []
        temp_lines = inFile.readlines()
        elem_list = []
        for line in temp_lines:
#       while i < num_atoms:
               	elem = atom_types[rand_list[i]]
               	elem_list.append(elem)
#              	print "in gum_machine elem is " + str(elem)
        	useful_line = temp_lines[i]
               	new_lines = useful_line.split()
               	lines.append(' ' + str(elem) + '    ' + new_lines[1] + '   ' + new_lines[2] + '   ' + new_lines[3] + '\n')
#              	print "in gum_machine lines is " + str(lines)
               	i = i + 1
#        lines.append('\n')

        outFile.writelines(lines)
        inFile.close()
        os.remove(xyzfile)
        shutil.move(tmp_out_name,xyzfile)

        return elem_list




def order_checker(ordered_list,list_of_ordered_lists):
        #checks to make sure the gum_machine doesn't make copies of the same structure


        #### check for similarity to others in list
        different_check = 0
        for i in range(0,len(list_of_ordered_lists)):
               # "i is: " + str(i)
                object_i = list_of_ordered_lists[i]
                #print "object_i is " + str(object_i)
                for j in range(0,num_atoms):
                        #print "j is: " + str(j)
                        #print "ordered_list[j] is " + str(ordered_list[j]) + " object_i[j] is " + str(object_i[j])
                        if ordered_list[j] != object_i[j]:
                        #       print "different. break"
                                different_check = different_check + 1
                                break
                        elif ordered_list[j] == object_i[j]:
                                pass
                        #       print ordered_list[j] + " is equal to " +  object_i[j]

        if different_check == len(list_of_ordered_lists):
                #print "the ordering in unique"
        #       print "ordering: " + str(ordered_list)
        #       print "others: " + str(list_of_ordered_lists)
                return "fine order"

        elif different_check < len(list_of_ordered_lists):
                #print "the ordering is NOT unique"
        #       print "ordering: " + str(ordered_list)
        #       print "others: " + str(list_of_ordered_lists)
		pass

        if len(list_of_ordered_lists) >= num_distinguishables:
                print "no unique orderings left, should create a new geometry"
                return "new geometry"

        if len(list_of_ordered_lists) < num_distinguishables and different_check < len(list_of_ordered_lists):
                print "there are still unique orderings left. should use gum_machine again."
                return "gum machine again"




def get_energy(log_file):
        filename = log_file
        inFile = open(filename,'r')
        outPath = filename[:filename.find('.log')] + '_energy_temp'
        outFile = open(outPath,'w')

        lines = []
        text = inFile.read()
        lines.append(text[text.rfind('SCF Done'):])

        outFile.writelines(lines)
        outFile.close()
        inFile.close()

        inFile2 = open(outPath, 'r')
        outPath2 = filename[:filename.find('.log')] + '_energy'
        outFile2 = open(outPath2,'w')

        lines = []
        text = inFile2.read()
        lines.append(text[:text.find('cycles')])

        outFile2.writelines(lines)
        outFile2.close()
        inFile2.close()

        inFile3 = open(outPath2,'r')
        raw_lines = inFile3.readlines()
        new_lines = []

        # loop through raw_lines and strip everything except what we want on each line
        for line in raw_lines:
                new_line = line.lstrip()
                new_line = new_line.split()
                if len(new_line) > 0:
                        new_lines.append(new_line[4])

        write_string = ''

        #build the string
        for line in new_lines:
                write_string += line

        energy = float(write_string)
        #remove extraneous files
        os.remove(outPath)
        os.remove(outPath2)

        return energy


def get_xyz(log_file):
        #print "log_file is: " + log_file
        inFile = open(log_file,'r')
        outPath = log_file[:-4] + '_tmp.xyz'
        outFile = open(outPath,'w')

        lines = []
        text = inFile.read()
        lines.append(text[text.rfind('Standard orientation'):text.rfind('Rotational constants')])

        outFile.writelines(lines)
        outFile.close()
        inFile.close()

        #now edit this file to get rid of unwanted text
        inFile2 = open(outPath,'r')
        outPath2 = outPath[:outPath.find('_tmp.xyz')] + '_out.xyz'
        #print "outPath2 is " + str(outPath2)
        outFile2 = open(outPath2,'w')

        raw_lines = inFile2.readlines()
        temp_lines = []
        new_lines = []
        i = 0
        #loop through raw lines and remove first five lines and last line, then store in temp lines
        for line in raw_lines:
                if i > 4:
                        if line.find('--') == -1:
                                temp_lines.append(line)
                i = i + 1
        #loop through temp lines and strip everything except what we want on each line (i.e. just the element and cartesians)
        for line in temp_lines:
                new_line = line.lstrip()
                new_line = new_line.split()
                if len(new_line) > 0:
                        new_lines.append("%s    %s    %s    %s  " % (new_line[1],new_line[3],new_line[4],new_line[5]))

        inFile2.close()
        #build the string to write to the _out.xyz file
        write_string = ''
        for line in new_lines:
                write_string += line + '\n'

        #write it, dude
        #print "write_string is : " + write_string
        #print "outFile2 is : " + str(outFile2)
        outFile2.write(write_string)
        outFile2.close

        #remove extraneous files
#       os.remove(outPath)

        #return the name of the _out.xyz file
        return outPath2


def xyz_morph(xyzfile,size,move,too_close):
#       print "inside xyz_morph, my_id: " + str(my_id)
        #perturbs the geometry, size sets the size of the moves, move is used to allow multiple different perturbation options

        #this functions calls the three functions:
        # 1. xyz_perturber(x,y,z,size,box_x,box_y,box_z)
        # 2. too_close_checker(morphed_xyzfile)
        #print "for DEBUG!!! line 431 my_id is: " + str(my_id)

        inFile = open(xyzfile,'r')
        temp_lines = inFile.readlines()

        if size == "std":
                size = 1.0
        elif size != "std":
                size = float(size)
        if move == "std":
                move = 0
        elif move != "std":
                move = int(move)
        if too_close == "std":
                too_close = 1.2
        elif too_close != "std":
                too_close = float(too_close)

        ## passes will indicate whether the perturned geometry is accepted or if we should continue looping to try a new set of coords
        passes = 1
        while passes == 1:
 #              print "DEBUG!!! line 489. in while passes loop, my_id is %d" % my_id
                outPath = xyzfile[:-4] + '_morphed.xyz'
                outFile = open(outPath,'w')
                #print "temp_lines is: " + str(temp_lines)
                new_lines = []

                #loop through temp_lines, grab coords, perturb them, check them, then write
                #grab coords
                for line in temp_lines:
                        new_line = line.lstrip()
                        new_lines = new_line.split()
                        if len(new_lines) > 0:
                                #print "new_lines is: " + str(new_lines)
                                atom_numb = new_lines[0]
                                x_coord = new_lines[1]
                                #print "x_coord is: " + x_coord
                                x = float(x_coord)
                                y_coord = new_lines[2]
                                y = float(y_coord)
                                z_coord = new_lines[3]
                                z = float(z_coord)
                                #print "(atom_numb,x,y,z): (%s,%f,%f,%f)" %(atom_numb,x,y,z)
                        #perturb coord
                        #box_x box_y and box_z are global variables grabbed from INS_gum.txt
                        (a,b,c) = xyz_perturber(x,y,z,size,box_x,box_y,box_z)
                        #print "DEBUG!!! line 530 exited xyz_perturber"
                        #build the string to write to the file
                        write_string = ''
                        a = str(a)
                        b = str(b)
                        c = str(c)
                        write_string = "  %s    %s    %s    %s  \n" % (atom_numb,a,b,c,)
                        #write it, Chump.
   #                    print "write_string is: " + write_string
                        outFile.write(write_string)

                #close the files
                outFile.close()

                #check for atoms too close, if it passes, the function returns something other than 1 for the variable passes
                #print "about to enter too_close_checker"
                passes = too_close_checker(outPath,too_close)
               # print "exited too close checker, passes is: " + str(passes)
                #return the perturbed xyz file
        inFile.close()
        return outPath


def xyz_perturber(x,y,z,size,box_X,box_Y,box_Z):
        #size effects how large the perturbation is
        #absolute value of box_x etc. effects the boundaries of the box in angrtons where atoms can exist

        while 1:
                #convert to spherical coords
                r = math.sqrt(x**2 + y**2 + z**2)
                phi = math.atan2(y,x)
                theta = math.atan2(z,sqrt(x**2+y**2))
                #morph them
                delta_theta = vonmises(math.pi,0)
                new_theta = theta + delta_theta
                delta_phi = vonmises(math.pi,0)/2
                new_phi = phi + delta_phi
                morphed_r = r + (standard_normal(1)/2)*size
                new_r = morphed_r[0]
                #convert back to cartesians
                morphed_x = new_r*outer(cos(new_phi), sin(new_theta))
                new_x = morphed_x[0][0]
                morphed_y = new_r*outer(sin(new_phi), sin(new_theta))
                new_y = morphed_y[0][0]
                new_z = new_r*cos(new_phi)
                #check if within boundary box
                if math.fabs(new_x)<=math.fabs(box_x) and math.fabs(new_y)<=math.fabs(box_y) and math.fabs(new_z)<=math.fabs(box_z):
                        break
                else:
                        pass

        #return the new coordinate
        return(new_x,new_y,new_z)


def too_close_checker(perturbedxyzfile,too_close_distance):

        #print "in too close checker. my_id: " + str(my_id)
        file2open = perturbedxyzfile
        its2close = too_close_distance
        its2close_number = 0
        #initialize bad_struct, bad_struct = 1 if the structure is bad. it will change to 0 if the structure passes this check.
        bad_struct = 1

        inFile=open(file2open,'r')
        temp_lines=inFile.readlines()
        for line in temp_lines:
                new_line = line.lstrip()
                new_line = new_line.split()
                if len(new_line) > 0:
                        x1 = new_line[1]
                        y1 = new_line[2]
                        z1 = new_line[3]
                #       print str(x1) + " " + str(y1) + " " + str(z1)
                for line2 in temp_lines:
                        new_line2 = line2.lstrip()
                        new_line2 = new_line2.split()
                        if len(new_line2) > 0:
                                x2 = new_line2[1]
                                y2 = new_line2[2]
                                z2 = new_line2[3]
                        #       print str(x2)  
                                distance = math.sqrt(math.pow(float(x2)-float(x1),2) + math.pow(float(y2)-float(y1),2) + math.pow(float(z2)-float(z1),2))
                        #       print "distance is: " + str(distance)
                                if distance < its2close:
                        #               print "it's too close. distance is: " + str(distance)
                                        its2close_number = its2close_number + 1
                if its2close_number > 1:
    #                    print "bad geometry, there are atoms too close to this atom"
                        its2close_number = 0
                        bad_struct = 1
                        break
                else:
#                       print "no atoms are too close to this atom"
                        its2close_number = 0
                        bad_struct = 0
        if bad_struct == 1:
                pass
#       print "In too_close_checker(n) its2close_number is:  " + str(its2close_number)
                #print "There are atoms too close together!!!"
        elif bad_struct == 0:
                pass
        #       print "In too_close_checker(n) bad_struct is: " + str(bad_struct)
        #       print "NO ATOMS too close. HOORAH!!!"
        else:
               # print "There was an error in too_close_checker(n). bad_struct != 0 or 1"
                pass
        # the value of bad_struct returned to xyz_morph is used as a conditional statement to create a new set of coords or not
        # bad_struct == 0 is a good structure, bad_struct == 1 means atoms are too close to each other
        return bad_struct





#######################################################################
######## MAIN COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOODE!!!  ##########
#######################################################################

#######################################################################
#### Step A: get info and set and/or initialize global variables
#######################################################################
home = os.getcwd()
comm = MPI.COMM_WORLD
my_id = comm.rank
num_procs = comm.Get_size()
generation = 0
inFile = open("INS_gum.txt",'r')
temp_lines = inFile.readlines()
l=0
atom_types = []
atom_amounts = []
num_atom_types = 0
num_atoms = 0
degen_threshhold = 0.0001
for line in temp_lines:
        #print "l is: " + str(l)
	if l == 0:
		useful_lines = temp_lines[l]
		new_lines = useful_lines.split()
		Title = ''
		for i in range(1,len(new_lines)):
			Title = Title + new_lines[i] + ' '
		l = l + 1
        if l == 1:
                useful_line = temp_lines[l]
                new_lines = useful_line.split()
                len_new_lines = len(new_lines)
                k = 1
                for k in range(1,len(new_lines)):
                        if k%2 != 1:
                                atom_amounts.append(int(new_lines[k]))
                                num_atoms = num_atoms + int(new_lines[k])
                        if k%2 == 1:
                                num_atom_types = num_atom_types + 1
                                for x in range(0,int(new_lines[k+1])):
                                        atom_types.append(new_lines[k])
                        k = k + 1
                l = l + 1
	if l == 2:
                useful_lines = temp_lines[l]
                new_lines = useful_lines.split()
                charge = new_lines[1]
                l = l + 1
        if l == 3:
                useful_lines = temp_lines[l]
                new_lines = useful_lines.split()
                mult = new_lines[1]
                l = l + 1
        if l == 4:
                useful_line = temp_lines[l]
                new_lines = useful_line.split()
                nSeeders = int(new_lines[1])
                l = l + 1
        if l == 5:
                useful_line = temp_lines[l]
                new_lines = useful_line.split()
                box_x = int(new_lines[1])
                #print str(box_x)
                box_y = int(new_lines[2])
                box_z = int(new_lines[3])
                l = l + 1
	if l == 6:
		useful_lines = temp_lines[l]
		new_lines = useful_lines.split()
		too_close = float(new_lines[1])
        if l == 7:
                useful_line = temp_lines[l]
                new_lines = useful_line.split()
                #print str(new_lines)
                MaxGenerations = int(new_lines[1])
                l = l + 1
        else:
                l = l + 1


denominator = 1
for amount in atom_amounts:
	denominator = denominator*math.factorial(amount)
num_distinguishables = math.factorial(num_atoms)/denominator

copies_path = home + '/gjf_log_copies'

comm.Barrier()
##################################################################
####Step A: complete
##################################################################



##################################################################
####Step B: write some basic info to output file
####        create/initialize some relevant directories and lists
##################################################################
if my_id == 0:
	Program_Output = open("Search_Output.txt",'w')
	Program_Output.write("********** Gumball Machine Hopping Search Output **********\n\nTitle: %s\nHomomorphs: %d\n\n" % (Title,num_distinguishables))
	Program_Output.close()
	gummed_winners_list = []
	morphed_winners_list = []
	Generation_Winners_List = []
	#directory to copy gjfs and logs to
	copies_path = home + '/gjf_log_copies'
	if not os.path.isdir(copies_path):
		os.mkdir(copies_path)

comm.Barrier()
#################################################################
####Step B: complete
#################################################################


#################################################################
####Step C: create num_proc random geometries and relax them
####        the nSeeder best will become the first geometric
####        seeds to be gum_machined 
#################################################################

#### each proc creates a random seed and runs it
box_size = int(num_atoms)/2
xyz_file = random_seed(box_size,box_size,box_size)
## change name of xyz_file to indicate its a trial seed
shutil.copy2(xyz_file,"init_seed%d.xyz"%my_id)
## run them, grab the energy 
gum_machine("init_seed%d.xyz"%my_id)
gjf_name = make_gjf("init_seed%d.xyz"%my_id,0)
log_name = run_gjf(gjf_name)
energy = get_energy(log_name)
## make them into trial_geom objects, set their energy values 
init_seed_object = trial_geom()
init_seed_object.set_names(gjf_name)
init_seed_object.set_energy(energy)
init_seed_object.set_geom_type(my_id)

comm.Barrier()
#print "DEBUG!!! past comm.Barrier() on line 773"
## proc 0 gathers all the energies from the other procs then evaluates which are the winners
if my_id != 0:
#	print "DEBUG!!! in if my_id!=0 on 776"
	comm.send(init_seed_object, dest=0, tag=(my_id*1+1))

if my_id == 0:
#	print "DEBUG!!! in if my_id==0 on 780"
	trial_seeds_list = []
	for i in range(1,num_procs):
		seed_object_i = comm.recv(source=i, tag=(i*1+1))
		trial_seeds_list.append(seed_object_i)
	
	#proc 0 adds its own to the list
	trial_seeds_list.append(init_seed_object)		

	seed_energies_list = []
	for seed_object in trial_seeds_list:
		seed_energies_list.append(seed_object.trial_energy)
	sorted_seeds_energies = sorted(seed_energies_list)

	seed_winners_list = []
	seed_winners_energies_list = []
	for i in range(len(sorted_seeds_energies)):
		#print "i is " + str(i)
		# we want to append seed_winners_list with the uniquely lowest energy isomers until it is the length of nSeeders
		if len(seed_winners_list) < int(nSeeders):
			if i == 0:
				for trial_seed in trial_seeds_list:
					if trial_seed.trial_energy == sorted_seeds_energies[i]:
#						print "about to append line 803"	
						seed_winners_list.append(trial_seed)
						seed_winners_energies_list.append(sorted_seeds_energies[i])
						break
			else:
				#check for degeneracy/uniqueness before appending
#				print "in else loop, i is %d" % i
				energy_i = sorted_seeds_energies[i]
				for energies in seed_winners_energies_list:
					energy_k = energies
					if math.fabs(energy_i-energy_k) > degen_threshhold:
						degen = 0
					else:
						degen = 1
						print "degeneracy found"
				if degen == 0:
					print "no degeneracy found, check geom_type. line 825"
					for trial_seed in trial_seeds_list:
						print "line 827.   looking for the correct object"
						if trial_seed.trial_energy == energy_i:
							print "coorect object found."
							print "about to check geom type. line 826"
							type_check = 0
							for seed_winner in seed_winners_list:
								if trial_seed.trial_geom_type != seed_winner.trial_geom_type:
									print "different geom type than other winner"
									type_check = 0
								else:
									print "same geom type as another winner, don't append"
									type_check = 1
									break
							if type_check == 0:
								seed_winners_list.append(trial_seed)
								seed_winners_energies_list.append(energy_i)
							else:
								continue

				elif i == (len(sorted_seeds_energies)-1):
#					print "degenerate, but no more unique structures, append"
					if trial_seed.trial_energy == energy_i:
						seed_winners_list.append(trial_seed)
						seed_winners_energies_list.append(energy_i)
						break

				else:
					pass
		
	
#	print "DEBUG 805.    trial_seeds_list is " + str(trial_seeds_list)
#	print "DEBUG 805.    sorted_seeds_energies is " + str(sorted_seeds_energies)	
	print "DEBUG 805.    seed_winners_list is " + str(seed_winners_list)	
			
	

	#copy the seeds to the winners list that is used in the main loop
	Generation_winners = []
	for winner in seed_winners_list:
		Generation_winners.append(winner)

	Program_Output = open("Search_Output.txt", 'a')
	winners_string = ''
	for winner in Generation_winners:
		winners_string = winners_string + "  %s  " % winner.trial_log
	Program_Output.write("*****  Initial Seed Winners are: %s \n" % winners_string)
	Program_Output.close()

comm.Barrier()
#########################################
####Step C: complete
#########################################

#####################################
########## MAIN LOOP ################
#####################################

for generation in range(MaxGenerations):
	print "beginning of main loop. generation is: " + str(generation)
###############################################################################################
######Step 1: gum_machine the previous Generations winners and create gjfs,
######         send out the gjfs to the appropriate procs, run them
###############################################################################################

	if my_id == 0:
		print "Step 1 if my_id == 0 line 870"
		## each winner should be gumballed num_proc/nSeeders times
		## if length of Generation_winners is less than nSeeders, copy the first item and append
		while len(Generation_winners) < nSeeders:
			Generation_winners.append(Generation_winners[0])
		gummed_objects_list = []
		for i in range(nSeeders):
			#get the Cartesian coords from the winner
			xyzFile = get_xyz(Generation_winners[i].trial_log)
			#print "xyzFile is: " + str(xyzFile)
			for k in range(num_procs/nSeeders):
				#create trial_geom object and set some of its info
				trial_object = trial_geom()
				trial_object.set_geom_type(i)
				gum_machine(xyzFile)
				ID_name = "gen%d_geom%d_num%d.xyz" % (generation,i,k)
				shutil.copy2(xyzFile,ID_name)
				gjf_name = make_gjf(ID_name,0)
				trial_object.set_names(gjf_name)
				gummed_objects_list.append(trial_object)

		#clean out the Generation_winners list
		Generation_winners = []

		## send out the gjfs
		for i in range(1,num_procs):
			comm.send(gummed_objects_list[i], dest=i, tag=i*2+1)
	
		## proc 0 makes its own trial_from runs the gjf, sets other trial_geom data and appends it to a list of sub-generation objects
		trial_object = trial_geom()
		gjf_name = gummed_objects_list[0].trial_gjf
		"DEBUG!!! 893 gjf_name is: " + gjf_name
		trial_object.set_names(gjf_name)
		trial_object.set_geom_type(0)
		log_file =  run_gjf(gummed_objects_list[0].trial_gjf)
		energy = get_energy(log_file)
		trial_object.set_energy(energy)
			
		sub_gen_objects_list1 = []	
		sub_gen_objects_list1.append(trial_object)	
			
		

	if my_id != 0:
		print "step 1 if my_id != 0 line 910"
		## receive gjfs, run them
		trial_object = comm.recv(source=0, tag=my_id*2+1)
		trial_object.set_genID(generation,my_id)
		log_file = run_gjf(trial_object.trial_gjf)
		energy = get_energy(log_file)
		trial_object.set_energy(energy)

	comm.Barrier()	
#########################################################################################
####Step 1: complete
#########################################################################################
	

#########################################################################################
####Step 2: collect the gum_machined logs, choose the sub-Generational winners
####        perturb them, send them out, run them
#########################################################################################

	if my_id != 0:
		print "beginning of Step 2 if my_id!=0"
		comm.send(trial_object, dest=0, tag=my_id*3+1)

	## append trial_geom objects from other procs to the sub_gen_logs_list1
	if my_id == 0:
		print "beginning of step 2 if my_id == 0"
		for i in range(1,num_procs):
			trial_object_i= comm.recv(source=i, tag=i*3+1)
			sub_gen_objects_list1.append(trial_object_i)

	        gummed_objects_energies_list = []
        	for gummed_object in sub_gen_objects_list1:
                	gummed_objects_energies_list.append(gummed_object.trial_energy)
        		sorted_gummed_objects_energies = sorted(gummed_objects_energies_list)

	        sub_gen_1_winners_list = []
        	sub_gen_1_winners_energies_list = []
	        for i in range(len(sorted_gummed_objects_energies)):
        	        if len(sub_gen_1_winners_list) < int(nSeeders):
				print "line 975"
                	        if i == 0:
					print "line 977"
                        	        for trial_object in sub_gen_objects_list1:
                                	        if trial_object.trial_energy == sorted_gummed_objects_energies[i]:
                                        	        print "about to append line 803"        
                                                	sub_gen_1_winners_list.append(trial_object)
                                             	  	sub_gen_1_winners_energies_list.append(sorted_gummed_objects_energies[i])

	                        else:
        	                        #check for degeneracy/uniqueness before appending
	                                print "in else loop, i is %d" % i
        	                        energy_i = sorted_gummed_objects_energies[i]
                	                for energies in sub_gen_1_winners_energies_list:
                        	                energy_k = energies
                                	        if math.fabs(energy_i-energy_k) > degen_threshhold:
                                        	        degen = 0
	                                        else:
        	                                        degen = 1
	               	                                print "degeneracy found"
                        	        if degen == 0:
                                	        for trial_object in sub_gen_objects_list1:
                                        	        if trial_object.trial_energy == energy_i:
                                        	       	        print "about to check geom_type. line 984"
								type_check = 0
								for sub_gen_1_winner in sub_gen_1_winners_list:
									"DEBUG line 992. trial_object.trial_geom_type is %s " % str(trial_object.trial_geom_type)
									"DEBUG line 993. sub_gen_1_winner.trial_geom_type is %s " % str(sub_gen_1_winner.trial_geom_type)
									if trial_object.trial_geom_type != sub_gen_1_winner.trial_geom_type:
										type_check = 0
										print "different type, possibly append, keep checking"
									else:
										print "matched geom type, do not append"
										type_check = 1
										break
								if type_check == 0:
									print "in if type_check == 0 line 1001"
                                                        		sub_gen_1_winners_list.append(trial_object)
                                                      			sub_gen_1_winners_energies_list.append(energy_i)
                                                        		break
								else:
									print "in else on line 1014."	
									continue
	                                elif i == (len(sorted_gummed_objects_energies)-1):
	                                        print "degenerate, but no more unique structures, append"
        	                                if trial_object.trial_energy == energy_i:
                	                                sub_gen_1_winners_list.append(trial_object)
                        	                        sub_gen_1_winners_energies_list.append(energy_i)
                                	                break

	                                else:
        	                                pass


		print "sorted_gummed_object_energies is: " + str(sorted_gummed_objects_energies)
		print "sub_gen_1_winners_energies_list is: " + str(sub_gen_1_winners_energies_list)
		print "sub_gen_1_winners_list is: " + str(sub_gen_1_winners_list)

		Program_Output = open("Search_Output.txt",'a')
		winners_string = ''
		for winner in sub_gen_1_winners_list:
			winners_string = winners_string + "  %s  " % winner.trial_log
		Program_Output.write("gum machined winners are: %s  \n" % winners_string)
		Program_Output.close()

	comm.Barrier()


	## we now have the winners from the gum balling. now perturb them and send them out
	if my_id == 0:
		print "in my_id == 0 at 990"
                ## each winner should be gumballed num_proc/nSeeders times
		## make sure sub_gen_1_winners_list is of length nSeeders. If not, copy the first item of the list and append
		while len(sub_gen_1_winners_list) < nSeeders:
			sub_gen_1_winners_list.append(sub_gen_1_winners_list[0])
                perturbed_objects_list = []
                for i in range(nSeeders):
                        #get the Cartesian coords from the winner
                        xyzFile = get_xyz(sub_gen_1_winners_list[i].trial_log)
                        print "xyzFile is: " + str(xyzFile)
                        for k in range(num_procs/nSeeders):
                                #create trial_geom object and set some of its info
                                trial_object = trial_geom()
                                trial_object.set_geom_type(i+nSeeders)
                                morphed_xyz_path = xyz_morph(xyzFile,"std","std",too_close)
                                ID_name = "gen%d_geom%d_num%d.xyz" % (generation,i+nSeeders,k)
                                shutil.copy2(morphed_xyz_path,ID_name)
                                gjf_name = make_gjf(ID_name,0)
                                trial_object.set_names(gjf_name)
                                perturbed_objects_list.append(trial_object)

                ## send out the new trial objects
                for i in range(1,num_procs):
                        comm.send(perturbed_objects_list[i], dest=i, tag=i*4+1)

                ## proc 0 makes its own trial_from runs the gjf, sets other trial_geom data and appends it to a list of sub-generation objects
                trial_object = trial_geom()
                gjf_name = perturbed_objects_list[0].trial_gjf
                "DEBUG!!! 893 gjf_name is: " + gjf_name
                trial_object.set_names(gjf_name)
                log_file =  run_gjf(perturbed_objects_list[0].trial_gjf)
                energy = get_energy(log_file)
                trial_object.set_energy(energy)

                sub_gen_objects_list2 = []
                sub_gen_objects_list2.append(trial_object)



        if my_id != 0:
		print "in if my_id != 0 at 1027"
                ## receive gjfs, run them
                trial_object = comm.recv(source=0, tag=my_id*4+1)
                trial_object.set_genID(generation,my_id)
		print "about to run gjf"
                log_file = run_gjf(trial_object.trial_gjf)
                energy = get_energy(log_file)
                trial_object.set_energy(energy)
		print "proc %d end of Step2" % my_id

        comm.Barrier()
#############################################################################################
####Step 2: complete
#############################################################################################



#############################################################################################
####Step 3: collect the perturbed structures, choose the overall Generational winners
#############################################################################################

	if my_id != 0:
		print "beginning of step 3"
        	comm.send(trial_object, dest=0, tag=my_id*5+1)

	## append trial_geom objects from other procs to the sub_gen_logs_list1
	if my_id == 0:
        	for i in range(1,num_procs):
                	trial_object_i= comm.recv(source=i, tag=i*5+1)
                	sub_gen_objects_list2.append(trial_object_i)

                	perturbed_objects_energies_list = []
	                for perturbed_object in sub_gen_objects_list2:
        	                perturbed_objects_energies_list.append(perturbed_object.trial_energy)
                	        sorted_perturbed_objects_energies = sorted(perturbed_objects_energies_list)
		

#		print "DEBUG!! line 1068. sorted_perturbed_objects_energies is: " + str(sorted_perturbed_objects_energies)
	        sub_gen_2_winners_list = []
        	sub_gen_2_winners_energies_list = []
	        for i in range(len(sorted_perturbed_objects_energies)):
        	        if len(sub_gen_2_winners_list) < int(nSeeders):
                	        if i == 0:
                        	        for trial_object in sub_gen_objects_list2:
                                	        if trial_object.trial_energy == sorted_perturbed_objects_energies[i]:
         	                              	        print "about to append line 803. i was 0"        
                                                	sub_gen_2_winners_list.append(trial_object)
                                                	sub_gen_2_winners_energies_list.append(sorted_perturbed_objects_energies[i])

		              	else:
	                                #check for degeneracy/uniqueness before appending
	                                print "in else loop, i is %d" % i
        	                        energy_i = sorted_perturbed_objects_energies[i]
                        	        for energies in sub_gen_2_winners_energies_list:
                	                        energy_k = energies
                                	        if math.fabs(energy_i-energy_k) > degen_threshhold:
                                        	        degen = 0
                                       		else:
                                                	degen = 1
                                            		print "degeneracy found"
                                if degen == 0:
					print "in degen = 0"
                                        for trial_object in sub_gen_objects_list2:
                                                if trial_object.trial_energy == energy_i:
                                                        print "about to append line 832"
                                                        sub_gen_2_winners_list.append(trial_object)
                                                        sub_gen_2_winners_energies_list.append(energy_i)
                                                        break

                                elif degen == 1:
					print "degeneracy found"
					if i == (len(sorted_perturbed_objects_energies)-1):
                                        	print "degenerate, but no more unique structures, append"
                                        	if trial_object.trial_energy == energy_i:
							print "line 1136. found object, now check its geom_type"
							type_check = 0
							for winner in sub_gen_2_winners_list:
								if trial_object.trial_geom_type != winner.trial_geom_type:
									print "not a geom_type match. maybe append, continue checking"
									type_check = 0
								else:
									print "same geom type, don't append"
									type_check = 1
									break
							if type_check == 0:
								print "no geom type matches, append winners list"
                                                		sub_gen_2_winners_list.append(trial_object)
                                                		sub_gen_2_winners_energies_list.append(energy_i)
                                                		break

                                			else:
								print "in else. geom_type was the same. don't append. line 1152"
                                        			pass


	        print "sorted_perturbed_object_energies is: " + str(sorted_perturbed_objects_energies)
       		print "sub_gen_2_winners_energies_list is: " + str(sub_gen_2_winners_energies_list)
        	print "sub_gen_2_winners_list is: " + str(sub_gen_2_winners_list)
       
		Program_Output = open("Search_Output.txt", 'a')
		winners_string = ''
		for winner in sub_gen_2_winners_list:
			winners_string = winners_string + "  %s  " % winner.trial_log
		Program_Output.write("Perturbed Winners are: %s  \n" % winners_string)
		Program_Output.close()

		sub_generations_winners_energies_list = sub_gen_2_winners_energies_list + sub_gen_1_winners_energies_list
		print "sub_gen_winners_energies is: " + str(sub_generations_winners_energies_list)
		sub_generations_winners_sorted_energies_list = sorted(sub_generations_winners_energies_list)
		sub_generation_winners_objects_list = sub_gen_2_winners_list + sub_gen_1_winners_list
		print "sub_gen_winners_objects_lst is: " + str(sub_generation_winners_objects_list)

		Generation_winners = []
		Generation_winners_energies_list = []
		for i in range(len(sub_generations_winners_sorted_energies_list)):
                	        if len(Generation_winners) < int(nSeeders):
                        	        if i == 0:
                                	        for trial_object in sub_generation_winners_objects_list:
		                             	        if trial_object.trial_energy == sub_generations_winners_sorted_energies_list[i]:
                	                               	        print "about to append line 803"        
                                                        	Generation_winners.append(trial_object)
                                                        	Generation_winners_energies_list.append(sub_generations_winners_energies_list[i])


                                	else:
                                        	#check for degeneracy/uniqueness before appending
    		    #                               print "in else loop, i is %d" % i
                	                        energy_i = sub_generations_winners_sorted_energies_list[i]
                        		        for energies in Generation_winners_energies_list:
                                        	        energy_k = energies
                                                	if math.fabs(energy_i-energy_k) > degen_threshhold:
                                                        	degen = 0
                                        		else:
                                                 		degen = 1
#                                            	  		 print "degeneracy found"
                        	   		if degen == 0:
                                			for trial_object in sub_generation_winners_objects_list:
                                        			if trial_object.trial_energy == energy_i:
#                                               		        print "about to append line 832"
                                                        		Generation_winners.append(trial_object)
                                                  	        	Generation_winners_energies_list.append(energy_i)
                                                        		break

						elif i == (nSeeders - 1):
#                               		        print "degenerate, but no more unique structures, append"
                                       			 if trial_object.trial_energy == energy_i:
                                                		Generation_winners.append(trial_object)
                                               			Generation_winners_energies_list.append(energy_i)
                                                		break

		print "Generation_winners is: " + str(Generation_winners)
		
		Program_Output = open("Search_Output.txt",'a')
		winners_string = ''
		for winner in Generation_winners:
			winners_string = winners_string + '  %s  ' %winner.trial_log
		Program_Output.write("*****  Generation winners are: %s  \n" % winners_string)
		Program_Output.close()

	comm.Barrier()
################################################################################
####Step 3: complete
################################################################################

print "all generations complete."



