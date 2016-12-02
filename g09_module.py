import os
import sys
import subprocess


def make_gjf(xyzFile,mod):
	#mod can be 0 or 1 to choose an alternate route
	inFile = open("INS_g09.txt",'r')
	g09_lines = inFile.readlines()
	g = 0
	for line in g09_lines:
		if g == 1:
			split_lines = g09_lines[g].split()
			charge = split_lines[1]
			g += 1
		if g == 2:
			split_lines = g09_lines[g].split()
			mult = split_lines[1]
			g += 1
		if g == 3:
			Route = g09_lines[g].lstrip('Route:')
			Route = Route.lstrip()
			g += 1
		if g == 4:
			altRoute = g09_lines[g].lstrip('altRoute:')
			altRoute = altRoute.lstrip()
			g += 1
		if g == 5:
			split_lines = g09_lines[g].split()
			mem = split_lines[1]
			g += 1
		if g == 6:
			split_lines = g09_lines[g].split()
			nproc = split_lines[1]
			g += 1
		else:
			g += 1
	inFile.close()


	inFile2 = open(xyzFile,'r')
	gjf_name = xyzFile[:-4] + ".gjf"
	outFile = open(gjf_name,'w')

	lines = []
	lines.append('%mem=' + mem + '\n')
	lines.append('%nproc=' + nproc + '\n')
	if mod == 0:
		lines.append('# ' + Route + '\n')
	elif mod == 1:
		lines.append('# ' + altRoute + '\n')
	lines.append('structure ' + xyzFile + '\n\n')
	lines.append(charge + ' ' + mult + '\n')
	while 1:
		line = inFile2.readline()
		if not line:
			break
		lines.append(' ' + str(line))
	lines.append('\n')

	outFile.writelines(lines)
	
	inFile2.close()
	outFile.close()
	



def make_PBS_script(jobname):
	#works for CCV computer
	inFile = open("INS_g09.txt",'r')
	g09_lines = inFile.readlines()
	g = 0
	for line in g09_lines:
		if g == 8:
			split_lines = g09_lines[g].split()
			walltime = split_lines[1]
			g += 1
		if g == 9:
			split_lines = g09_lines[g].split()
			PBS_procs = split_lines[1]
			g += 1
		else:
			g += 1
	inFile.close()

	jobScript = open(jobname + ".job",'w')
	lines = []
	lines.append("#!/bin/csh \n#PBS -N \"%s\"\n#PBS -l %s\n"%(jobname,PBS_procs))
	lines.append("#PBS -j oe\n#PBS -l walltime=%s\n#PBS -o %s.out\n\n"%(walltime,jobname))
	lines.append("cd $PBS_O_WORKDIR\n\n")
	lines.append("module load gaussian/g09\n\n")
	lines.append("setenv DIR \"/gpfs/scratch/$USER/%s\"\n"%jobname)
	lines.append("mkdir -p $DIR\n")
	lines.append("cp $PBS_O_WORKDIR/%s.gjf $DIR\n"%jobname)
	lines.append("cd $DIR\n")
	lines.append("g09 < %s.gjf >& $PBS_O_WORKDIR/%s.log\n\n"%(jobname,jobname))
	lines.append("rm *.chk *.rwf *.int *.scr *.d2e *.int\n")	

	jobScript.writelines(lines)

	jobScript.close()


def run_gjf(name,run_mod):
	#run mod decides if to submit to PBS queue or run locally
	if run_mod == 0:
		#use PBS queue
		#NEED TO SET UP TO CHECK OUTPUT AND RERUN ACCORDINGLY!!!
		script_name = name + ".job"
		commands = "qsub " + script_name
		subprocess.Popen(commands, shell=True).wait()
	elif run_mod == 1:
		#use local resources
		print "not yet set up to run 09 on local resources. error."
		sys.exit()


def check_log(gjf_log):
	#checks gaussian log file for errors
	#returns coordinates
	#returns energy

	inFile = open(gjf_log,'r')
	
        okay = -1
        text = inFile.read()
        search = 'Normal termination of Gaussian'
        index = text.find(search)
        okay = index
        if okay != -1:
		#print 'Normal termination of Gaussian found'
                check = 0
        else:
                search = 'termination request processed by link 9999'
                index = text.find(search)
                okay = index
                if okay != -1:
                        check = 1
		else:
			search = 'galloc:  could not allocate memory.'
			index = text.find(search)
			okay = index
			if okay != -1:
				check = 3
			else:
				search = 'Error in internal coordinate system.'
				index = text.find(search)
				okay = index
				if okay != -1:
					check = 4
                		else:
                        		print 'gaussian error during optimization'
                      			check = 2
	inFile.close()
	
	#get energy if appropriate else return generic energy of 0.0
	if check == 0 or check == 1:
		#print "gjf_log is: " + str(gjf_log)
		inFileA = open(gjf_log,'r')
		text = inFileA.read()
		outPath = gjf_log[:-4]+'_energy_temp'
		outFile = open(outPath,'w')
		lines = []
		x = text[text.rfind('SCF Done'):]
		#print "x is: " + str(x)
		lines.append(text[text.rfind('SCF Done'):])
	
		#print "lines is " + str(lines)
			
		outFile.writelines(lines)
		outFile.close()
		
		inFile2 = open(outPath,'r')
		outPath2 = gjf_log[:gjf_log.find('.log')]+'_energy'
		outFile2 = open(outPath2,'w')

		lines = []
		text = inFile2.read()
		lines.append(text[:text.find('cycles')])

		outFile2.writelines(lines)
		outFile2.close()
		inFile2.close()

		inFile3 = open(outPath2,'r')
		raw_lines = inFile3.readlines()
		##print "raw_lines is: " + str(raw_lines)
		new_lines = []

		#loop through raw_lines and strip everything except what we want
		for line in raw_lines:
			new_line = line.lstrip()
			new_line = new_line.split()
			if len(new_line) > 0:
				new_lines.append(new_line[4])

		#print "new_lines is: " + str(new_lines)
		write_string = ''

		#build the string
		for line in new_lines:
			write_string += line
	
		energy = float(write_string)

		#remove extraneous files
		inFile3.close()
		os.remove(outPath)
		os.remove(outPath2)
		inFileA.close()
	
	else:
		energy = 0.0
		

	#get new coordinates, return coordinate file
	inFileB = open(gjf_log,'r')	
	outPath = gjf_log[:-4]+'_tmp.xyz'
	outFile = open(outPath,'w')
	
	lines = []
	text = inFileB.read()
	lines.append(text[text.rfind('Standard orientation'):text.rfind('Rotational constants')])

	outFile.writelines(lines)
	outFile.close()
	
	inFile2 = open(outPath,'r')
	outPath2 = outPath[:outPath.find('_tmp.xyz')]+'.xyz'
	outFile2 = open(outPath2,'w')

	raw_lines = inFile2.readlines()
	temp_lines = []
	new_lines = []
	i = 0
	#loop through raw lines and remove first five lines and last line, then store in tmp file
	for line in raw_lines:
		if i > 4:
			if line.find('--') == -1:
				temp_lines.append(line)
		i += 1
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
	

	#write it dude
	outFile2.write(write_string)
	outFile2.close()

	os.remove(outPath)

	#outPath2 is name of the xyz file
	xyzFile = outPath2

	inFileB.close()

	return [check,energy,xyzFile]


