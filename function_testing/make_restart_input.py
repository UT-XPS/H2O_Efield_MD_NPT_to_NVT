import csv

## This is a tool for creating the new .dat file for an NVT run,
## using the results of a previous NPT run. Specifically, it reads
## input from three files:
##
## water.dat     => the original .dat file
## finalpos.dump => the positions and velocities of the atoms in
##                  the final timestep
## log.lammps    => contains the thermo output
##
## 1) Using the data in log.lammps, the program calculates the
##    average x_len, y_len and z_len of the system during the
##    timesteps defined by the input variables first_step and
##    last_step.
## 2) It reads everything in water.dat upto the start of the
##    'Atoms' section. It edits the box size definitions to
##    match the average x_len, y_len and z_len determined
##    previously. (In most simulations, these should be the
##    same.) All of these lines are added to output_list
##    (initiated as an empty list).
## 3) It reads the lines defining the atoms in the finalpos.dump
##    file, and creates two lists of lines. The first list is the
##    new list of Atoms. The positions of the atoms in scaled
##    coordinates are converted into positions in distance units
##    by using the values of x_len, y_len and z_len. This list
##    is added to output_list. The second list is the new list
##    of atom velocities.
## 4) It reads the sections containing the definitions of bonds
##    and angles from the original .dat file, and adds them to
##    output_list.
## 5) The list containing atom velocities is added to output_list
##    and the lines in output_list are written to a new .dat file.

def make_restart_input(old_dat_path,old_dat_delimm,finalpos_path,finalpos_delimm,log_path,log_delimm,output_path,first_step,last_step):
	output_list = []
	ave_lengths = get_ave_lenghts(log_path,log_delimm,first_step,last_step)
	new_header = get_new_header(old_dat_path,old_dat_delimm,ave_lengths)
	output_list += new_header
	new_atoms = get_new_atoms(finalpos_path,finalpos_delimm,ave_lengths)
	output_list += ["Atoms\n","\n"]
	output_list += new_atoms
	bonds_and_angles = get_bonds_angles(old_dat_path,old_dat_delimm)
	output_list += ["\n"]
	output_list += bonds_and_angles
	velocities = get_new_velocities(finalpos_path,finalpos_delimm)
	output_list += ["\n","Velocities\n","\n"]
	output_list += velocities
	output(output_path,output_list)

def get_lines(file_path,delimm):
	file = open(file_path,'r')
	reader = csv.reader(file,delimiter=delimm)
	lines = []
	for line in reader:
		lines += [line]
	return lines

def get_ave_lenghts(file_path,delimm,first_step,last_step):
	log_lines = get_lines(file_path,delimm)
	used_lines = []
	use = False
	for line in log_lines:
		if use == False:
			if str(first_step) in line:
				nospace_line = []
				for element in line:
					if element != '':
						nospace_line += [element]
				if nospace_line[0] == str(first_step):
					use = True
		if use == True:
			nospace_line = []
                        for element in line:
                                if element != '':
                                        nospace_line += [element]
			if nospace_line[0] != 'SHAKE':
				if int(nospace_line[0]) > last_step:
					use = False
					print "WARNING! Last step not found"
		if use == True:
			nospace_line = []
			for element in line:
				if element != '':
					nospace_line += [element]
			if nospace_line[0] not in ['SHAKE','1']:
				used_lines += [nospace_line]
				if int(nospace_line[0]) == last_step:
					use = False
	xvals = []
	yvals = []
	zvals = []
	for line in used_lines:
		xvals += [float(line[8])]
		yvals += [float(line[9])]
		zvals += [float(line[10])]
	ave_x = average(xvals)
	ave_y = average(yvals)
	ave_z = average(zvals)
	ave_lengths = [ave_x,ave_y,ave_z]
	return ave_lengths

def average(list):
	sum = 0
	for element in list:
		sum += element
	average = sum/float(len(list))
	return average

def get_new_header(file_path,delimm,ave_lengths):
	old_dat_lines = get_lines(file_path,delimm)
	new_header_lines = []
	use = True
	for line in old_dat_lines:
		if len(line) > 0:
			if line[0] == 'Atoms':
				use = False
		if use == True:
			new_header_lines += [line]
	index = 0
	for line in new_header_lines:
		if 'xlo xhi' in line:
			new_header_lines[index] = [format(0.0,'.6f'),format(ave_lengths[0],'.6f'),'xlo xhi']
                if 'ylo yhi' in line:
                        new_header_lines[index] = [format(0.0,'.6f'),format(ave_lengths[1],'.6f'),'ylo yhi']
                if 'zlo zhi' in line:
                        new_header_lines[index] = [format(0.0,'.6f'),format(ave_lengths[2],'.6f'),'zlo zhi']
		index += 1
	new_header = put_delimiters_back(new_header_lines,delimm)
	return new_header

def put_delimiters_back(lines_list,delimm):
	delimlines_list = []
	for line in lines_list:
		newstring = ''
		for element in line:
			if newstring != '':
				newstring += delimm
			newstring += str(element)
		newstring += '\n'
		delimlines_list += [newstring]
	return delimlines_list

def get_new_atoms(file_path,delimm,ave_lengths):
	finalpos_lines = get_lines(file_path,delimm)
	atoms_reversed = []
	use = True
	for line in reversed(finalpos_lines):
		if line[0] == 'ITEM:':
			use = False
		if use == True:
			atom = line[:-4]
			atoms_reversed += [atom]
	atoms = []
	for atom in reversed(atoms_reversed[:]):
		atoms += [atom]
	xlen = ave_lengths[0]
	ylen = ave_lengths[1]
	zlen = ave_lengths[2]
	scaled_atoms = []
	for atom in atoms:
		scaled_atom = atom[:-3]
		scaled_x = float(atom[-3])*xlen
		scaled_y = float(atom[-2])*ylen
		scaled_z = float(atom[-1])*zlen
		scaled_atom += [format(scaled_x,'.6f'),format(scaled_y,'.6f'),format(scaled_z,'.6f')]
		scaled_atoms += [scaled_atom]
	new_atoms = put_delimiters_back(scaled_atoms,delimm)
	return new_atoms

def get_bonds_angles(file_path,delimm):
	old_dat_lines = get_lines(file_path,delimm)
	bonds_angles_lines = []
	use = False
	for line in old_dat_lines:
		if use == False:
			if len(line) != 0:
				if line[0] == 'Bonds':
					use = True
		if use == True:
			if len(line) != 0:
				if len(line[0]) != 0:
					if (line[0] != 'Bonds') and (line[0] != 'Angles'):
						if line[0][0].isalpha() == True:
							use = False
		if use == True:
			bonds_angles_lines += [line]
	bonds_and_angles = put_delimiters_back(bonds_angles_lines,delimm)
	return bonds_and_angles

def get_new_velocities(file_path,delimm):
	finalpos_lines = get_lines(file_path,delimm)
	velocities_reversed = []
	use = True
	for line in reversed(finalpos_lines):
		if line[0] == 'ITEM:':
			use = False
		if use == True:
                        nospace_line = []
                        for element in line:
                                if element != '':
                                        nospace_line += [element]
			velocity = [nospace_line[0]]+nospace_line[-3:]
			velocities_reversed += [velocity]
	velocities = []
	for velo in reversed(velocities_reversed[:]):
		velocities += [velo]
	new_velocities = put_delimiters_back(velocities,delimm)
	return new_velocities

def output(file_path,output_list):
	file = open(file_path,"w")
	file.writelines(output_list)
	file.close()
