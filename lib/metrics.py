import mdtraj as md
import MDAnalysis as mda
#from MDAnalysis.analysis import contacts
import MDAnalysis.analysis.rms
#import matplotlib as plt
#from matplotlib import pyplot as plt
#import os
#import sys
import numpy as np
import csv
import pandas
import statistics

########################################################## Compute the least square methods on the data list provided
########################################################## called by other metrics functions
def slope(values_metric, mothfolder):
	
	# create a dict from a list
	data = {k: v for v, k in enumerate(values_metric)}
	print(data, type(data))
	time= [] # needed for interpolation
	dist = []
	nume = []
	denume = []
	for key, value in data.items():
		dist.append(float(key))
		time.append(float(value))
	time_max = max(time)
	for key, value in data.items():
		if time_max == value:
			distance = key
	#print(distance)

	mean_t = np.mean(time) # needed for interpolation
	mean_d = np.mean(dist)
	for key, value in data.items():
		product = (float(value) - mean_t)*(float(key) - mean_d)
		nume.append(product)
		simple_diff = (float(value) - mean_t)**2
		denume.append(simple_diff)
	slope = float(np.sum(nume))/float(np.sum(denume))
	
	# write the last distance in a crude log file
	with open(mothfolder+'/crude_last_values.log' , 'a') as logF:
		logF.write('%s\n' %str(distance))
	
	logF.close()
		
	return slope
	

########################################################## Used by distance() function below
def compute_center_of_mass(traj, select=None):
    """Compute the center of mass for each frame.
    Parameters
    ----------
    traj : Trajectory
        Trajectory to compute center of mass for
    select : str, optional, default=all
        a mdtraj.Topology selection string that
        defines the set of atoms of which to calculate
        the center of mass, the default is all atoms
    Returns
    -------
    com : np.ndarray, shape=(n_frames, 3)
         Coordinates of the center of mass for each frame
    """
# create a np array
    com = np.empty((traj.n_frames, 3))

# get alll the atoms in the system
    if select is None:
        masses = np.array([a.element.mass for a in traj.top.atoms])
        masses /= masses.sum()

        xyz = traj.xyz

# get masses of selected atoms
    else:
        atoms_of_interest = traj.topology.select(select)

        masses = np.array([traj.top.atom(i).element.mass for i in atoms_of_interest])
        masses /= masses.sum()

# compute COM
        xyz = traj.xyz[:, atoms_of_interest]

    for i, x in enumerate(xyz):
        com[i, :] = x.astype('float64').T.dot(masses)
    return com





#################################################### Compute distances between COMs
def distance(par, sel_1, sel_2, n, mothfolder):


# Debug: log with crude values for score computation
	def logAll(distances,mean_distance,last_distance,distMetric):
		with open(mothfolder+'/distance_crude.log','a') as distF:
			distF.write('distances: %s, mean: %s, last: %s, score: %s\n' %(" ".join(map(str, distances)),str(mean_distance), str(last_distance), str(distMetric)))
			
	if par['MDEngine'] == 'ACEMD':
		xtc = '%s_%s_wrapped.xtc' %(par['Output'], str(n))
		
		if par['Forcefield'] == 'CHARMM':
			topology = '%s.psf' %par['Topology']
		elif par['Forcefield'] == 'AMBER':
			topology = '%s.prmtop' %par['Topology']
	

	elif par['MDEngine'] == 'GROMACS':
		xtc = '%s_%s.xtc' %(par['Output'], str(n))
		topology = '%s.gro' %par['Topology']

		
	traj  = md.load(xtc, top=topology)

# let's decide what selections to use according to the number of metrics we wan to compute
	#if par['NumberCV'] == '1' and par['Metric_1'] == 'Distance':
	c1 = compute_center_of_mass(traj, select=sel_1)
	c2 = compute_center_of_mass(traj, select=sel_2)
		
	#elif par['NumberCV'] == '2' and par['Metric_1'] == 'Distance':
	#	c1 = compute_center_of_mass(traj, select=selection_list[0])
	#	c2 = compute_center_of_mass(traj, select=selection_list[1])	
		
	#elif par['NumberCV'] == '2' and par['Metric_2'] == 'Distance':
	#	c1 = compute_center_of_mass(traj, select=selection_list[2])
	#	c2 = compute_center_of_mass(traj, select=selection_list[3])	

	distances = []

# compute distance between elements sam eposiiton in 2 different lists	
	for a,b in zip(c1, c2):
		D = np.linalg.norm(a-b)*10
		distances.append(D)


	n = len(distances)
	mean_distance = sum(distances)/n
	last_distance = distances[-1]
	
	if par['NumberCV'] == '1': # we return the walker score or slope already
		if par['Slope'] == 'NO': # we return the walker score
			distMetric = (mean_distance*last_distance)**0.5
			logAll(distances,mean_distance,last_distance,distMetric)
			return distMetric

		elif par['Slope'] == 'YES': #we return the walker slope
			distMetric = slope(distances, mothfolder)
			return distMetric

	elif par['NumberCV'] == '2': # we return all the metric values and the last distance
		return distances, last_distance
		
##################################################
def contacts(par, sel_1, sel_2, n, mothfolder):
	#import MDAnalysis as mda
	from MDAnalysis.analysis import contacts

# Debug: log with crude values for score computation
	def logAll(timeseries,mean_contacts,last_contacts,distMetric):
		with open(mothfolder+'/contacts_crude.log','a') as distF:
			distF.write('contacts: %s, mean: %s, last: %s, score: %s\n' %(" ".join(map(str, timeseries)),str(mean_contacts), str(last_contacts), str(distMetric)))


	if par['MDEngine'] == 'ACEMD':
		xtc = '%s_%s_wrapped.xtc' %(par['Output'], str(n))
	
		if par['Forcefield'] == 'CHARMM':
			psf = '%s.psf' %par['Topology']
		elif par['Forcefield'] == 'AMBER':
			psf = '%s.prmtop' %par['Topology']

	elif par['MDEngine'] == 'GROMACS':
		xtc = '%s_%s.xtc' %(par['Output'], str(n))
		psf = '%s.gro' %par['Topology']



	u = mda.Universe(psf, xtc)

# let's decide what selections to use according to the number of metrics we wan to compute
	#if par['NumberCV'] == '1' and par['Metric_1'] == 'Contacts':
	sel_1 = u.select_atoms(sel_1)
	sel_2 = u.select_atoms(sel_2)
		
	#elif par['NumberCV'] == '2' and par['Metric_1'] == 'Contacts':
	#	sel_1 = u.select_atoms(selection_list[0])
	#	sel_2 = u.select_atoms(selection_list[1])	
		
	#elif par['NumberCV'] == '2' and par['Metric_2'] == 'Contacts':
	#	sel_1 = u.select_atoms(selection_list[2])
	#	sel_2 = u.select_atoms(selection_list[3])


	timeseries = []
	for ts in u.trajectory:
# calculate distances between sel_1 and sel_2
		dist = contacts.distance_array(sel_1.positions, sel_2.positions)
# determine which distances <= radius
		n_contacts = contacts.contact_matrix(dist, 3.5).sum()
#timeseries.append([ts.frame, n_contacts])
		timeseries.append(n_contacts)

	n = len(timeseries)
	mean_contacts = sum(timeseries)/n
	last_contacts = timeseries[-1]


	if par['NumberCV'] == '1': # we return the walker score or slope already
		if par['Slope'] == 'NO': # we return the walker score
			distMetric = (mean_contacts*last_contacts)**0.5	
			logAll(timeseries,mean_contacts,last_contacts,distMetric)
			return distMetric

		elif par['Slope'] == 'YES': # we return the walker slope
			distMetric = slope(timeseries, mothfolder)
			return distMetric
		
	elif par['NumberCV'] == '2': # we return all the metric values and the last distance
		return timeseries, last_contacts

################################################
def RMSD(par, sel_1, sel_2, n, mothfolder):
	#import MDAnalysis
	#import MDAnalysis.analysis.rms

# Debug: log with crude values for score computation
	def logAll(data,mean_rmsd,last_rmsd,distMetric):
		with open(mothfolder+'/rmsd_crude.log','a') as distF:
			distF.write('RMSD: %s, mean: %s, last: %s, score: %s\n' %(" ".join(map(str, data)), str(mean_rmsd), str(last_rmsd), str(distMetric)))


	if par['MDEngine'] == 'ACEMD':
		xtc = '%s_%s_wrapped.xtc' %(par['Output'], str(n))

		if par['Forcefield'] == 'CHARMM':
			psf = '%s.psf' %par['Topology']
		elif par['Forcefield'] == 'AMBER':
			psf = '%s.prmtop' %par['Topology']

	elif par['MDEngine'] == 'GROMACS':
		xtc = '%s_%s.xtc' %(par['Output'], str(n))
		psf = '%s.gro' %par['Topology']
		



	pdb = par['REFERENCE']
	
	u = MDAnalysis.Universe(psf,xtc)
	ref = MDAnalysis.Universe(pdb,pdb) 


# let's decide what selections to use according to the number of metrics we wan to compute
	#if par['NumberCV'] == '1' and par['Metric_1'] == 'RMSD':
	R = MDAnalysis.analysis.rms.RMSD(u, ref, select="%s" %sel_1, groupselections=["%s" %sel_2])                                    

	#elif par['NumberCV'] == '2' and par['Metric_1'] == 'RMSD':
	#	R = MDAnalysis.analysis.rms.RMSD(u, ref, select="%s" %selection_list[0], groupselections=["%s" %selection_list[1]])
	
	#elif par['NumberCV'] == '2' and par['Metric_2'] == 'RMSD':
	#	R = MDAnalysis.analysis.rms.RMSD(u, ref, select="%s" %selection_list[2], groupselections=["%s" %selection_list[3]])   



	R.run()
	rmsd = R.rmsd.T   # transpose makes it easier for plotting
	#time = rmsd[1]
	data = list(rmsd[3])


	n = len(data)
	mean_rmsd = sum(data)/n
	last_rmsd = data[-1]


	if par['NumberCV'] == '1': # we return the walker score or slope already
		if par['Slope'] == 'NO': # we return the walker score
			distMetric = (mean_rmsd*last_rmsd)**0.5	
			logAll(data,mean_rmsd,last_rmsd,distMetric)
			return distMetric
			
		elif par['Slope'] == 'YES': # we return the walker slope
			distMetric = slope(data, mothfolder)
			return distMetric

	elif par['NumberCV'] == '2': # we return all the metric values and the last distance
		return data,last_rmsd 


############################################
def HB_score(par, selection_list, n, mothfolder):

# to get resnam and atom involved in hydrogen bonds
	def label(hbond):
		hbond_label = '%s:%s' % (traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
		return hbond_label

# write HB for each frame	
	def logHB(Frame2hbond):
		with open(mothfolder+'/hydrogen_bonds.log', 'a') as hbs:
			print(Frame2hbond, file=hbs)
# debug info about the paqrametrs for HBì_score
	def logScores_detailed(index, value, numHB_frame,HBscore):
		with open(mothfolder+'/HB_scores_crude.log', 'a') as scoreF:
			scoreF.write('number wat mols: %s, number HB: %s, HB_score: %s\n' %(str(value), str(numHB_frame), str(HBscore)))		

########

	xtc = '%s_%s_wrapped.xtc' %(par['Output'], str(n))
	if par['Forcefield'] == 'CHARMM':
		psf = '%s.psf' %par['Topology']
	elif par['Forcefield'] == 'AMBER':
		psf = '%s.prmtop' %par['Topology']

	traj = md.load(xtc, top=psf)

	ligand = traj.topology.select(par['ligand_HB'])
	protein = traj.topology.select("protein")
# figure out how many frames you loaded, this is how many frames we will look for hbonds in
	n_frames = len(traj)
	#print(n_frames)

# This set will give us all of the unique hbonds that are made with the ligand, without repeats
# We will want to have this later so we make it not to avoid repeating hbond calculation 
	all_hbonds_set = set()
# This list will store all of the hbonds made per frame
	hbonds_each_frame = []

# We want to create a dictionary containing every frame and the ligand hbonds which occur in that frame
	Frame2hbond = {}
	for frame in range(n_frames):
    # The dictionary "words" are the frame number
		Frame2hbond[frame] = []
    # We are doing the hbond analysis frame by frame
		hbonds = md.baker_hubbard(traj[frame])
        #print(hbonds)
		hbonds_each_frame.append(hbonds)
        #print(hbonds_each_frame)
    # We only care about the hbonds if they involve the ligand 
		for hbond in hbonds:
			if ((hbond[0] in ligand) and (hbond[2] in protein) or (hbond[2] in ligand) and (hbond[0] in protein)): #ligand is donating or accepting 
				hbond_label= label(hbond) # get the atom names of the ligand involved 
				#all_hbonds_set.add(tuple(hbond))
				all_hbonds_set.add(hbond_label)
            # The dictionary "definitions" are all the hbonds in that frame
                        #Frame2hbond[frame].append(tuple(hbond))
				Frame2hbond[frame].append(hbond_label)

	print(Frame2hbond)
	logHB(Frame2hbond)
# let's get just the atom names of the ligand involved without repetitions
	lig_hetatm = []
	for frame,hbs in Frame2hbond.items():
		for hb in hbs:
			atom = hb.split(':')
			if 'ZMA' in atom[0]:
				lig_hetatm.append(atom[0].split('-')[1])
			elif 'ZMA' in atom[1]:
				lig_hetatm.append(atom[1].split('-')[1])
	lig_hetatm = list(set(lig_hetatm))
	LIG_HETATM = ' '.join(lig_hetatm)
	print(LIG_HETATM)

# get number of HB in each frame
	numHB_frame = []
	for frame,hbs in Frame2hbond.items():
		numHB_frame.append(len(hbs))
	print(numHB_frame)

# let's get contacts between water and ligand's hetatm involved 
	from MDAnalysis.analysis import contacts
	u = mda.Universe(psf, xtc)
	
	if len(LIG_HETATM) > 0: # there are hydrogen bonds
		sel_hetatm_lig = "%s and name %s" %(par['ligand_HB'], LIG_HETATM)
		
	else: # no hydrogen bonds, then use heteroatoms
		sel_hetatm_lig = "%s and (name O* or name N*) " %par['ligand_HB']
		
	if par['Forcefield'] == 'CHARMM':
		sel_water_oxy = "resname TIP3 and name O*"
	elif par['Forcefield'] == 'AMBER':
		sel_water_oxy = "resname WAT and name O*"
	lig = u.select_atoms(sel_hetatm_lig)
	wat = u.select_atoms(sel_water_oxy)

	timeseries = []
	for ts in u.trajectory:
# calculate distances between group_a and group_b
		dist = contacts.distance_array(lig.positions, wat.positions)
# determine which distances <= radius
		n_contacts = contacts.contact_matrix(dist, 5).sum()
		timeseries.append(n_contacts)

	print(timeseries)

# part re the score
	score_per_frame = []
	for index, value in enumerate(timeseries):
		print("DEBUG:", index, value, numHB_frame[index])
		if numHB_frame[index] > 0: # there are HB fot this frame
			HBscore = value**(1/numHB_frame[index])
		else: #  no HB fot this frame
			HBscore = 10
# some logging for debugging etc
		score_per_frame.append(HBscore)
		logScores_detailed(index, value, numHB_frame[index],HBscore)
	return score_per_frame, score_per_frame[-1]




