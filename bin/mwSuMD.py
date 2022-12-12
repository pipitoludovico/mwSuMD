#!/usr/bin/env python
import sys
import os
import __main__
#from moleculekit.molecule import Molecule

PATH="/".join(sys.argv[0].split("/")[:-1])
print(PATH)
LIBPATH=PATH.replace("bin", "lib")
sys.path.append(LIBPATH)
PARPATH=PATH.replace("bin", "toppar_c36_jul20")
print(LIBPATH)
print(PATH)
print(sys.path)
print(PARPATH)

from SuMD_1_CV import *
from SuMD_2_CV import *



#########################################	
def getpid():
	mypid = os.getpid()
	pidFile = open(".mypid","w")
	pidFile.write(str(mypid))
	pidFile.close()
	os.system('date >> .dates.txt')
	
	
#########################
def check_kill_process():
	with open(".mypid","r") as f:
		for line in f:
			pid = line.split()[0].strip()
	os.system('kill %s' %pid)
	quit()	
	
	
##########################
def logsettings(par, selection_list):
	logF = open("settings.txt", "w")	
	for sel in 	selection_list:
		print(sel,selection_list.index(sel) )
		logF.write('Metric_%s	%s\n' %(str(selection_list.index(sel)), sel))
	print(par, file=logF)
	logF.close()


##########################################
def inputFile():

	
	try:
		input_file = sys.argv[1]
	except:
		print('Input file for SuMD simulation required')
		quit()

	if input_file == 'kill':
		check_kill_process()


	par = {}
	selection_list = []

	with open(input_file, "r") as infile:
		for line in infile:	
			if line.startswith('#'):
				continue 


			if line.startswith('MDEngine'):
				par['MDEngine'] = line.split('=')[1].strip()
				
			if line.startswith('GROMACS_MDP'):
				par['GROMACS_MDP'] = line.split('=')[1].strip()
				if len(par['GROMACS_MDP']) > 0:
					par['GROMACS_MDP'] = os.path.abspath(par['GROMACS_MDP'])					
				else:
					par['GROMACS_MDP'] = None
				
			if line.startswith('GROMACS_TOP'):
				par['GROMACS_TOP'] = line.split('=')[1].strip()
				if len(par['GROMACS_TOP']) > 0:
					par['GROMACS_TOP'] = os.path.abspath(par['GROMACS_TOP'])					
				else:
					par['GROMACS_TOP'] = None				
				
			if line.startswith('Forcefield'):
				par['Forcefield'] = line.split('=')[1].strip()

			if line.startswith('Restart'):
				par['Restart'] = line.split('=')[1].strip()
				
			if line.startswith('Wrap'):
				par['Wrap'] = line.split('=')[1].strip()

			if line.startswith('Topology'):
				par['Topology'] = line.split('=')[1].strip()
				par['Topology'] = os.path.abspath('%s' %par['Topology']+'.pdb')
				pathname, extension = os.path.splitext(par['Topology'])
				par['Topology'] = pathname
				
				
			if line.startswith('Output'):
				par['Output'] = line.split('=')[1].strip()
				
			if line.startswith('NumberCV'):
				par['NumberCV'] = line.split('=')[1].strip()
				
			if line.startswith('REFERENCE'):
				par['REFERENCE'] = line.split('=')[1].strip()
				if len(par['REFERENCE']) > 0:
					par['REFERENCE'] = os.path.abspath(par['REFERENCE'])					
				else:
					par['REFERENCE'] = None

			if line.startswith('PLUMED'):
				par['PLUMED'] = line.split('=')[1].strip()
				if len(par['PLUMED']) > 0:
					par['PLUMED'] = os.path.abspath(par['PLUMED'])
				else:
					par['PLUMED'] = None

			if line.startswith('ligand_HB'):
				par['ligand_HB'] = line.split('=')[1].strip()
				if len(par['ligand_HB']) == 0:
					par['ligand_HB'] = None

			if line.startswith('Metric_1'):
				par['Metric_1'] = line.split('=')[1].strip()
				
			if line.startswith('Cutoff_1'):
				par['Cutoff_1'] = line.split('=')[1].strip()
				
			if line.startswith('Slope'):
				par['Slope'] = line.split('=')[1].strip()
			
			if line.startswith('Metric_2'):
				par['Metric_2'] = line.split('=')[1].strip()
				if len(par['Metric_2']) == 0:
					par['Metric_2']	 = None	

			if line.startswith('Cutoff_2'):
				par['Cutoff_2'] = line.split('=')[1].strip()
				if len(par['Cutoff_2']) == 0:
					par['Cutoff_2']	 = None	
					
			if line.startswith('Walkers'):
				par['Walkers'] = line.split('=')[1].strip()

			if line.startswith('Timewindow'):
				par['Timewindow'] = line.split('=')[1].strip()
				
			if line.startswith('Timestep'):
				par['Timestep'] = line.split('=')[1].strip()

			if line.startswith('Savefreq'):
				par['Savefreq'] = line.split('=')[1].strip()

			if line.startswith('GPU_ID'):
				par['GPU_ID'] = line.split('=')[1].strip()

			if line.startswith('Parameters'):
				if len(line.split('=')[1].strip()) > 0:
					par['Parameters'] = line.split('=')[1].split(';')
					parfiles = []
					for i in par['Parameters']:
						i = i.strip()
						parfiles.append(i)
					par['Parameters'] = parfiles
					for i in range(len(par['Parameters'])):
						par['Parameters'][i] = os.path.abspath(par['Parameters'][i])
					
				else:
					par['Parameters'] = None


			if line.startswith('Sel_'):
				selection_list.append(line.split('=')[1].strip())
				
			if line.startswith('Transition_1'):
				par['Transition_1'] = line.split('=')[1].strip().lower()
				
			if line.startswith('Transition_2'):
				par['Transition_2'] = line.split('=')[1].strip().lower()
				if len(par['Transition_2']) == 0:
					par['Transition_2']	 = None	
				

				
	return par, selection_list



#######

def main():
	getpid()
	par, selection_list = inputFile()
	logsettings(par, selection_list) 


	if par['NumberCV'] == '1':
		SuMD_1_CV(par, selection_list, PARPATH)
		
	elif par['NumberCV'] == '2':
		SuMD_2_CV(par, selection_list, PARPATH)
	
	else:
		print('Check the number of metrics in the setting input file')
		quit()

	
		
if __main__:
	main()		
			



