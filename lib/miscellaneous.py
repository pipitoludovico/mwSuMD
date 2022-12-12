import os
import csv 
#from moleculekit.molecule import Molecule
import glob


################### At what cycle number mwSuMD was stopped?
def restart(mothfolder):
	import glob
	
	n = len(glob.glob('./trajectories/*'))
	os.system("rm -r ./tmp/* ")

	logF = open(mothfolder+'/walkerSummary.log' , 'a')
	logF.write('##################################    Simulation Restarted    ###########################\n') 
	logF.close()

	return n

############################
def make_folders():
	os.system('mkdir trajectories')
	os.system('mkdir  restarts')
	os.system('mkdir tmp')

####################
def acemdInput(par, n, r, coor, xsc, vel, PARPATH):
	
	savefreq = int((int(par['Savefreq'])*1000)/int(par['Timestep']))
	
	eq_file = open('input_%s_%s.inp' %(str(n),str(r)), 'w')
	
	if par['Forcefield'] == 'CHARMM':
	
		txt = ['restart off\n', 'minimize        0\n', 'run            %sps\n' %par['Timewindow'], 'timeStep        %s\n' %par['Timestep'], 	'structure               %s.psf\n' %par['Topology'], \
		'coordinates             %s.pdb\n' %par['Topology'], 'temperature     310\n', 'PME             on\n', 'cutoff          9.0\n','switchDistance  7.5\n', \
		'thermostat      on\n', 'thermostatDamping       0.1\n', 'thermostatTemperature   310\n', 'barostat                off\n',  \
		'trajectoryFile          %s_%s.xtc\n' %(par['Output'], str(n)), 'trajectoryPeriod               %s\n' %str(savefreq), \
		'binCoordinates		%s.coor\n' %coor, \
		'extendedSystem          %s.xsc\n' %xsc, \
		'binVelocities           %s.vel\n' %vel, \
		'parameters              %s/par_all36_prot.prm\n' %PARPATH, \
		'parameters              %s/par_all36_lipid.prm\n' %PARPATH]
	
		for e in txt:
			eq_file.write(e)		

		if par['Parameters'] != None:
			for e in par['Parameters']:
				eq_file.write('parameters		%s\n' %e)
				
		if par['PLUMED'] != None:
			 eq_file.write('plumedFile		%s\n' %par['PLUMED'])
	
		eq_file.close()

	elif par['Forcefield'] == 'AMBER':
		txt = ['restart off\n', 'minimize        0\n', 'run            %sps\n' %par['Timewindow'], 'timeStep        %s\n' %par['Timestep'], 	'parmfile               %s.prmtop\n' %par['Topology'], \
		'coordinates             %s.pdb\n' %par['Topology'], 'temperature     310\n', 'PME             on\n', 'cutoff          9.0\n','switchDistance  7.5\n', \
		'thermostat      on\n', 'thermostatDamping       0.1\n', 'thermostatTemperature   310\n', 'barostat                off\n',  \
		'trajectoryFile          %s_%s.xtc\n' %(par['Output'], str(n)), 'trajectoryPeriod               %s\n' %str(savefreq), \
		'binCoordinates		%s.coor\n' %coor, \
		'extendedSystem          %s.xsc\n' %xsc, \
		'binVelocities           %s.vel\n' %vel]


		for e in txt:
			eq_file.write(e)
			
		if par['PLUMED'] != None:
			 eq_file.write('plumedFile		%s\n' %par['PLUMED'])
			 
		eq_file.close()

#################### Write acemd input file and run it
def runMD(par, n, mothfolder, PARPATH):

	os.chdir('tmp')
	
	if par['MDEngine'] == 'ACEMD':
		for r in range(1, int(par['Walkers'])+1):
			os.system('mkdir walker_%s' %str(r))
			if n == 0:
				os.chdir('walker_%s' %str(r))
				acemdInput(par, n, r, par['Topology'], par['Topology'], par['Topology'], PARPATH)
				print('acemd3 --device %s input_%s_%s.inp' %(par['GPU_ID'],str(n),str(r)))
				os.system('acemd3 --device %s input_%s_%s.inp' %(par['GPU_ID'],str(n),str(r)))
				wrap(par, n)		
				os.chdir(mothfolder+'/tmp')
			
			elif n > 0:
				os.chdir('walker_%s' %str(r))
				if par['PLUMED'] != None:
					os.system('cp %s/restarts/HILLS .' %mothfolder)
					os.system('cp %s/restarts/grid.dat .' %mothfolder)
					os.system('cp %s/restarts/COLVAR .' %mothfolder)
				acemdInput(par, n, r, '../../restarts/%s_%s' %(par['Output'], str(n-1)), '../../restarts/%s_%s' %(par['Output'],str(n-1)), '../../restarts/%s_%s' %(par['Output'],str(n-1)), PARPATH)
				print('acemd3 --device %s input_%s_%s.inp' %(par['GPU_ID'],str(n),str(r)))
				os.system('acemd3 --device %s input_%s_%s.inp' %(par['GPU_ID'],str(n),str(r)))	
				wrap(par, n)		
				os.chdir(mothfolder+'/tmp')


	if par['MDEngine'] == 'GROMACS':
		for r in range(1, int(par['Walkers'])+1):
			os.system('mkdir walker_%s' %str(r))
			if n == 0:
				os.chdir('walker_%s' %str(r))
				os.system('gmx grompp -f %s -c %s.gro -t %s.cpt -p %s -o %s_%s.tpr' %(par['GROMACS_MDP'], par['Topology'], par['Topology'], par['GROMACS_TOP'], par['Output'], str(n)))
				if par['PLUMED'] != None:
					os.system('gmx mdrun -plumed %s -deffnm %s_%s' %(par['PLUMED'], par['Output'],str(n)))
				else:
					os.system('gmx mdrun -deffnm %s_%s' %(par['Output'],str(n)))
				os.chdir(mothfolder+'/tmp')


			elif n > 0:
				os.chdir('walker_%s' %str(r))
				### add here commands to handle plumed if condition plumed satisfied
				os.system('gmx grompp -f %s -c ../../restarts/%s_%s.gro -t ../../restarts/%s_%s.cpt -p %s -o %s_%s.tpr' %(par['GROMACS_MDP'], par['Output'], str(n-1), par['Output'], str(n-1), par['GROMACS_TOP'], par['Output'], str(n)))
				if par['PLUMED'] != None:
					os.system('gmx mdrun -plumed %s -deffnm %s_%s' %(par['PLUMED'], par['Output'],str(n)))
				else:
					os.system('gmx mdrun -deffnm %s_%s' %(par['Output'],str(n)))
				os.chdir(mothfolder+'/tmp')				

	
	os.chdir(mothfolder)

########################### Eachtime window is wrapped after production before metrics computation
def wrap(par, n):
	from moleculekit.molecule import Molecule
	if par['Forcefield'] == 'CHARMM':
		psf = '%s.psf' %par['Topology']

	elif par['Forcefield'] == 'AMBER':
		psf = '%s.prmtop' %par['Topology']

	xtc = '%s_%s.xtc' %(par['Output'], str(n))

	if glob.glob(xtc):
		fname = os.path.splitext(os.path.basename(xtc))[0]

		mol = Molecule(psf)
		mol.read(xtc)
		mol.wrap()
		mol.wrap(par['Wrap'])
		mol.write('%s_wrapped.xtc' %fname)
	else:
		return

######################### Write log file with minimum info
def logStep(par, n, walkers_metrics, best_walker, mothfolder, walkers_metrics_1=[], walkers_metrics_2=[]):

	walkers_metrics.append(best_walker)
	
	with open(mothfolder+'/walkerSummary.log' , 'a') as logF:
		logF.write('#####     Step number %s     ######\n' %str(n))

		if par['NumberCV'] == '1':
		#with open(mothfolder+'/walkerSummary.log' , 'a') as logF:
			write = csv.writer(logF) 
			write.writerow(walkers_metrics)
		
		elif par['NumberCV'] == '2':		
		#with open(mothfolder+'/walkerSummary.log' , 'a') as logF:
			write = csv.writer(logF) 
			write.writerow(walkers_metrics_1)
			write.writerow(walkers_metrics_2)
			write.writerow(walkers_metrics)
		
		

################### Handle the restart files and the xtc storage
def saveStep(par, best_walker, n, mothfolder):
	
	for r in range(1, int(par['Walkers'])+1):
		if r == best_walker:
			os.chdir('tmp/walker_%s' %str(r))
			if par['MDEngine'] == 'ACEMD':
				os.system('cp %s_%s_wrapped.xtc %s/trajectories/%s_%s_wrapped.xtc '%(par['Output'], str(n), mothfolder, par['Output'], str(n)))
				os.system('mv  %s/restarts/%s_%s.coor %s/restarts/previous.coor' %(mothfolder, par['Output'], str(n-1), mothfolder))
				os.system('mv  %s/restarts/%s_%s.xsc %s/restarts/previous.xsc' %(mothfolder,par['Output'], str(n-1), mothfolder))
				os.system('mv  %s/restarts/%s_%s.vel %s/restarts/previous.vel' %(mothfolder,par['Output'], str(n-1),mothfolder))
				os.system('cp output.coor %s/restarts/%s_%s.coor '%(mothfolder, par['Output'], str(n)))
				os.system('cp output.xsc %s/restarts/%s_%s.xsc '%(mothfolder, par['Output'], str(n)))
				os.system('cp output.vel %s/restarts/%s_%s.vel '%(mothfolder, par['Output'], str(n)))

				
			

			elif par['MDEngine'] == 'GROMACS':
				os.system('cp %s_%s.xtc %s/trajectories/%s_%s.xtc '%(par['Output'], str(n), mothfolder, par['Output'], str(n)))
				os.system('mv  %s/restarts/%s_%s.gro %s/restarts/previous.gro' %(mothfolder, par['Output'], str(n-1), mothfolder))
				os.system('mv  %s/restarts/%s_%s.cpt %s/restarts/previous.cpt' %(mothfolder,par['Output'], str(n-1), mothfolder))
				os.system('cp %s_%s.cpt %s/restarts/%s_%s.cpt '%(par['Output'], str(n), mothfolder, par['Output'], str(n)))
				os.system('cp %s_%s.gro %s/restarts/%s_%s.gro '%(par['Output'], str(n), mothfolder, par['Output'], str(n)))


			if par['PLUMED'] != None:
				os.system('cp HILLS  %s/restarts/ '%mothfolder)
				os.system('cp COLVAR  %s/restarts/ '%mothfolder)
				os.system('cp grid.dat  %s/restarts/ '%mothfolder)

			os.chdir(mothfolder)
			os.system('rm -r tmp/walker_*')
	
	
	
	
	
