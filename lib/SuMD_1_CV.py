from miscellaneous import *
from metrics import *


#################### Compute single metric and return list of values per frame
def metricCompute(par, selection_list, n, mothfolder):
	
	os.chdir('tmp')
	walkers_metrics = []
	for r in range(1, int(par['Walkers'])+1):
		os.chdir('walker_%s' %str(r))
		if (glob.glob('*.xtc') and glob.glob('output.*')) or (glob.glob('*.xtc') and glob.glob('*.gro')):
    # do stuff if no file ending in .true exists
# First metric is a distance 
			if par['Metric_1'] == 'Distance':
# get the score of the walker
				distMetric = distance(par, selection_list[0], selection_list[1], n, mothfolder)
# add it to the list 
				walkers_metrics.append(distMetric)
				os.chdir(mothfolder+'/tmp')
			
			if par['Metric_1'] == 'Contacts':
				distMetric = contacts(par, selection_list[0], selection_list[1], n, mothfolder)
				walkers_metrics.append(distMetric)		
				os.chdir(mothfolder+'/tmp')
			
			if par['Metric_1'] == 'RMSD':
				distMetric = RMSD(par, selection_list[0], selection_list[1], n, mothfolder)
				#print(distMetric)
				walkers_metrics.append(distMetric)		
				os.chdir(mothfolder+'/tmp')			
		else:
			os.chdir('..')
			walkers_metrics.append(None)	
	os.chdir(mothfolder)
	return walkers_metrics
	
###################	Compute score on list of single metric values
def bestWalker(par, walkers_metrics):
	if all(v is None for v in walkers_metrics):
		return (0,0)	


	if par['Transition_1'] == 'positive': # we want the metric to increase
		value = max([i for i in walkers_metrics if i is not None])
		#max_value = max(walkers_metrics) # so we take the maximum value
		index = walkers_metrics.index(value)
		
	elif par['Transition_1'] == 'negative': # we want the metric to dencrease
		value = min([i for i in walkers_metrics if i is not None])
		#max_value = min(walkers_metrics) # so we take the minimum value
		index = walkers_metrics.index(value)
	
	return (index +1,value)
					

###############################################################
###############################################################
def SuMD_1_CV(par, selection_list, PARPATH):
	
	
	mothfolder = os.getcwd()

	if par['Restart'] == 'YES':
		n = restart(mothfolder)
		c = 0

	if par['Restart'] == 'NO':
		n = 0
		c = 0
		make_folders()
		
	if par['Transition_1'] == 'positive':

		if par['Slope'] == 'NO':
			max_value = 0
			while max_value < float(par['Cutoff_1']):
			# run acemd
				runMD(par, n, mothfolder, PARPATH)
			# get list of single metric values for each walker
				walkers_metrics = metricCompute(par, selection_list, n, mothfolder)
				print(walkers_metrics)
			# compute the score and decide the best walker
				best_walker, max_value = bestWalker(par, walkers_metrics)
				if best_walker == 0:
					continue		
			# write log file
				logStep(par, n, walkers_metrics, best_walker, mothfolder)
			# save trj and restart files
				saveStep(par, best_walker, n, mothfolder)
				n += 1

		elif par['Slope'] == 'YES':
			max_cycles = 1/(int(par['Timewindow'])/10**5) # run for 1 microsecond and then stop
			while c < max_cycles:
			# run acemd
				runMD(par, n, mothfolder, PARPATH)
			# get list of single metric values for each walker
				walkers_metrics = metricCompute(par, selection_list, n, mothfolder)
				print(walkers_metrics)
			# compute the score and decide the best walker
				best_walker, max_value = bestWalker(par, walkers_metrics)
				if best_walker == 0:
					continue
			# write log file
				logStep(par, n, walkers_metrics, best_walker, mothfolder)
			# save trj and restart files
				saveStep(par, best_walker, n, mothfolder)
				n += 1			
				

	if par['Transition_1'] == 'negative':
		if par['Slope'] == 'NO':
			max_value = 1000000
			while max_value > float(par['Cutoff_1']):
				runMD(par, n, mothfolder, PARPATH)
				walkers_metrics = metricCompute(par, selection_list, n, mothfolder)
				best_walker, max_value = bestWalker(par, walkers_metrics)
				logStep(par, n, walkers_metrics, best_walker, mothfolder)
				saveStep(par, best_walker, n, mothfolder)
				n += 1
	
		elif par['Slope'] == 'YES':
			max_cycles = 1/(int(par['Timewindow'])/10**5) # run for 1 microsecond and then stop
			while c < max_cycles:
				runMD(par, n, mothfolder, PARPATH)
				walkers_metrics = metricCompute(par, selection_list, n, mothfolder)
				best_walker, max_value = bestWalker(par, walkers_metrics)
				logStep(par, n, walkers_metrics, best_walker, mothfolder)
				saveStep(par, best_walker, n, mothfolder)
				n += 1
	
