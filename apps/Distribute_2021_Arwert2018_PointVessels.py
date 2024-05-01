import multiprocessing
import os
import subprocess
import itertools
import time
import random

import numpy as np

executable = '/mi/share/scratch/bull/ChasteStuff/chaste_build_Jan11/projects/JB_Chaste/apps/Exe_2021_Arwert2018_PointVessels'

chaste_test_dir = os.environ.get('CHASTE_TEST_OUTPUT')
path_to_output = os.path.join(chaste_test_dir,'Arwert2018','ParameterSweeps','II_PointVessels_iv')

if not(os.path.isfile(executable)):
    raise Exception('Could not find executable: ' + executable)

command_line_args = [' --ID ', ' --ACCD_T ',' --ACCD_S ',' --CS_M2CSF ',' --CS_M2CXCL12 ', ' --CS_TC2EGF ', ' --MPI ', ' --DTBV ',' --NOBV ',' --HMECSF ',' --MPOE ', ' --TGF_TFPS ',' --IT ']
params_list = ['simulation_id', 'averageCellCycleDuration_tumour','averageCellCycleDuration_stroma', 'chi_macrophageToCSF', 'chi_macrophageToCXCL12', 'chi_tumourCellToEGF', 'macrophagePhenotypeIncrementPerHour', 'distanceToBloodVessel','numberOfBloodVessels','halfMaximalExtravasationCsf1Conc','maximumProbOfExtravasationPerHour', 'tgfThresholdForPhenotypeSwitch', 'iterationNumber']

today = time.strftime('%Y-%m-%dT%H%M')

# helper function = choose a random combination of parameters to sample
def random_combination(iterable,r):
    "Random selection from itertools.combinations(iterable,r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n),r))
    out = tuple(pool[i] for i in indices)
    return out


# Param ranges (in lists, for itertools product)
accd_t_Range = [24,24]
accd_s_Range = [32,32]
cs_m2csf_Range = [0.1,5]
cs_m2cxcl12_Range = [0.1,5]
cs_tc2egf_Range = [0.1,5]
mpi_Range = [0.01,0.01]
dtbv_Range = [17,17]#[12,20]
nobv_Range = [35,36]#[20,50]
hmecsf_Range = [0.1,0.5]
mpoe_Range = [0.01,0.1]
tgf_tfps_Range = [0,2]

random.seed(23032021) # Seed for reproducibility
numSimulations = 48*10

def main():
    run_simulations()


# Create a list of commands and pass them to separate processes
def run_simulations():
    print('Running')
    # Make a list of calls to a Chaste executable
    command_list = []
    
    if not os.path.exists(path_to_output):
    	os.makedirs(path_to_output)
    
    params_file = open(path_to_output + '/params_file.csv', 'w')
    params_file.write(','.join(params_list) + '\n')

    base_command = 'nice ' + executable

    for i in range(numSimulations):
        param_set = [0 for v in range(13)]
        param_set[0] = i # 'simulation_id'
        param_set[1] = np.random.uniform(accd_t_Range[0],accd_t_Range[1])# 'averageCellCycleDuration_tumour'
        param_set[2] = np.random.uniform(accd_s_Range[0],accd_s_Range[1])#'averageCellCycleDuration_stroma'
        param_set[3] = np.random.uniform(cs_m2csf_Range[0],cs_m2csf_Range[1])# 'chi_macrophageToCSF'
        param_set[4] = np.random.uniform(cs_m2cxcl12_Range[0],cs_m2cxcl12_Range[1])# 'chi_macrophageToCXCL12'
        param_set[5] = np.random.uniform(cs_tc2egf_Range[0],cs_tc2egf_Range[1])# 'chi_tumourCellToEGF'
        param_set[6] = np.random.uniform(mpi_Range[0],mpi_Range[1])# 'macrophagePhenotypeIncrementPerHour'
        param_set[7] = np.random.uniform(dtbv_Range[0],dtbv_Range[1])# 'distanceToBloodVessel'
        param_set[8] = np.random.randint(nobv_Range[0],nobv_Range[1])# 'numberOfBloodVessels'
        param_set[9] = np.random.uniform(hmecsf_Range[0],hmecsf_Range[1])# 'halfMaximalExtravasationCsf1Conc'
        param_set[10] = np.random.uniform(mpoe_Range[0],mpoe_Range[1])# 'maximumProbOfExtravasationPerHour'
        param_set[11] = np.random.uniform(tgf_tfps_Range[0],tgf_tfps_Range[1])# 'tgfThresholdForPhenotypeSwitch'
        param_set[12] = 1# 'iterationNumber'


        param_set = tuple(param_set)
        idx = i

        params_file.write(",".join(map(str, param_set)) + '\n')
        command = base_command
        command += ' --ID ' + str(idx)

        for arg in range(1,len(param_set)):
            print(arg)
            print(command)
            command += command_line_args[arg] + str(param_set[arg])

        command_list.append(command)

    params_file.close()

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Py: Starting simulations with " + str(count) + " processes")

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Pass the list of bash commands to the pool
    # Wait at most one week
    pool.map(execute_command, command_list).get(6048000)
    #pool.map_async(execute_command, command_list).get(6048000)


# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    main()
