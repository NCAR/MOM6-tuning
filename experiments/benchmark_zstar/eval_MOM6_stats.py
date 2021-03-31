from cheyenne.read_amplxe_results import read_amplxe_results
import sys
import os
import time
import xarray as xr
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def get_MOM_continuity_PPM_time(logPath):

    totalTime = 0.0

    routines = read_amplxe_results(os.path.join(logPath,"vtune_report.txt"))

    for r in [routines[r] for r in routines]:
        if r.name.startswith("mom_continuity_ppm"):
            totalTime += float(r.self)

    return totalTime


def check_consistency(ds_list_ensemble, ds_experiment):

    EPS = 3.0

    fig, axs = plt.subplots(3,figsize=(10, 6))
    fig.suptitle('Statistical consistency check')
    #axs[0].plot(x, y)
    #axs[1].plot(x, -y)
    
    min_consistency_score = 0.7
    N = len(ds_list_ensemble)
    PASSED = True
    
    # total mass change
    da_ensemble = [ds['Mass_anom'] for ds in ds_list_ensemble]
    da_exp = ds_experiment['Mass_anom']
    for da in da_ensemble:
        axs[0].plot(da, color='gray', alpha=.4)
    #mean
    da_mean = sum(da_ensemble)/N
    axs[0].plot(da_mean, color='blue', alpha=.6)
    #std
    da_std = np.sqrt( sum([(da-da_mean)**2 for da in da_ensemble]) / N )   
    axs[0].plot(da_mean+EPS*da_std, color='orange')
    axs[0].plot(da_mean-EPS*da_std, color='orange')
    axs[0].set_title(da.long_name + " ("+da.units+")")
    #test
    da_test = np.logical_and(da_exp<(da_mean+EPS*da_std), da_exp>(da_mean-EPS*da_std))
    consistency_score = sum(da_test)/len(da_test)
    if consistency_score>min_consistency_score:
        axs[0].plot(da_exp, color='green')
    else:
        axs[0].plot(da_exp, color='red')
        PASSED = False

    # total mass change
    da_ensemble = [ds['Salt_anom'] for ds in ds_list_ensemble]
    da_exp = ds_experiment['Salt_anom']
    for da in da_ensemble:
        axs[1].plot(da, color='gray', alpha=.4)
    #mean
    da_mean = sum(da_ensemble)/N
    axs[1].plot(da_mean, color='blue', alpha=.6)
    #std
    da_std = np.sqrt( sum([(da-da_mean)**2 for da in da_ensemble]) / N )   
    axs[1].plot(da_mean+EPS*da_std, color='orange')
    axs[1].plot(da_mean-EPS*da_std, color='orange')
    axs[1].set_title(da.long_name + " ("+da.units+")")
    #test
    da_test = np.logical_and(da_exp<(da_mean+EPS*da_std), da_exp>(da_mean-EPS*da_std))
    consistency_score = sum(da_test)/len(da_test)
    if consistency_score>min_consistency_score:
        axs[1].plot(da_exp, color='green')
    else:
        axs[1].plot(da_exp, color='red')
        PASSED = False
            
    # total mass change
    da_ensemble = [ds['Heat_anom'] for ds in ds_list_ensemble]
    da_exp = ds_experiment['Heat_anom']
    for da in da_ensemble:
        axs[2].plot(da, color='gray', alpha=.4)
    #mean
    da_mean = sum(da_ensemble)/N
    axs[2].plot(da_mean, color='blue', alpha=.6)
    #std
    da_std = np.sqrt( sum([(da-da_mean)**2 for da in da_ensemble]) / N )   
    axs[2].plot(da_mean+EPS*da_std, color='orange')
    axs[2].plot(da_mean-EPS*da_std, color='orange')
    axs[2].plot(da_exp, color='red')
    axs[2].set_title(da.long_name + " ("+da.units+")")
    #test
    da_test = np.logical_and(da_exp<(da_mean+EPS*da_std), da_exp>(da_mean-EPS*da_std))
    consistency_score = sum(da_test)/len(da_test)
    if consistency_score>min_consistency_score:
        axs[2].plot(da_exp, color='green')
    else:
        axs[2].plot(da_exp, color='red')
        PASSED = False
        
    timestr = time.strftime("%Y%m%d-%H%M%S")
    fig.savefig('check_{}_{}.png'.format('PASS' if PASSED else 'FAIL', timestr))
    return PASSED


if __name__ == "__main__":

    logPath = sys.argv[1]

    ensemble_dir = "/glade/work/jdvanover/quarantine/ECT_gen/MOM6-tuning/experiments/ensemble/benchmark_zstar.empty"
    nens = 25
    ens_ocnstats = []
    for i in range(1,nens+1):
        ens_ocnstats.append(xr.open_dataset(ensemble_dir+"/out_"+str(i)+"/ocean.stats.nc"))
    ds_experiment = xr.open_dataset( "./ocean.stats.nc")
    PASSED = check_consistency(ens_ocnstats, ds_experiment)

    if PASSED:
        clock = get_MOM_continuity_PPM_time(logPath)
    else:
        clock = -1.0
    print(clock)