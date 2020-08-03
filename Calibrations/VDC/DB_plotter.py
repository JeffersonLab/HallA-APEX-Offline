#!/usr/bin/python


#######################################################
# Python Script to plot vdc t0 offsets
# with chosen arm and run number
#  John Williamson
#  30/7/2020
#######################################################


import sys

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import numpy as np


# function designed to read in DB with t0 values into
# need to select run number, arm and version (slope or gaus calc of t0?)
def read_t0s(run, arm, version):

    # vairable to hold lists of t0s from all planes
    #t0s = []
    t0s = {'u1':[],'u2':[],'v1':[],'v2':[]}

    # plane switches
    plane_switch = None
    
    
    
    
    with open(f'DB/db_{arm}_{version}.vdc.{run}.dat','r') as file:
        for line in file:

            if line.startswith(f'{arm}.vdc.u1'):
                plane_switch = 'u1'
                continue
            if line.startswith(f'{arm}.vdc.u2'):
                plane_switch = 'u2'
                continue
            if line.startswith(f'{arm}.vdc.v1'):
                plane_switch = 'v1'
                continue
            if line.startswith(f'{arm}.vdc.v2'):
                plane_switch = 'v2'
                continue

            line_t0s = line.split()

            for t0 in line_t0s:
                t0s[plane_switch].append(float(t0))
                
            #print(plane_switch)
                
        
    return t0s
        
        
    
    

def DB_plotter(run, arm = 'L'):

    print(f'{arm}-arm and run {run}')
    
    t0s_gaus = read_t0s(run,arm,'gaus')

    t0s_slope = read_t0s(run,arm,'slope')


    # creat list of vdc planes

    planes = ['u1','u2','v1','v2']

    figs = {}
    axs = {}
    
    for idx,plane in enumerate(planes):
        print(plane)

        figs[idx] = plt.figure()
        axs[idx] = figs[idx].add_subplot(111)

        axs[idx].set_title(f'{arm}-arm {plane} t0 calibration (run {run})')
    
        axs[idx].plot(t0s_slope[plane],'r+', label='slope')
        axs[idx].plot(t0s_gaus[plane],'b+', label='gaus')

        axs[idx].legend(loc=4) # loc=4 gives lower right placement of legend
 
        major_ticks = np.arange(8,368,16)
        axs[idx].set_xticks(major_ticks, minor=True)

        for tic in axs[idx].xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
        
        for tic in axs[idx].xaxis.get_minor_ticks():
            tic.tick1On = tic.tick2On = False

        axs[idx].set_xlabel("Wire number")
        axs[idx].set_ylabel("t0 channel number")
    
        axs[idx].xaxis.grid(True, which='minor')
        #    ax.xaxis.grid(which='major', alpha = 16)


        figs[idx].savefig(f"plots/comparison/{arm}_{run}_{plane}_t0s_comp.pdf")
        print(f'type of fig is {type(figs[idx])}')
        

    print(f'type of figs = {type(figs)}')

    with PdfPages(f"plots/comparison/{arm}_{run}_t0s_comp.pdf") as pdf:
        for fig in figs:
            pdf.savefig(figs[fig])
            print(f'type of fig is {type(fig)}')
            print(f'fig = {figs[fig]}')
            #pdf.savefig(fig)


def main():
    run = sys.argv[1]
    arm = sys.argv[2]

    DB_plotter(run,arm)


main()
