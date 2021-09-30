#!/usr/bin/env python3
## =============================================================================
#$ -N recons-GRAND
#$ -P P_trend
#$ -j y
#$ -cwd
#$ -notify
#$ -l h_rt=48:00:00
#$ -l s_rss=1.G
#$ -l s_fsize=0.1G
#$ -l sps=1
## =============================================================================
import os
import shutil
import sys
import time
import numpy as np
sys.path.insert(0, '/sps/hep/trend/vdecoene/recons_cca')
import mod_recons_tools as recons
import hdf5fileinout as hdf5io
import mod_fun as modfun
########################################################################################
# Reconstruction modules direcory
ReconsDir = "/sps/hep/trend/vdecoene/WavefrontStudy/"
# TREND direcory
TRENDDir = "/sps/hep/trend/vdecoene/WavefrontStudy/TREND_wavefront/"
# Output direcory for reconstruction results.
#OutDir = "/sps/hep/trend/vdecoene/recons_cca/recons_outputs/"
OutDir = "/sps/hep/trend/vdecoene/WavefrontStudy/recons_wavefront/"
########################################################################################
# # Configure the environment.
TAG = os.getenv("JOB_ID")
TMPDIR = '/tmp/vdecoene/'
###############################################################
###############################################################

###############################################################
###############################################################

if __name__ == "__main__":

    #create the reconstruction directory in the USER space as "OutDir + run#.txt" (-.txt)
    reconsdir = OutDir + sys.argv[1].split(".txt")[0]
    if not os.path.exists(reconsdir):
    	os.makedirs(reconsdir)
    #############################
    #############################

    SimulationList = []
    #Retrieve the run#.txt containing the list of path to reconstructions
    for line in open(ReconsDir+sys.argv[1], 'r'):
        SimulationList.append(line.replace('\n', ''))                               #Simulation list to reconstruct argv[1] -> filename of file containing simulation *.hdf5 files
    AntennasThreshold = float(sys.argv[2])                                          #Amplitude threshold for T1 trigger at antenna level
    AmplitudeThreshold = float(sys.argv[3])                                         #Antennas number for T2 trigger at the array level
    TablesFlag = int(sys.argv[4])                                                   #Creates coord_antennas.txt and Rec_coicntable.txt tables
    TREND_flag = int(sys.argv[5])                                                   #Run the reconstruction with TREND (Tools for Reconstruction in Energy, Nature and Direction)

    #############################
    #SELECTION
    #############################

    if len(SimulationList) == 0 :
        print("No simulations found -> system exit")
        sys.exit()

    if TablesFlag:

        #Retrieve all simulation parameters for reconstruction:
        #   1) ShowersParametersAll (shower geometry, energy, xmax, etc...) -> for analysis
        #   3) AntennasParametersAll (IDs, X, Y, Z, Peaktimes, PeakAmplitudes) -> for reconstruction methods
        ShowersParametersAll, AntennasParametersAll = recons.GetAllReconstructionParameters(SimulationList, AntennasThreshold, AmplitudeThreshold, flag_P2PInfo=False, Verbatim=True)

        if len(ShowersParametersAll) > 0:
            #Write txt tables for reconstruction with TREND (Tools for Reconstruction in Energy, Nature and Direction) software
            a = recons.WriteAntennaPositionsTable(reconsdir+"/coord_antennas.txt", AntennasParametersAll)
            b = recons.WriteReconsTable(reconsdir+"/Rec_coinctable.txt", AntennasParametersAll)
            c = recons.WriteToTxt(reconsdir+"/shower_parameters_"+str(AmplitudeThreshold)+"_"+str(AntennasThreshold)+".txt", ShowersParametersAll)
        else:
            print("No simulations suitable for reconstruction -> system exit")
            sys.exit()

    else:
        print("Using old coord_antennas.txt and Rec_coinctable.txt tables for next steps")

    #############################
    #RECONSTRUCTION
    #############################

    if TREND_flag:

        #Run TREND (Tools for Reconstruction in Energy, Nature and Direction)
        # 1) Plane Wave Front Reconstruction
        os.system(TRENDDir+"bin/recons 0 "+reconsdir)
        # 2) Spherical Reconstruction
        os.system(TRENDDir+"bin/recons 1 "+reconsdir)
        # 3) LDF Reconstruction
        #os.system(TRENDDir+"bin/recons 2 "+reconsdir)

    else:
        print("Using old reconstruction results for next steps")

    print("Reconstruction done")
