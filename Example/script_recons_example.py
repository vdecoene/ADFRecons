import os
import sys
import glob
import numpy as np
#sys.path.insert(0, 'path_to_modules') #Make sure to access these modules
import hdf5fileinout as hdf5io
import mod_recons_tools as recons
################################################################################
#Exemple script for reconstruction procedure
################################################################################
#Paths
TREND_path = "" #Path to the TREND code
Script_path = os.getcwd()
#############################
#############################

SimulationList = glob.glob(sys.argv[1]  +'*/*.hdf5')                            #Simulation list to reconstruct argv[1] -> directory path
AntennasThreshold = float(sys.argv[2])                                          #Antennas number for T2 trigger at the array level
AmplitudeThreshold = float(sys.argv[3])                                         #Amplitude threshold for T1 trigger at antenna level
TablesFlag = int(sys.argv[4])                                                   #Creates coord_antennas.txt and Rec_coicntable.txt tables
TREND_flag = int(sys.argv[5])                                                   #Run the reconstruction with TREND (Tools for Reconstruction in Energy, Nature and Direction)
Analysis_flag = int(sys.argv[6])                                                #Run the analysis

#############################
#SELECTION
#############################

if len(SimulationList) == 0 :
    print("No simulations found -> system exit")
    sys.exit()

##Remove the *.NoTraces.hdf5 files
CleanSimList = []
for sim in SimulationList:
    if sim.find('NoTraces') == -1:
        CleanSimList.append(sim)
SimulationList = CleanSimList

if TablesFlag:

    #Retrieve all simulation parameters for reconstruction:
    #   1) ShowersParametersAll (shower geometry, energy, xmax, etc...) -> for analysis
    #   3) AntennasParametersAll (IDs, X, Y, Z, Peaktimes, PeakAmplitudes) -> for reconstruction methods
    ShowersParametersAll, AntennasParametersAll = recons.GetAllReconstructionParameters(SimulationList, AntennasThreshold, AmplitudeThreshold, flag_P2PInfo=False, Verbatim=True)

    if len(ShowersParametersAll) > 0:
        #Write txt tables for reconstruction with TREND (Tools for Reconstruction in Energy, Nature and Direction) software
        a = recons.WriteAntennaPositionsTable("coord_antennas.txt", AntennasParametersAll)
        b = recons.WriteReconsTable("Rec_coinctable.txt", AntennasParametersAll)
        c = recons.WriteToTxt("shower_parameters_"+str(AmplitudeThreshold)+"_"+str(AntennasThreshold)+".txt", ShowersParametersAll)
    else:
        print("No simulations suitable for reconstruction -> system exit")
        sys.exit()

    #Move all files to the TREND direcory
    os.system("cp coord_antennas.txt "+TREND_path+" && cp Rec_coinctable.txt "+TREND_path)

else:
    print("Using old coord_antennas.txt and Rec_coinctable.txt tables for next steps")

#############################
#RECONSTRUCTION
#############################

if TREND_flag:

    #Run TREND (Tools for Reconstruction in Energy, Nature and Direction)
    # 1) Plane Wave Front Reconstruction
    os.system(TREND_path+"bin/recons 0 "+Script_path)
    # 2) Spherical Reconstruction
    os.system(TREND_path+"bin/recons 1 "+Script_path)
    # 3) LDF Reconstruction
    os.system(TREND_path+"bin/recons 2 "+Script_path)

else:
    print("Using old reconstruction results for next steps")

#############################
#ANALYSIS
#############################

if Analysis_flag:

    #Run the analysis

    if TablesFlag==0 :

        # 0) Load shower parameters
        EventName, Azimuth, Zenith, Energy, Primary, XmaxDistance, SlantXmax, x_Xmax, y_Xmax, z_Xmax, AntennasNumber = np.loadtxt("./shower_parameters_"+str(AmplitudeThreshold)+"_"+str(AntennasThreshold)+".txt").T
    else :
        EventName, Azimuth, Zenith, Energy, Primary, XmaxDistance, SlantXmax, x_Xmax, y_Xmax, z_Xmax, AntennasNumber = ShowersParametersAll.T
        if len(EventName)==1:
            EventName = EventName[0]; Azimuth = Azimuth[0]; Zenith = Zenith[0]; Energy = Energy[0]; Primary = Primary[0]; XmaxDistance = XmaxDistance[0]; SlantXmax = SlantXmax[0]; x_Xmax =x_Xmax[0]; y_Xmax=y_Xmax[0]; z_Xmax=z_Xmax[0]; AntennasNumber = AntennasNumber[0]

    # 1) Plane Wave Front Analysis
    IDsRec, _, ZenithRec, _, AzimuthRec, _, Chi2, Signif = np.loadtxt("./Rec_plane_wave_recons.txt").T

    AzimErrors, ZenErrors = recons.ComputeAngularErrors(AzimuthRec, ZenithRec, Azimuth, Zenith)                 #Compute angles errors
    AngularDistances = recons.ComputeAngularDistance(AzimuthRec, ZenithRec, Azimuth, Zenith)                    #Compute angular errors projected on sphere

    print()
    print("****Plane Wave Reconstruction****")
    print("Mean angular error = ", np.mean(AngularDistances), " STD = ", np.std(AngularDistances))
    PlaneRecStats = np.array([EventName, Azimuth, Zenith, Energy, Primary, XmaxDistance, SlantXmax, x_Xmax, y_Xmax, z_Xmax , AntennasNumber, ZenithRec, AzimuthRec, Chi2, Signif, AzimErrors, ZenErrors, AngularDistances]).T
    e = recons.WriteToTxt("plane_wave_recons_stats.txt", PlaneRecStats)

    # 3) LDF Analysis
    IDsRec, _, ZenithRec, _, AzimuthRec, _, Chi2, _, WidthRec, AmpRec = np.loadtxt("./Rec_adf_recons.txt").T

    AzimErrors, ZenErrors = recons.ComputeAngularErrors(AzimuthRec, ZenithRec, Azimuth, Zenith)                 #Compute angles errors
    AngularDistances = recons.ComputeAngularDistance(AzimuthRec, ZenithRec, Azimuth, Zenith)                    #Compute angular errors projected on sphere

    print()
    print("****ADF Reconstruction****")
    print("Mean angular error = ", np.mean(AngularDistances), " STD = ", np.std(AngularDistances))

    CerenkovRecStats = np.array([EventName, Azimuth, Zenith, Energy, Primary, XmaxDistance, SlantXmax, x_Xmax, y_Xmax, z_Xmax , AntennasNumber, ZenithRec, AzimuthRec, Chi2, WidthRec, AmpRec, AzimErrors, ZenErrors, AngularDistances]).T
    g = recons.WriteToTxt("adf_recons_stats.txt", CerenkovRecStats)

    # 2) Spherical Analysis
    IDsRec, _, Chi2, _, XSourceRec, YSourceRec, ZSourceRec, TSourceRec = np.loadtxt("./Rec_sphere_wave_recons.txt").T

    if np.isscalar(Zenith) :
        InjectionHeight = 1.e5                                                  #upper atmosphere 100km
        ShowerCoreHeight = 1086.#2900.                                          #depends on simulation...
    else:
        InjectionHeight = np.repeat(1.e5, len(Zenith))
        ShowerCoreHeight = np.repeat(1086., len(Zenith))#np.repeat(2900., len(Zenith))

    XError = x_Xmax - XSourceRec
    YError = y_Xmax - YSourceRec
    ZError = z_Xmax - ZSourceRec

    print()
    print("****Spherical Wave Reconstruction****")
    print("Mean X error = ", np.mean(XError), " STD = ", np.std(XError))
    print("Mean Y error = ", np.mean(YError), " STD = ", np.std(YError))
    print("Mean Z error = ", np.mean(ZError), " STD = ", np.std(ZError))

    LongitudinalError, LateralError = recons.ComputeSourceError_Long_Lat(Azimuth, Zenith, x_Xmax, y_Xmax, z_Xmax, XSourceRec, YSourceRec, ZSourceRec)

    print("Mean longitudinal error = ", np.mean(LongitudinalError), " STD = ", np.std(LongitudinalError))
    print("Mean lateral error = ", np.mean(LateralError), " STD = ", np.std(LateralError))

    GrammageRecons, GrammageError, LongitudinalDistance_Xmax, LongitudinalDistance_Source = recons.ComputeSourceErrorGrammage_alternative_method(Azimuth, Zenith, x_Xmax, y_Xmax, z_Xmax, AzimuthRec, ZenithRec, XSourceRec, YSourceRec, ZSourceRec, InjectionHeight, ShowerCoreHeight, XmaxDistance)
    print("Mean Grammage error = ", np.mean(GrammageError), " STD = ", np.std(GrammageError))

    SphereRecStats = np.array([EventName, Azimuth, Zenith, Energy, Primary, XmaxDistance, SlantXmax, x_Xmax, y_Xmax, z_Xmax , AntennasNumber, ZenithRec, AzimuthRec, Chi2, XSourceRec, YSourceRec, ZSourceRec, TSourceRec, GrammageRecons, XError, YError, ZError, LongitudinalError, LateralError, GrammageError, LongitudinalDistance_Xmax, LongitudinalDistance_Source]).T
    f = recons.WriteToTxt("Sphere_wave_recons_stats.txt", SphereRecStats)


else:
    print("No analysis performed")
