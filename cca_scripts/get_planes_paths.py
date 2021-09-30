import os
import glob
import sys
import numpy as np
############################################################
def WriteToTxt(file_name_, InputsArray_):
    '''
    Write to txt file any array
    '''
    # TODO: handle errors
    _file = open(file_name_, 'w')
    np.savetxt(_file, InputsArray_, fmt='%s')
    _file.close()

    return 0

############################################################

SimFolderPath = sys.argv[1]

#Matias Librairy is organised as SimFolder/hdf5_file -> we only want hdf5_files
#We want to gather simulations with same geometry and the the two farthest planes i.e. 133km and 200km
#File name: Stshp_Stshp_Iron_3.98_87.1_0.0_1_200000.0

PrimaryList = ["Proton", "Iron"] #To Complete
EnergyList = ["3.98", "3.16", "2.51", "2.0", "1.58", "1.26", "1.0", "0.794",  "0.158", "0.126", "0.0631", "0.0501", "0.0398", "0.0316", "0.0251", "0.02"]
ZenithList = ["71.6", "74.8", "77.4", "79.5", "81.3", "82.7", "83.9", "85.0", "85.8", "86.5", "87.1"] #To Complete
AzimuthList = ["0.0", "90.0", "180.0"] #To Complete
PlanList = ["132747.0", "200000.0"] #To Complete

SimList = []
#Ugly but fuck you
for prim in PrimaryList:
    for en in EnergyList:
        for zen in ZenithList:
            for azim in AzimuthList:
                    SimList.append("Stshp_Stshp_"+prim+"_"+en+"_"+zen+"_"+azim+"_1_")
i = 0
for sim in SimList:
    SimPaths = [SimFolderPath+sim+PlanList[0]+"/"+sim+PlanList[0]+".hdf5", SimFolderPath+sim+PlanList[1]+"/"+sim+PlanList[1]+".hdf5"]
    print(SimPaths)
    WriteToTxt("run"+str(i)+".txt", np.array(SimPaths))
    i+=1
