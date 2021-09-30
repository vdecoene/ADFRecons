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
NRuns = int(sys.argv[2])

#Matias Librairy is organised as SimFolder/hdf5_file -> we only want hdf5_files
SimList = glob.glob(SimFolderPath+"*/*.hdf5")
#Remove the *.NoTraces.hdf5 files
CleanSimList = []
for sim in SimList:
    if sim.find('NoTraces') == -1:
        CleanSimList.append(sim)

TotalSim = len(CleanSimList)
SubRuns = np.floor(TotalSim/NRuns)
LastSims = TotalSim - SubRuns*NRuns
print("TotalSim = ", TotalSim, " Nruns = ", NRuns, " Subruns = ", SubRuns, " LastSims = ", LastSims)
if (int(input("Proceed ? 0/1 : ")) == 0):
    print("Stopping")
    sys.exit()

for i in range(NRuns):
    print("Run #: ", i, " -> Simulations from ", i*int(SubRuns), " to ", (i+1)*int(SubRuns))
    InputsArray = CleanSimList[i*int(SubRuns):(i+1)*int(SubRuns)]
    #print(InputsArray)
    WriteToTxt("run"+str(i)+".txt", np.array(InputsArray))
