import sys
import os

for i in range(int(sys.argv[1])):
    #only require 22 muV/m for threshold for the wavefront study
    print("qsub cca_recons.py run"+str(i)+".txt 5 44 1 1")
    os.system("qsub cca_recons.py run"+str(i)+".txt 5 44 1 1")
