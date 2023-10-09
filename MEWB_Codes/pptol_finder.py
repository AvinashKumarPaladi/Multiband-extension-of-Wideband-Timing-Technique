import psrchive as pr  # psrchive
import numpy as np  # numpy
import ppspline as pps  # Import ppspline.py
import ppalign as ppa  # Import ppalign.py
import pptoas as ppt  # Import pptoas.py
from pplib import write_TOAs  # the function to write TOAs from pplib.py
import tempo_utils as tu  # Import tempo_utils to do the analysis
import matplotlib.pyplot as plt  # pyplot for simple plots
from sklearn.decomposition import PCA
import os
import sys
import time

psr_name = sys.argv[1] #'J2145' name of pulsar
portfile = sys.argv[2] #'files/J2145.200.ts.meta' #metafile containing template epoch.
bw = int(sys.argv[3])
neigen = int(sys.argv[4])
tol= [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

start = time.time()

d1=str(psr_name)+'output_tol_neig'+str(neigen)
try: os.mkdir(d1)
except OSError: pass
#d2='plots_'+tol
#try: os.mkdir(d2)
#except OSError: pass

#try:os.mkdir('tmp')
#except OSError: pass


#outfile = str(d1)+'/'+str(tol)+'_'+str(outfile1)

    
for i in tol:
    outfile = str(d1)+"/"+portfile+".n_"+str(neigen)+".T_"+str(i)+".spl"
    cmd1 = "ppspline.py -d "+portfile+" -o "+outfile+" -N prof -T "+str(i)+" -s -n "+str(neigen)+" -S 0 -t 2 --plots"
    print("COMMAND:"+cmd1)
    os.system(cmd1)
    #copy/rename the residual plots file so that we can retain all of them. feel free to modify the location. 
    cmd2="mv "+str(portfile)+".spl.resids.png "+str(d1)+"/"+str(portfile)+".n_"+str(neigen)+".T_"+str(i)+".spl.resids.png"
    print("COMMAND:"+str(cmd2))
    os.system(cmd2)
    #copy/rename the residual plots file so that we can retain all of them. feel free to modify the location. 
    cmd3="mv "+str(portfile)+".spl.freq.png "+str(d1)+"/"+str(portfile)+".n_"+str(neigen)+".T_"+str(i)+".spl.freq.png"
    print("COMMAND:"+str(cmd3))
    os.system(cmd3)
    #copy/rename the residual plots file so that we can retain all of them. feel free to modify the location. 
    cmd4="mv "+str(portfile)+".spl.proj.png "+str(d1)+"/"+str(portfile)+".n_"+str(neigen)+".T_"+str(i)+".spl.proj.png" 
    print("COMMAND:"+str(cmd4))
    os.system(cmd4)
    #copy/rename the residual plots file so that we can retain all of them. feel free to modify the location. 
    cmd5="mv "+str(portfile)+".spl.eigvec.png "+str(d1)+"/"+str(portfile)+".n_"+str(neigen)+".T_"+str(i)+".spl.eigvec.png"
    print("COMMAND:"+str(cmd5))
    os.system(cmd5)


end = time.time()
print('runtime: ', (end - start)/60, 'minutes')
