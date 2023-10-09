from __future__ import absolute_import, division, print_function
import numpy as np
import time
import os
import sys
import psrchive

#### Input files
Band3 = sys.argv[1]      # Meta file with path to the directory with Band3 fits files
Band5 = sys.argv[2]      # Meta file with path to the directory with Band5 fits files
mpath = sys.argv[3]      # Path to write the merged files with / in the end
pname = sys.argv[4]      # Pulsar Name
bw = int(sys.argv[5])    # Bandwidth of data

def merge_fits():
    
    temp_dir_b3 = pname+'.b3.'+str(bw)+'.temp'               # Creating a temporary directory to copy the band 3 data
    try : os.mkdir(temp_dir_b3)
    except OSError : pass

    temp_dir_b5 = pname+'.b5.'+str(bw)+'.temp'               # Creating a temporary directory to copy the band 5 data
    try : os.mkdir(temp_dir_b5)
    except OSError : pass
    
    file1 = open(Band3,'r')
    file2 = open(Band5,'r')

    file1_lines = file1.readlines()
    file1_lines = [i.strip() for i in file1_lines]
    file2_lines = file2.readlines()
    file2_lines = [i.strip() for i in file2_lines]
    
    for i in range(len(file1_lines)) : 
        cmd_b3 = 'cp -r '+file1_lines[i]+' '+temp_dir_b3       # Making a copy of the band 3 data to the temporary directory
        os.system(cmd_b3)
    
    for i in range(len(file2_lines)) : 
        cmd_b5 = 'cp -r '+file2_lines[i]+' '+temp_dir_b5       # Making a copy of the band 3 data to the temporary directory
        os.system(cmd_b5)
    
    metafile_b3 = pname+'.b3.'+str(bw)+'.temp.meta'                        # Creating new band 3 meta file
    metafile_b5 = pname+'.b5.'+str(bw)+'.temp.meta'                        # Creating new band 5 meta file
    
    cmd_meta_b3 = 'ls '+temp_dir_b3+'/* > '+metafile_b3
    os.system(cmd_meta_b3)
    
    cmd_meta_b5 = 'ls '+temp_dir_b5+'/* > '+metafile_b5
    os.system(cmd_meta_b5)
    
    file1.close()
    file2.close()
    
    
    # Correcting for be delay
    file1 = open(metafile_b3,'r')
    file2 = open(metafile_b5,'r')
    
    file1_lines = file1.readlines()
    file1_lines = [i.strip() for i in file1_lines]
    file2_lines = file2.readlines()
    file2_lines = [i.strip() for i in file2_lines]

    for i in range(len(file1_lines)):
        a = file1_lines[i]
        a1 = file1_lines[i].split('_')
        aa = float(a1[1])
        arch1 = psrchive.Archive_load(str(a))
        nsub1,npol1,nchan1,nbin1 = arch1.get_data().shape
        
        for j in range(len(file2_lines)):
            b = file2_lines[j]
            b1 = file2_lines[j].split('_')
            bb = float(b1[1])
            z = abs(aa-bb)          
            
            if z < 1 :
                arch2 = psrchive.Archive_load(b)
                nsub2,npol2,nchan2,nbin2 = arch2.get_data().shape
                print(f"merging {a} with {b}")
                be_b3 = arch1.get_backend_delay()
                be_b5 = arch2.get_backend_delay()
                Ps_b5 = np.array([sub.get_folding_period() for sub in arch2], dtype=np.double)[0]
                be_offset = (be_b3 - be_b5)/Ps_b5
                rotate = "pam -m -r "+str(be_offset%1)+" "+str(b)
                os.system(rotate)

    file1.close()
    file2.close()
    
    # Using the new meta files to render the original files unchanged
    file1 = open(metafile_b3,'r')
    file2 = open(metafile_b5,'r')

    file1_lines = file1.readlines()
    file1_lines = [i.strip() for i in file1_lines]
    file2_lines = file2.readlines()
    file2_lines = [i.strip() for i in file2_lines]
    
    for i in range(len(file1_lines)):
        
        a = file1_lines[i]
        a1 = file1_lines[i].split('_')
        aa = float(a1[1])                       # Extracting the MJD
        arch1 = psrchive.Archive_load(str(a))
        nsub1,npol1,nchan1,nbin1 = arch1.get_data().shape

        for j in range(len(file2_lines)):
            b = file2_lines[j]
            b1 = file2_lines[j].split('_')
            bb = float(b1[1])                   # Extracting the MJD
            z = abs(aa-bb)

            if z < 1: #checking if the epochs are within 24 hours
                print(z)                
                arch2 = psrchive.Archive_load(b)
                nsub2,npol2,nchan2,nbin2 = arch2.get_data().shape
                print(f'{a} and {b}')
                
                
                if nbin1 != nbin2:
                    c = min(nbin1,nbin2)
                    print(f"changing nbin to {c}")
                    arch1.bscrunch_to_nbin(c)
                    arch2.bscrunch_to_nbin(c)
                    arch1.unload()
                    arch2.unload()
                
                if bw == 200 : 
                    print(f"changing B3 nchan to {64} and B5 nchan to {32} for the %d"%bw+" MHz data...")
                    arch1.fscrunch_to_nchan(64)
                    arch2.fscrunch_to_nchan(32)
                    arch1.unload()
                    arch2.unload()
                
                elif bw == 100 : 
                    print(f"changing B3 nchan to {32} and B5 nchan to {16} for the %d"%bw+" MHz data...")
                    arch1.fscrunch_to_nchan(32)
                    arch2.fscrunch_to_nchan(16)
                    arch1.unload()
                    arch2.unload()

            else:
                continue
               
    file1.close()
    file2.close()
    
    # Using the new meta files to render the original files unchanged
    file1 = open(metafile_b3,'r')
    file2 = open(metafile_b5,'r')
    
    file1_lines = file1.readlines()
    file1_lines = [i.strip() for i in file1_lines]
    file2_lines = file2.readlines()
    file2_lines = [i.strip() for i in file2_lines]

    for i in range(len(file1_lines)):
        a = file1_lines[i]
        a1 = file1_lines[i].split('_')
        aa = float(a1[1])

        for j in range(len(file2_lines)):
            b = file2_lines[j]
            b1 = file2_lines[j].split('_')
            bb = float(b1[1])
            z = abs(aa-bb)          
            
            if z < 1 :
                print(f"merging {a} with {b}")
                merge = "psradd -R -o "+str(mpath)+str(pname)+"_"+str(aa)+"_3+5.fits " +str(a)+ ' '+str(b)
                os.system(merge)

    file1.close()
    file2.close()
    
    rm_cmd = 'rm -r '+temp_dir_b3+' '+temp_dir_b5+' '+metafile_b3+' '+metafile_b5    
    os.system(rm_cmd)

if __name__ == "__main__" :
    
    merge_fits() 

