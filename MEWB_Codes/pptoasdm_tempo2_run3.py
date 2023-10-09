#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import sys
import os
import numpy as np
from libstempo.libstempo import *
import libstempo as T

#### Input files
parfile = sys.argv[1]
timfile = sys.argv[2]
toa_output = sys.argv[3]
dm_output = sys.argv[4]
psr_name = sys.argv[5]

# parfile= 'J1909-3744.wb.inital.par' 
# timfile= 'J1909-3744.cc.tim'
# toa_output = "temp.res"
# dm_output = "temp.dm"
# psr_name = "temp"


#Make a new temperory par file and add DMX parameters
psr_temp = T.tempopulsar(
    parfile= parfile, 
    timfile= timfile,
    dofit=False
)
os.system('cp ./'+str(parfile)+' ./temp.par')
data=psr_temp.stoas
ntoa=len(data)

for i in range(ntoa):
    number_str = str(i+1)
    zero_filled_number = number_str. zfill(4) 
    os.system("echo 'DMX_"+str(zero_filled_number)+"\t0.0\t1\t0.01' >> temp.par")
    os.system("echo 'DMXEP_"+str(zero_filled_number)+"\t"+str(data[i].round(5))+"' >> temp.par")
    os.system("echo 'DMXR1_"+str(zero_filled_number)+"\t"+str(data[i].round(5)-0.025)+"' >> temp.par")
    os.system("echo 'DMXR2_"+str(zero_filled_number)+"\t"+str(data[i].round(5)+0.025)+"' >> temp.par")

parfile="temp.par"
psr = T.tempopulsar(
    parfile= parfile, 
    timfile= timfile,
    dofit=False
)


# Extract T2EFAC, T2EUAD, DMEFAC from parfile
wnp_parfile = {}
wnp_parfile["wnp_id"]=[]
wnp_parfile["flag_id"]=[]
wnp_parfile["flag"]=[]
wnp_parfile["wnp_val"]=[]

text_file = open(parfile, "r")
lines = text_file.readlines()
for iline,line in enumerate(lines):
        if ('T2EFAC' in line) or ('T2EQUAD' in line) or ('DMEFAC' in line):
            toks = line.strip().split()
            wnp_parfile["wnp_id"].append(toks[0])
            wnp_parfile["flag_id"].append(toks[1].replace("-",""))
            wnp_parfile["flag"].append(toks[2])
            wnp_parfile["wnp_val"].append(float(toks[3]))

wnp_parfile["wnp_id"]=np.array(wnp_parfile["wnp_id"])
wnp_parfile["flag_id"]=np.array(wnp_parfile["flag_id"])
wnp_parfile["flag"]=np.array(wnp_parfile["flag"])
wnp_parfile["wnp_val"]=np.array(wnp_parfile["wnp_val"])

wnp={}
wnp["T2EFAC"]= np.ones((psr.nobs))
wnp["T2EQUAD"]= np.zeros((psr.nobs))
wnp["DMEFAC"]= np.ones((psr.nobs))
for i in range(len(wnp_parfile["wnp_id"])):
    wnp[wnp_parfile["wnp_id"][i]][np.where(psr.flagvals(wnp_parfile["flag_id"][i]) ==                                            wnp_parfile["flag"][i])[0]]= wnp_parfile["wnp_val"][i]

    
#Defining GLS function
def glsfit(psr,ntoa,wnp, renormalize=True):
    """Solve local GLS problem using scipy.linalg.cholesky.
    Update psr[...].val and psr[...].err from solution.
    If renormalize=True, normalize each design-matrix column by its norm."""

    mask = psr.deleted == 0
    
#     res, err = psr.residuals(removemean=False)[mask], psr.toaerrs[mask]
#     M = psr.designmatrix(updatebats=False, fixunits=True, fixsigns=True)[mask, :]

    efac = wnp["T2EFAC"]
    print("t2efac:",efac)
    equad = wnp["T2EQUAD"]
    dmefac = wnp["DMEFAC"]
    print("dmefac",dmefac)
    
    delta_t = psr.residuals(removemean=False)
    delta_t = np.array([list(delta_t)])
    sigma_t = 1.0e-6*psr.toaerrs
    sigma_t = ((efac*sigma_t)**2+equad**2)
    dm = psr.flagvals('pp_dm')
    dm = dm.astype(np.float64)
    sigma_dm = psr.flagvals('pp_dme')
    sigma_dm = sigma_dm.astype(np.float64)
    sigma_dm = (dmefac*sigma_dm)**2
    fid_dm = psr['DM'].val
    delta_dm = dm - fid_dm
    delta_dm = np.array([list(delta_dm)])
#     print("ok")
#     print(delta_t.shape,sigma_t.shape,delta_dm.shape,sigma_dm.shape)
    
    sigma = np.concatenate((sigma_t,sigma_dm))
    res = np.concatenate((delta_t,delta_dm),axis=1)
    res = np.array(res[0])
    C = numpy.diag((sigma))
    
    M_temp = psr.designmatrix()
    I = np.identity(ntoa)
    zero_mat = np.zeros((ntoa,np.shape(M_temp)[1]-ntoa))
    M2_temp = np.concatenate((zero_mat,I),axis=1)
    M = np.concatenate((M_temp,M2_temp),axis=0)
#     print("M_temp,I,zero_mat,M2_temp,M")
#     print(M_temp.shape,I.shape,zero_mat.shape,M2_temp.shape,M.shape)
    

    if renormalize:
        norm = numpy.sqrt(numpy.sum(M ** 2, axis=0))
        M /= norm
    else:
        norm = numpy.ones_like(M[0, :])

    mtcm = numpy.dot(M.T, numpy.dot(numpy.linalg.inv(C), M))
    mtcy = numpy.dot(M.T, numpy.dot(numpy.linalg.inv(C), res))
    xvar = numpy.linalg.inv(mtcm)
    c = scipy.linalg.cho_factor(mtcm)
    xhat = scipy.linalg.cho_solve(c, mtcy)

    
    sol = psr.vals()
    psr.vals(sol + xhat[1:] / norm[1:])
    psr.errs(numpy.sqrt(numpy.diag(xvar)[1:]) / norm[1:])
    
    p1 = delta_t - numpy.dot(M_temp,(xhat/norm))
    n1 = len(xhat) - ntoa
    p2 = xhat[n1:] / norm[n1:] - delta_dm
    N_t = numpy.linalg.inv(numpy.diag((sigma_t)))
    N_dm = numpy.linalg.inv(numpy.diag((sigma_dm)))
    chi2 = np.dot(p1,np.dot(N_t,p1.T)) + np.dot(p2,np.dot(N_dm,p2.T))
    redchi2 = chi2/(2*ntoa-len(xhat))
    
    return redchi2
#     return 1

#Running GLS and saving new parfile
redchi2 = glsfit(psr,ntoa, wnp, renormalize=True)
psr.savepar(psr_name)


#tempo_chi2 = (psr.residuals()**2)/()
#print(\n)
# print("Total postfit chi2: ", tempo_chi2)
# print("Degree of freedom: ", ndof)
# print("reduced chi2: ", tempo_chi2/ndof)
# print("Postfit weighted rms: ", rms_t)

efac = wnp["T2EFAC"]
equad = wnp["T2EQUAD"]
dmefac = wnp["DMEFAC"]
    
#Writing residuals to a file
epoch = np.reshape(np.array((psr.toas().astype(np.float128))),(-1,1))
toa = np.reshape(np.array(( 1.0e6*psr.residuals().astype(np.float128))),(-1,1))
sigma_toa = np.sqrt((efac*psr.toaerrs.astype(np.float128))**2+equad**2)
toae = np.reshape(np.array((sigma_toa)),(-1,1))
mat_int1 = np.concatenate((epoch,toa,toae),axis=1)
np.savetxt(toa_output,mat_int1)#,fmt='%.6e')

#Writing DM to a file
epoch = np.reshape(np.array((psr.toas().astype(np.float128))),(-1,1))
dM = np.reshape(np.array(( psr.flagvals('pp_dm').astype(np.float128))),(-1,1))
sigma_dm = np.sqrt((efac*psr.flagvals('pp_dme').astype(np.float128))**2+equad**2)
dMe = np.reshape(np.array((sigma_dm)),(-1,1))
mat_int2 = np.concatenate((epoch,dM,dMe),axis=1)
np.savetxt(dm_output,mat_int2)#,fmt='%.6e')

#Print median DM
median_dm = np.median(dM)
print ('median DM : %f'  %median_dm)

#Print WRMS
# print(toa,toae)
# print(type(toa))
# print(type(toae))
#x_mean = np.sum((1/toae**2)*toa)/np.sum(1/toae**2)
#WRMS = np.sum((1/toae**2)*((toa-x_mean)**2))/np.sum(1/toae**2)
#print('median DM :'+str(np.sqrt(WRMS)))

#Print chi2
# print(toa,toae)
# print(type(toa))
# print(type(toae))
x_mean = np.sum((1/toae**2)*toa)/np.sum(1/toae**2)
WRMS = np.sum((1/toae**2)*((toa-x_mean)**2))/np.sum(1/toae**2)
print('correct WRMS :'+str(np.sqrt(WRMS)))


print('correct redchi2 :'+str(redchi2[0][0]))
# # tempo_chi2 = 
# #print(\n)
# # print("Total postfit chi2: ", tempo_chi2)
# # print("Degree of freedom: ", ndof)
# # print("reduced chi2: ", tempo_chi2/ndof)
# # print("Postfit weighted rms: ", rms_t)

# #Writing residuals to a file
# f1=open(toa_output,'w')
# epoch = list(map(float,psr.toas()))
# toa = list(map(float, 1.0e6*psr.residuals()))
# toae = list(map(float, psr.toaerrs))

# for i in range(len(epoch)):
#     xgrid_int=np.column_stack((epoch[i], toa[i],toae[i]))
#     mat_int = np.matrix(xgrid_int)
#     np.savetxt(f1,mat_int)#,fmt='%.6e')

# #Writing DM to a file
# f2=open(dm_output,'w')
# epoch = list(map(float,psr.toas()))
# dM = list(map(float, psr.flagvals('pp_dm')))
# dMe = list(map(float, psr.flagvals('pp_dme')))

# for i in range(len(epoch)):
#     xgrid_int=np.column_stack((epoch[i], dM[i], dMe[i]))
#     mat_int = np.matrix(xgrid_int)
#     np.savetxt(f2,mat_int)#,fmt='%.6e')

# #Print median DM
# median_dm = np.median(dM)
# print ('median DM : %f'  %median_dm)


# In[ ]:




