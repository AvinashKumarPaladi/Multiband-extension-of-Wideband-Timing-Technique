#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np

# b3100 = np.genfromtxt("J1909-3744.band3.cc.100.meta",dtype="str")
# b4100 = np.genfromtxt("J1909-3744.band4.cc.100.meta",dtype="str")
# b5100 = np.genfromtxt("J1909-3744.band5.cc.100.meta",dtype="str")
# b3200 = np.genfromtxt("J1909-3744.band3.cc.200.meta",dtype="str")
# b5200 = np.genfromtxt("J1909-3744.band5.cc.200.meta",dtype="str")


# In[ ]:


import sys
b3100 = np.genfromtxt(sys.argv[1],dtype="str")
b5100 = np.genfromtxt(sys.argv[2],dtype="str")
#b3200 = np.genfromtxt(sys.argv[3],dtype="str")
#b5200 = np.genfromtxt(sys.argv[4],dtype="str")


# In[2]:


b3100_mjd = np.char.split(b3100,sep="_")
b3100_mjd = [el[1] for el in b3100_mjd]

# b4100_mjd = np.char.split(b4100,sep="_")
# b4100_mjd = [el[1] for el in b4100_mjd]

b5100_mjd = np.char.split(b5100,sep="_")
b5100_mjd = [el[1] for el in b5100_mjd]

#b3200_mjd = np.char.split(b3200,sep="_")
#b3200_mjd = [el[1] for el in b3200_mjd]

#b5200_mjd = np.char.split(b5200,sep="_")
#b5200_mjd = [el[1] for el in b5200_mjd]


# In[3]:


#b3200_35=[]
#b5200_35=[]
#b3200_35_mjd=[]
#b5200_35_mjd=[]
#i=0
#j=0
#while (i < len(b3200_mjd)) and (j < len(b5200_mjd)):
#    if -0.1<(float(b3200_mjd[i])-float(b5200_mjd[j]))<0.1:
#        b3200_35.append(b3200[i])
#        b5200_35.append(b5200[j])
#        b3200_35_mjd.append(float(b3200_mjd[i]))
#        b5200_35_mjd.append(float(b5200_mjd[j]))
#        i=i+1
#        j=j+1
#    elif (float(b3200_mjd[i])-float(b5200_mjd[j]))>0.1:
#        j=j+1
#    elif (float(b3200_mjd[i])-float(b5200_mjd[j]))<-0.1:
#        i=i+1
            


# In[4]:



#print((np.array(b3200_35_mjd)-np.array(b5200_35_mjd)))


# In[5]:


b3100_35=[]
b5100_35=[]
b3100_35_mjd=[]
b5100_35_mjd=[]
i=0
j=0
while (i < len(b3100_mjd)) and (j < len(b5100_mjd)):
    if -0.1<(float(b3100_mjd[i])-float(b5100_mjd[j]))<0.1:
        b3100_35.append(b3100[i])
        b5100_35.append(b5100[j])
        b3100_35_mjd.append(float(b3100_mjd[i]))
        b5100_35_mjd.append(float(b5100_mjd[j]))
        i=i+1
        j=j+1
    elif (float(b3100_mjd[i])-float(b5100_mjd[j]))>0.1:
        j=j+1
    elif (float(b3100_mjd[i])-float(b5100_mjd[j]))<-0.1:
        i=i+1
print(np.array(b3100_35_mjd)-np.array(b5100_35_mjd))


# In[11]:


np.savetxt(sys.argv[1]+".cc.b35",b3100_35,fmt="%s")
np.savetxt(sys.argv[2]+".cc.b35",b5100_35,fmt="%s")
#np.savetxt(sys.argv[3]+".cc.b35",b3200_35,fmt="%s")
#np.savetxt(sys.argv[4]+".cc.b35",b5200_35,fmt="%s")


# In[ ]:




