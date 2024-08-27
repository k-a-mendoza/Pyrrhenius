#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
np.set_printoptions(edgeitems=3, infstr='inf', linewidth=75, nanstr='nan',
                   precision=8, suppress=False, threshold=3, formatter=None)

import pyrrhenius.database as phsd

database_location = 'database/publication_database.csv'
ecdatabase = phsd.Database(database_location)


# In[2]:


nv17HD  = ecdatabase.get_model('nv_17_ol[010]') # Novella et al. 2017's HD diffusion equation
print(nv17HD)

seo3dry =  ecdatabase.get_model('SEO3_ol') # Constable et al. 2006's dry olivine equation
print(seo3dry)

combined = nv17HD+seo3dry
print(combined)


# In[3]:


model = ecdatabase.get_model('SEO3_ol')
model.get_conductivity(T=1000, P=1.0, logfo2=10**-11)


# In[4]:


T = np.ones(4)*700 # in degrees K
print(T.shape)
model.get_conductivity(T=T, P=1.0, logfo2=10**-11)


# In[5]:


T = np.ones((4,4))*700 # in degrees K
print(T.shape)
model.get_conductivity(T=T, P=1.0, logfo2=10**-11)


# In[6]:


try:
    model.get_conductivity(T=T)
except AssertionError as e:
    print(e)


# In[7]:


model.print_required_parameters()

