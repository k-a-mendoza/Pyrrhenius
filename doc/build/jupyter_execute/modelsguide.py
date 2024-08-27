#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
np.set_printoptions(edgeitems=3, infstr='inf', linewidth=75, nanstr='nan',
                   precision=8, suppress=False, threshold=3, formatter=None)

import pyrrhenius.database as phsd

ecdatabase = phsd.Database()


# In[2]:


nv17HD  = ecdatabase.get_model('nv_17_ol[010]') # Novella et al. 2017's hydrogen diffusion equation
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


# In[8]:


import pyrrhenius.database as phsd

ecdatabase = phsd.Database()

ecdatabase.get_model_list_for_phase('olivine')


# In[9]:


import matplotlib.pyplot as plt
import numpy as np
import pyrrhenius.database as phsd
import pyrrhenius.utils as pyhutils

ecdatabase = phsd.Database()

models = ['fei_20_ol[100]','fei_20_ol[010]', 'fei_20_ol[001]']


T = np.linspace(400,1800,num=120) # temperature in kelvin

fig, ax = plt.subplots()
linear_major_ticks = np.asarray([2000,1400,1100,900,800,700,600,500,400])
pyhutils.format_ax_arrhenian_space(ax,linear_major_ticks=linear_major_ticks,xlim=[5,10])
Cw = 100 # 100 ppm water
P = 3 # 3 GPa
for model_id in models:
    ecmodel = ecdatabase.get_model(model_id)
    c = ecmodel.get_conductivity(T=T,Cw=Cw,P=P)
    ax.plot(1e4/T,c,label=model_id)
ax.legend()
ax.set_ylim([1e-4,1e-1])
fig.show()


# In[10]:


ecdatabase = phsd.Database()
ecdatabase.create_isotropic_models()
ecdatabase.get_model_list_for_phase('olivine')[-12:]


# In[11]:


fig, ax = plt.subplots()
linear_major_ticks = np.asarray([2000,1400,1100,900,800,700,600,500,400])
pyhutils.format_ax_arrhenian_space(ax,linear_major_ticks=linear_major_ticks,xlim=[5,10])
Cw = 100 # 100 ppm water
P = 3 # 3 GPa
for model_id in models:
    ecmodel = ecdatabase.get_model(model_id)
    c = ecmodel.get_conductivity(T=T,Cw=Cw,P=P)
    ax.plot(1e4/T,c,label=model_id)

isotropic_model_id = 'isotropic_model:fei_20_ol[100]+fei_20_ol[010]+fei_20_ol[001]'
ecmodel = ecdatabase.get_model(isotropic_model_id)
c = ecmodel.get_conductivity(T=T,Cw=Cw,P=P)
ax.plot(1e4/T,c,label=model_id,linestyle='--',color='blue')
ax.legend()
ax.set_ylim([1e-4,1e-1])
fig.show()


# In[12]:


Cw = 100 # 100 ppm water
P = 3 # 3 GPa
T = np.linspace(1000,1500,num=10)
physiokwargs = {'Cw':Cw,'P':P,'T':T}
conductivity = ecmodel.get_conductivity(**physiokwargs)
conductivity


# In[13]:


conductivity = ecmodel.get_conductivity(averaging='max_aniso',**physiokwargs)
conductivity


# In[14]:


conductivity = ecmodel.get_conductivity(averaging='min_aniso',**physiokwargs)
conductivity


# In[15]:


conductivity = ecmodel.get_conductivity(crystal_direction='[100]',**physiokwargs)
conductivity


# In[16]:


conductivity = ecmodel.get_conductivity(crystal_direction=0.5,**physiokwargs)
conductivity


# In[17]:


import numpy as np

# Provide an array of anisotropic factors for batch querying
factors = np.random.uniform(0,1,size=len(T))
print(f'Anisotropic Factors: {factors}')
conductivity = ecmodel.get_conductivity(crystal_direction=factors,**physiokwargs)
conductivity


# In[18]:


fig, ax = plt.subplots()
linear_major_ticks = np.asarray([2000,1400,1100,900,800,700,600,500,400])
pyhutils.format_ax_arrhenian_space(ax,linear_major_ticks=linear_major_ticks,xlim=[5,10])

# Physiochemical states
Cw = 100 #(100 ppm water)
P = 3 # 3 GPa
T = np.linspace(1000,1500,num=10)
physiokwargs = {'Cw':Cw,'P':P,'T':T}


isotropic_model_id = 'isotropic_model:fei_20_ol[100]+fei_20_ol[010]+fei_20_ol[001]'
ecmodel = ecdatabase.get_model(isotropic_model_id)

for averaging in ['min_aniso','max_aniso','geometric']:
    c = ecmodel.get_conductivity(averaging=averaging,**physiokwargs)
    print(c)
    print(1e4/T)
    ax.plot(1e4/T,c,label=averaging,linewidth=2)

for factor in [0,0.3,0.5,0.6,1.0]:
    c = ecmodel.get_conductivity(crystal_direction=factor,**physiokwargs)
    ax.plot(1e4/T,c,label=f'f={factor}',linestyle='--',linewidth=1,marker='+')
ax.legend()
ax.set_ylim([1e-4,1e-1])
fig.show()


# In[19]:


from pyrrhenius.database import Database

ecdatabase = Database()

dry_model = ecdatabase.get_model('SEO3_ol')
novella_models =['nv_17_ol[100]','nv_17_ol[010]','nv_17_ol[001]']
fei_models = ['fei_20_ol_ionic[100]','fei_20_ol_ionic[010]', 'fei_20_ol_ionic[001]']


# In[20]:


from pyrrhenius.utils import calc_QFM
T = np.linspace(1000,2300,num=10) # in K
P = np.linspace(3,10,num=10) # in GPa
qfm = calc_QFM(T,P)
Cw = 300 # (in ppm)
physiochem = {'T':T,'P':P,'Cw':Cw,'logfo2':qfm}
models_to_composite = novella_models+fei_models

for id in models_to_composite:
    print('*'*20)
    print(id)
    ecmodel = ecdatabase.get_model(id)
    print(ecmodel.get_conductivity(**physiochem))
    new_ecmodel = ecmodel + dry_model
    print(new_ecmodel.get_conductivity(**physiochem))
    print(f'old ec model representation:{ecmodel}')
    print(f'new ec model representation:{new_ecmodel}')
    ecdatabase.register_new_model(new_ecmodel)


# In[21]:


ecdatabase.create_isotropic_models()
model_list = ecdatabase.get_model_list_for_phase('olivine')
model_list[-20:]


# In[22]:


from pyrrhenius.utils import calc_QFM
T = np.linspace(1000,3000,num=10) # in K
P = np.linspace(3,10,num=10) # in GPa
qfm = calc_QFM(T,P)
Cw = 300# (in ppm)
physiochem = {'T':T,'P':P,'Cw':Cw,'logfo2':qfm}


# In[23]:


fig, ax = plt.subplots()
linear_major_ticks = np.asarray([2000,1400,1100,900,800,700,600,500,400])
pyhutils.format_ax_arrhenian_space(ax,linear_major_ticks=linear_major_ticks,xlim=[5,10])


fei_ionic_isotropic  = 'isotropic_model:fei_20_ol_ionic[100]+fei_20_ol_ionic[010]+fei_20_ol_ionic[001]'
novella_isotropic    = 'isotropic_model:nv_17_ol[100]+nv_17_ol[010]+nv_17_ol[001]'
fei_ionic_compound   = 'isotropic_model:fei_20_ol_ionic[100]+SEO3_ol+fei_20_ol_ionic[010]+SEO3_ol+fei_20_ol_ionic[001]+SEO3_ol'
novella_wet_compound = 'isotropic_model:nv_17_ol[100]+SEO3_ol+nv_17_ol[010]+SEO3_ol+nv_17_ol[001]+SEO3_ol'
ecmodel = ecdatabase.get_model('SEO3_ol')
c = ecmodel.get_conductivity(**physiochem)
ax.plot(1e4/T,c,label='SEO3 Olivine',linewidth=2)

for label, model_id in zip(['isotropic fei','isotropic novella','iso fei w SEO3', 'iso Novella w SEO3'],
                    [fei_ionic_isotropic,novella_isotropic,fei_ionic_compound,novella_wet_compound ]):
    ecmodel = ecdatabase.get_model(model_id)
    c = ecmodel.get_conductivity(**physiochem)
    ax.plot(1e4/T,c,label=label,linewidth=2)


ax.legend()
ax.set_ylim([1e-4,1])
fig.show()

