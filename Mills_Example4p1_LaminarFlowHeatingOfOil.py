
# coding: utf-8

# In[3]:


import numpy as np

mdot = 0.007 # mass flow rate
pipe_id = 0.01 #pipe inner diameter
pipe_perim = np.pi*pipe_id #pipe perimeter
pipe_length = 1.5 #length of pipe

k_b = 0.138 #Thermal conductivity at Tb
rho_b = 848.0 #density at Tb
cp_b = 2160 #specific heat at Tb
mu_b =  0.0255 #viscosity at Tb
nu_b = 30.1e-6 # kinematic viscosity at Tb
Pr_b = 400.0

ReD_b = 4*mdot/(pipe_perim*mu_b)
print(ReD_b)


# In[7]:


def NusseltContTwPipeLaminarCirc(ReD, Pr, DbyL):
    return 3.66 + 0.065*(DbyL)*ReD*Pr/(1+0.04*(DbyL*ReD*Pr)**(2.0/3))


# In[8]:


NuD_b = NusseltContTwPipeLaminarCirc(ReD_b, Pr_b, pipe_id/pipe_length)
print(NuD_b)


# In[9]:


mu_s = 0.503 # viscosity of SAE50 at Ts = 300 K
n= -0.11
NuD_corrected = NuD_b*(mu_s/mu_b)**n

print(NuD_corrected)


# In[10]:


h_avg = NuD_corrected*k_b/pipe_id
print(h_avg)


# In[11]:


Tw = 300.0
Tb_in = 377.0
Tb_out = Tw + (Tb_in - Tw)*np.exp(-h_avg*np.pi*pipe_id*pipe_length/(mdot*cp_b))
print(Tb_out)


# In[15]:


#Modified problem - rectangular channel of same area
ch_thickness = (np.pi/80)**0.5*pipe_id #
ch_width = 20*ch_thickness
ch_Ac = ch_thickness*ch_width 
ch_P = 2*(ch_thickness+ch_width)
print('pipe areas', ch_thickness*ch_width, np.pi*pipe_id**2/4)

Dh = 4*ch_Ac/ch_P

Re_Dh_b = 4*mdot/(ch_P*mu_b)

print('Re_Dh_b =', Re_Dh_b)



# In[18]:


def NusseltContTwPipeLaminarFlatChannel(ReDh, Pr, DbyL):
    return 7.54 + 0.03*(DbyL)*ReDh*Pr/(1+0.016*(DbyL*ReDh*Pr)**(2.0/3))


# In[19]:


NuDh_b = NusseltContTwPipeLaminarFlatChannel(Re_Dh_b, Pr_b, Dh/pipe_length)
print(NuDh_b)


# In[21]:


NuDh_corrected = NuDh_b*(mu_s/mu_b)**n
print(NuDh_corrected)


# In[24]:


h_avg_fc = NuDh_corrected*k_b/Dh
print(h_avg_fc)


# In[25]:


Tb_out_fc = Tw + (Tb_in - Tw)*np.exp(-h_avg_fc*ch_P*pipe_length/(mdot*cp_b))
print(Tb_out_fc)

