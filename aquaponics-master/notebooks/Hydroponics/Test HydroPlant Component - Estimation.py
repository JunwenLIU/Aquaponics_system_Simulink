
# coding: utf-8

# In[3]:


import numpy as np
import matplotlib.pyplot as plt

from aquaponics import Aquaponics

def plot():
    get_ipython().run_line_magic('matplotlib', 'inline')
    plt.figure(figsize=(12,9))
    ax = plt.subplot(311)
    plt.plot(x, w_mhe, label='Dry Weight')
    plt.plot(x, y, 'r:', label='Data')
    plt.grid()
    plt.legend()
    plt.ylabel('Plant Weight (g)')

    plt.subplot(312, sharex=ax)
    plt.plot(x, dNup_mhe, label='dNup')
    plt.grid()
    plt.legend()
    plt.ylabel('Nitrogen Changes (mmol / day)')

    plt.subplot(313, sharex=ax)
    plt.plot(x, solar/1e6)
    plt.grid()
    plt.ylabel('Light Intensity (MJ)')

    #plt.xlim(0, tf)
    plt.xlabel('Time (days)')


# ## Import Data

# In[4]:


data = np.loadtxt('exp1.csv', delimiter=',')
x = data[:,0]
y = data[:,1]
# plt.plot(x, y)

solar_constant = 7.1e6  # Megajoules
solar_on = 0.25       # On at 6 AM
solar_off = 0.75      # Off at 6 PM
solar = np.zeros(len(x))
for i in range(len(x)):
    if solar_on < x[i] % 1 < solar_off:
        solar[i] = solar_constant

alpha0=0.150
Jmax_00 = 0.0498*24
K0 = 0.0525

a = Aquaponics('hydroplant', T0=20, N0=0.1, K=K0, 
                Jmax_0=Jmax_00, alpha=alpha0, kswitch=1)

m = a.get_model()

a.K.FSTATUS=0
a.K.STATUS=1
a.K.UPPER = 0.1
a.K.LOWER = 0.005

a.Jmax_0.FSTATUS=0
a.Jmax_0.STATUS=1
a.Jmax_0.LOWER = 0.02*24
a.Jmax_0.UPPER = 0.08*24

a.alpha.FSTATUS=0
a.alpha.STATUS=1
a.alpha.LOWER = 0.10
a.alpha.UPPER = 0.2

a.w.FSTATUS = 1
a.w.STATUS = 1
a.w.MEAS_GAP = 0.1

#6 time points in horizon
m.time = np.linspace(0,.1,2)


# In[5]:


# Allocate storage
w_mhe = np.zeros(len(x))
K_mhe = np.ones(len(x)) * K0
Jmax_0_mhe = np.ones(len(x)) * Jmax_00
alpha_mhe = np.ones(len(x)) * alpha0
dNup_mhe = np.zeros(len(x)) 

for i in range(0, len(x)):
    
    a.I.value = solar[i]
    a.w.MEAS = y[i]    
    a.solve(glamdring=True, imode=5, disp=False)
    
    # check if successful
    if m.options.APPSTATUS == 1:
        # retrieve solution
        w_mhe[i] = a.w.MODEL
        K_mhe[i] = a.K.NEWVAL
        Jmax_0_mhe[i] = a.Jmax_0.NEWVAL
        alpha_mhe[i] = a.alpha.NEWVAL
        dNup_mhe[i] = a.dNup.VALUE[-1]
    else:
        # failed solution
        w_mhe[i] = 0
        K_mhe[i] = 0
        Jmax_0_mhe[i] = 0
        alpha[i] = 0
        dNup_mhe[i] = 0

    print('MHE results: K: {}, Jmax_0: {}, alpha: {}'.format(
        K_mhe[i], Jmax_0_mhe[i], alpha_mhe[i]))
    


# In[6]:


plot()


# In[ ]:


plot()


# In[ ]:


a = Aquaponics('hydroplant', T0=20, N0=0.1, K=K0, 
                Jmax_0=Jmax_00, alpha=alpha0, kswitch=1)

m = a.get_model()

m.time = x
a.solve(glamdring=True, imode=7, disp=False)


# In[ ]:


plot()

