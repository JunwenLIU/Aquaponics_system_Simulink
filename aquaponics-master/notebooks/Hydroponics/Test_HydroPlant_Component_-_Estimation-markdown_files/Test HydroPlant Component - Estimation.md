

```python
import numpy as np
import matplotlib.pyplot as plt

from aquaponics import Aquaponics

def plot():
    %matplotlib inline
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
```

## Import Data


```python
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
```


```python
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
    
```

    MHE results: K: 0.0383855341, Jmax_0: 1.39389731, alpha: 0.145662984
    MHE results: K: 0.0463971094, Jmax_0: 1.31529235, alpha: 0.148303906
    MHE results: K: 0.0476669247, Jmax_0: 1.07240462, alpha: 0.150759132
    MHE results: K: 0.0510010176, Jmax_0: 1.13206634, alpha: 0.149764651
    MHE results: K: 0.0507377166, Jmax_0: 1.46187781, alpha: 0.149424234
    MHE results: K: 0.0519622151, Jmax_0: 1.38683966, alpha: 0.149828549
    MHE results: K: 0.0512883602, Jmax_0: 1.32093841, alpha: 0.149686349
    MHE results: K: 0.0514967108, Jmax_0: 1.28209859, alpha: 0.149725941
    MHE results: K: 0.0521183105, Jmax_0: 1.25770786, alpha: 0.149879027
    MHE results: K: 0.0520119576, Jmax_0: 1.24291089, alpha: 0.149843006
    MHE results: K: 0.0522402544, Jmax_0: 1.2329404, alpha: 0.14990418
    MHE results: K: 0.0523860618, Jmax_0: 1.22584585, alpha: 0.149947535
    MHE results: K: 0.0523841711, Jmax_0: 1.2213514, alpha: 0.149945378
    MHE results: K: 0.0517999832, Jmax_0: 1.21821512, alpha: 0.149733459
    MHE results: K: 0.0518066859, Jmax_0: 1.21608422, alpha: 0.149722177
    MHE results: K: 0.0524122528, Jmax_0: 1.21445943, alpha: 0.149950243
    MHE results: K: 0.0524189671, Jmax_0: 1.21306159, alpha: 0.149952037
    MHE results: K: 0.051982156, Jmax_0: 1.21223121, alpha: 0.149756244
    MHE results: K: 0.0522106738, Jmax_0: 1.21134849, alpha: 0.149851464
    MHE results: K: 0.0525125206, Jmax_0: 1.21076314, alpha: 0.149992201
    MHE results: K: 0.0524900973, Jmax_0: 1.21007629, alpha: 0.149981709
    MHE results: K: 0.0521508388, Jmax_0: 1.2096966, alpha: 0.149800483
    MHE results: K: 0.0520179457, Jmax_0: 1.20947685, alpha: 0.149720007
    MHE results: K: 0.0527804585, Jmax_0: 1.20778671, alpha: 0.150159659
    MHE results: K: 0.0524201653, Jmax_0: 1.2076685, alpha: 0.149940734
    MHE results: K: 0.0529273918, Jmax_0: 1.20630389, alpha: 0.150261078
    MHE results: K: 0.05200844, Jmax_0: 1.20719476, alpha: 0.149668965
    MHE results: K: 0.0523653255, Jmax_0: 1.20686112, alpha: 0.149912201
    MHE results: K: 0.0527060472, Jmax_0: 1.20609311, alpha: 0.150143678
    MHE results: K: 0.0528766044, Jmax_0: 1.2053302, alpha: 0.150267016
    MHE results: K: 0.0520690286, Jmax_0: 1.20611857, alpha: 0.149670247
    MHE results: K: 0.0523299789, Jmax_0: 1.20606329, alpha: 0.149859297
    MHE results: K: 0.0527018667, Jmax_0: 1.20530299, alpha: 0.150150791
    MHE results: K: 0.0528200163, Jmax_0: 1.20499226, alpha: 0.150260542
    MHE results: K: 0.0529512194, Jmax_0: 1.20444905, alpha: 0.150378353
    MHE results: K: 0.0529922465, Jmax_0: 1.20420981, alpha: 0.150424625
    MHE results: K: 0.0530071411, Jmax_0: 1.2038635, alpha: 0.150445923
    MHE results: K: 0.0528219887, Jmax_0: 1.20363657, alpha: 0.150287208
    MHE results: K: 0.0528806791, Jmax_0: 1.20327887, alpha: 0.150347517
    MHE results: K: 0.0528192878, Jmax_0: 1.20294943, alpha: 0.150309748
    MHE results: K: 0.0523197148, Jmax_0: 1.20357591, alpha: 0.149818785
    MHE results: K: 0.0525642877, Jmax_0: 1.20339091, alpha: 0.150058834
    MHE results: K: 0.0526995842, Jmax_0: 1.20310176, alpha: 0.150196115
    MHE results: K: 0.0527613879, Jmax_0: 1.20322027, alpha: 0.150267551
    MHE results: K: 0.052230516, Jmax_0: 1.20380264, alpha: 0.149708093
    MHE results: K: 0.0525087297, Jmax_0: 1.20349413, alpha: 0.150005697
    MHE results: K: 0.05268245, Jmax_0: 1.20316063, alpha: 0.150197914
    MHE results: K: 0.0523505361, Jmax_0: 1.20449763, alpha: 0.149844722
    MHE results: K: 0.0525547384, Jmax_0: 1.20407635, alpha: 0.150070654
    MHE results: K: 0.0523443533, Jmax_0: 1.20387104, alpha: 0.149811862
    MHE results: K: 0.0524646101, Jmax_0: 1.20395338, alpha: 0.149956818
    MHE results: K: 0.0525967757, Jmax_0: 1.20374387, alpha: 0.150121688
    MHE results: K: 0.0527285939, Jmax_0: 1.20321611, alpha: 0.150281852
    MHE results: K: 0.052731256, Jmax_0: 1.20317127, alpha: 0.150288868
    MHE results: K: 0.0523231798, Jmax_0: 1.20320758, alpha: 0.149768412
    MHE results: K: 0.0524160234, Jmax_0: 1.20362894, alpha: 0.149890351
    MHE results: K: 0.0525170346, Jmax_0: 1.20349081, alpha: 0.15002277
    MHE results: K: 0.0526203628, Jmax_0: 1.20319903, alpha: 0.150161338
    MHE results: K: 0.0525919256, Jmax_0: 1.20333268, alpha: 0.15012625
    MHE results: K: 0.0523283958, Jmax_0: 1.20308507, alpha: 0.149757038
    MHE results: K: 0.0524845474, Jmax_0: 1.20264955, alpha: 0.149980487
    MHE results: K: 0.0525461205, Jmax_0: 1.20270746, alpha: 0.150067173
    MHE results: K: 0.0526439989, Jmax_0: 1.20244949, alpha: 0.150209406
    MHE results: K: 0.0526733135, Jmax_0: 1.20209267, alpha: 0.150259242
    MHE results: K: 0.0525422276, Jmax_0: 1.20204045, alpha: 0.150059993
    MHE results: K: 0.0523777664, Jmax_0: 1.20228256, alpha: 0.149809542
    MHE results: K: 0.0525260441, Jmax_0: 1.20212869, alpha: 0.150034439
    MHE results: K: 0.0526304297, Jmax_0: 1.20196832, alpha: 0.150201151
    MHE results: K: 0.0524239713, Jmax_0: 1.20229053, alpha: 0.149889377
    MHE results: K: 0.0523474984, Jmax_0: 1.20329708, alpha: 0.149760296
    MHE results: K: 0.0523833819, Jmax_0: 1.20283442, alpha: 0.149805109
    MHE results: K: 0.0524479635, Jmax_0: 1.20281555, alpha: 0.149910721
    MHE results: K: 0.0525044719, Jmax_0: 1.20269928, alpha: 0.150004482
    MHE results: K: 0.05236439, Jmax_0: 1.20295085, alpha: 0.149775386
    MHE results: K: 0.052443541, Jmax_0: 1.20280785, alpha: 0.149905451
    MHE results: K: 0.0523777535, Jmax_0: 1.20250359, alpha: 0.14978285
    MHE results: K: 0.0524318943, Jmax_0: 1.20236468, alpha: 0.149887286
    MHE results: K: 0.0525075199, Jmax_0: 1.20224448, alpha: 0.150017169
    MHE results: K: 0.0524988368, Jmax_0: 1.20234534, alpha: 0.149999132
    MHE results: K: 0.0525153586, Jmax_0: 1.20209697, alpha: 0.150031686
    MHE results: K: 0.052411601, Jmax_0: 1.20198186, alpha: 0.149832997
    MHE results: K: 0.0524731017, Jmax_0: 1.20203375, alpha: 0.149947377
    MHE results: K: 0.052517791, Jmax_0: 1.20195205, alpha: 0.150030541
    MHE results: K: 0.0525227953, Jmax_0: 1.20200108, alpha: 0.150040716
    MHE results: K: 0.0525000721, Jmax_0: 1.20211558, alpha: 0.149997644
    MHE results: K: 0.052406759, Jmax_0: 1.20196274, alpha: 0.149813878
    MHE results: K: 0.0524630112, Jmax_0: 1.20191139, alpha: 0.149922851
    MHE results: K: 0.0524918105, Jmax_0: 1.2019159, alpha: 0.149979203
    MHE results: K: 0.0525006921, Jmax_0: 1.20192758, alpha: 0.149996941
    MHE results: K: 0.0524211525, Jmax_0: 1.20175085, alpha: 0.149836219
    MHE results: K: 0.0524259669, Jmax_0: 1.20193875, alpha: 0.149848635
    MHE results: K: 0.0524354752, Jmax_0: 1.20189828, alpha: 0.149866785
    MHE results: K: 0.0524693821, Jmax_0: 1.20186246, alpha: 0.14993442
    MHE results: K: 0.0524265648, Jmax_0: 1.20167603, alpha: 0.149840959
    MHE results: K: 0.0524946676, Jmax_0: 1.20148114, alpha: 0.149987379
    MHE results: K: 0.0525244144, Jmax_0: 1.20132288, alpha: 0.150051687
    MHE results: K: 0.0525103699, Jmax_0: 1.20147088, alpha: 0.150020287
    MHE results: K: 0.0525388951, Jmax_0: 1.2012972, alpha: 0.150084031
    MHE results: K: 0.0524634534, Jmax_0: 1.20127119, alpha: 0.14991569
    MHE results: K: 0.0524988564, Jmax_0: 1.2011453, alpha: 0.149995183
    MHE results: K: 0.0524927998, Jmax_0: 1.2013637, alpha: 0.149980081
    MHE results: K: 0.0524600079, Jmax_0: 1.20128771, alpha: 0.14990524
    MHE results: K: 0.0524497991, Jmax_0: 1.20125458, alpha: 0.149880577
    MHE results: K: 0.0524866246, Jmax_0: 1.20113069, alpha: 0.149966001
    MHE results: K: 0.0525028684, Jmax_0: 1.20104477, alpha: 0.150004355
    MHE results: K: 0.0525111492, Jmax_0: 1.20118404, alpha: 0.15002273
    MHE results: K: 0.0525281047, Jmax_0: 1.20106664, alpha: 0.150064187
    MHE results: K: 0.052528685, Jmax_0: 1.20098032, alpha: 0.150066569
    MHE results: K: 0.0524532867, Jmax_0: 1.20103205, alpha: 0.149883443
    MHE results: K: 0.0524361534, Jmax_0: 1.20140549, alpha: 0.149840064
    MHE results: K: 0.0524571891, Jmax_0: 1.20137381, alpha: 0.149890564
    MHE results: K: 0.0524683347, Jmax_0: 1.20132157, alpha: 0.149917327
    MHE results: K: 0.0524768368, Jmax_0: 1.2013409, alpha: 0.149938819
    MHE results: K: 0.0524547459, Jmax_0: 1.20121149, alpha: 0.149882357
    MHE results: K: 0.0524764881, Jmax_0: 1.20120307, alpha: 0.149936293
    MHE results: K: 0.0524859744, Jmax_0: 1.2012015, alpha: 0.149959697
    MHE results: K: 0.0524505947, Jmax_0: 1.20116013, alpha: 0.149868124
    MHE results: K: 0.0524550966, Jmax_0: 1.20112307, alpha: 0.149878773
    MHE results: K: 0.0524766537, Jmax_0: 1.20110315, alpha: 0.149934222



```python
plot()
```


![png](Test%20HydroPlant%20Component%20-%20Estimation_files/Test%20HydroPlant%20Component%20-%20Estimation_4_0.png)



```python
plot()
```


```python
a = Aquaponics('hydroplant', T0=20, N0=0.1, K=K0, 
                Jmax_0=Jmax_00, alpha=alpha0, kswitch=1)

m = a.get_model()

m.time = x
a.solve(glamdring=True, imode=7, disp=False)
```


```python
plot()
```
