# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 18:13:53 2021

@author: KEVIN-DORGANS
"""

import scipy as sp

NeuN_Array = []
MAP2_Array = []
AlphaTub_Array = []
NeuN_POSITIVE = []
CELL_AREA = []

PATHS = load_directory_content_and_sub__()[0]

for i in range(len(PATHS)):
    if ('csv' in PATHS[i])==True and ('ZSTACK_NeuN' in PATHS[i])==True:
        NeuN_Array_temp = pd.read_csv(PATHS[i])
        NeuN_Array.append(np.transpose(NeuN_Array_temp)[1::])
    elif ('csv' in PATHS[i])==True and ('ZSTACK_MAP2' in PATHS[i])==True:
        MAP2_Array_temp = pd.read_csv(PATHS[i])
        MAP2_Array.append(np.transpose(MAP2_Array_temp)[1::])
    elif ('csv' in PATHS[i])==True and ('ZSTACK_AlphaTub' in PATHS[i])==True:
        AlphaTub_Array_temp = pd.read_csv(PATHS[i])
        AlphaTub_Array.append(np.transpose(AlphaTub_Array_temp)[1::])
    elif ('csv' in PATHS[i])==True and ('ZSTACK_AREA' in PATHS[i])==True:
        CELL_AREA_temp = pd.read_csv(PATHS[i])
        CELL_AREA.append(np.array(CELL_AREA_temp)[:,1])

#RADIUS = NeuN_Array[0][0]
AlphaTub_Array = np.concatenate(AlphaTub_Array)
NeuN_Array = np.concatenate(NeuN_Array)
MAP2_Array = np.concatenate(MAP2_Array)
CELL_AREA = np.concatenate(CELL_AREA)

plt.figure(figsize=(6,4))
ax1 = plt.subplot(231)
ax2 = plt.subplot(232)
ax3 = plt.subplot(233)

ax4 = plt.subplot(234, sharey=ax1, sharex=ax1)
ax5 = plt.subplot(235, sharey=ax2, sharex=ax2)
ax6 = plt.subplot(236, sharey=ax3, sharex=ax3)

NeuN_Array_BackgroundSUB = [NeuN_Array[i]-np.nanmin(NeuN_Array[i]) for i in range(len(NeuN_Array))]
NeuN_Array_BackgroundSUB = [NeuN_Array_BackgroundSUB [i]/np.nanmax(NeuN_Array_BackgroundSUB ) for i in range(len(NeuN_Array_BackgroundSUB ))]


MEAN = np.nanmean(NeuN_Array_BackgroundSUB, axis=0)
RADIUS = np.linspace(0, 14, len(MEAN))
#ax1.plot(RADIUS, MEAN, color='purple')
#ax4.plot(RADIUS, MEAN, color='purple')

NeuN_ASTRO = []
NeuN_NEURO = []
MAP2_ASTRO = []
MAP2_NEURO = []
AlphaTub_ASTRO = []
AlphaTub_NEURO = []
AREA_NEURO = []
AREA_ASTRO = []
NeuNIntensity = []
MAP2Intensity = []

for i in range(len(NeuN_Array_BackgroundSUB)):
    NeuNIntensity.append(np.nanmean(NeuN_Array_BackgroundSUB[i][0:12]))
    if np.nanmean(NeuN_Array_BackgroundSUB[i][0:12])<0.2:
        ax1.plot(RADIUS, NeuN_Array_BackgroundSUB[i], lw=0.1, color='purple')
        NeuN_POSITIVE.append(False)
        NeuN_ASTRO.append(NeuN_Array_BackgroundSUB[i])
    else:
        ax4.plot(RADIUS, NeuN_Array_BackgroundSUB[i], lw=0.1, color='purple')
        NeuN_POSITIVE.append(True)
        NeuN_NEURO.append(NeuN_Array_BackgroundSUB[i])

MAP2_Array_BackgroundSUB = [MAP2_Array[i]-np.nanmin(MAP2_Array[i]) for i in range(len(MAP2_Array))]
MAP2_Array_BackgroundSUB = [MAP2_Array_BackgroundSUB [i]/np.nanmax(MAP2_Array_BackgroundSUB ) for i in range(len(MAP2_Array_BackgroundSUB ))]


MEAN = np.nanmean(MAP2_Array_BackgroundSUB, axis=0)
RADIUS = np.linspace(0, 14, len(MEAN))
#ax2.plot(RADIUS, MEAN, color='blue')
for i in range(len(NeuN_Array_BackgroundSUB)):
    MAP2Intensity.append(np.nanmax(MAP2_Array_BackgroundSUB[i][6:27]))
    if NeuN_POSITIVE[i]==False:
        ax2.plot(RADIUS, MAP2_Array_BackgroundSUB[i], lw=0.1, color='blue')
        MAP2_ASTRO.append(MAP2_Array_BackgroundSUB[i])
    else:
        ax5.plot(RADIUS, MAP2_Array_BackgroundSUB[i], lw=0.1, color='blue')
        MAP2_NEURO.append(MAP2_Array_BackgroundSUB[i])

AlphaTub_Array_BackgroundSUB = [AlphaTub_Array[i]-np.nanmin(AlphaTub_Array[i]) for i in range(len(AlphaTub_Array))]
AlphaTub_Array_BackgroundSUB = [AlphaTub_Array_BackgroundSUB[i]/np.nanmax(AlphaTub_Array_BackgroundSUB) for i in range(len(AlphaTub_Array_BackgroundSUB))]

MEAN = np.nanmean(AlphaTub_Array_BackgroundSUB, axis=0)
RADIUS = np.linspace(0, 14, len(MEAN))
#ax3.plot(RADIUS, MEAN)
AlphaTub_PEAK_DISTANCE = []
for i in range(len(NeuN_Array_BackgroundSUB)):
    AlphaTub_PEAK_DISTANCE.append( RADIUS[AlphaTub_Array_BackgroundSUB[i].tolist().index(np.nanmax(AlphaTub_Array_BackgroundSUB[i]))])
    if NeuN_POSITIVE[i]==False:
        ax3.plot(RADIUS, AlphaTub_Array_BackgroundSUB[i], lw=0.1, color='green')
        AlphaTub_ASTRO.append(AlphaTub_Array_BackgroundSUB[i])
        AREA_ASTRO.append(CELL_AREA[i])
    else:
        ax6.plot(RADIUS, AlphaTub_Array_BackgroundSUB[i], lw=0.1, color='green')
        AlphaTub_NEURO.append(AlphaTub_Array_BackgroundSUB[i])
        AREA_NEURO.append(CELL_AREA[i])

MEAN = np.nanmean(NeuN_ASTRO, axis=0)
SEM = sp.stats.sem(NeuN_ASTRO, axis=0)
ax1.plot(RADIUS, MEAN)
ax1.fill_between(RADIUS, MEAN+SEM, MEAN-SEM, alpha=0.1)

MEAN = np.nanmean(NeuN_NEURO, axis=0)
SEM = sp.stats.sem(NeuN_NEURO, axis=0)
ax4.plot(RADIUS, MEAN)
ax4.fill_between(RADIUS, MEAN+SEM, MEAN-SEM, alpha=0.1)

MEAN = np.nanmean(MAP2_ASTRO, axis=0)
SEM = sp.stats.sem(MAP2_ASTRO, axis=0)
ax2.plot(RADIUS, MEAN)
ax2.fill_between(RADIUS, MEAN+SEM, MEAN-SEM, alpha=0.1)

MEAN = np.nanmean(MAP2_NEURO, axis=0)
SEM = sp.stats.sem(MAP2_NEURO, axis=0)
ax5.plot(RADIUS, MEAN)
ax5.fill_between(RADIUS, MEAN+SEM, MEAN-SEM, alpha=0.1)

MEAN = np.nanmean(AlphaTub_ASTRO, axis=0)
SEM = sp.stats.sem(AlphaTub_ASTRO, axis=0)
ax3.plot(RADIUS, MEAN)
ax3.fill_between(RADIUS, MEAN+SEM, MEAN-SEM, alpha=0.1)

MEAN = np.nanmean(AlphaTub_NEURO, axis=0)
SEM = sp.stats.sem(AlphaTub_NEURO, axis=0)
ax6.plot(RADIUS, MEAN)
ax6.fill_between(RADIUS, MEAN+SEM, MEAN-SEM, alpha=0.1)


ax1.set_title('NeuN Astro')
ax2.set_title('MAP2 Astro')
ax3.set_title('Alpha-Tubulin Astro')
ax4.set_title('NeuN Neuro')
ax5.set_title('MAP2 Neuro')
ax6.set_title('Alpha-Tubulin Neuro')
plt.tight_layout()


plt.figure()
#plt.hist(CELL_AREA, histtype='stepfilled', color='black', alpha=0.5)
plt.hist(AREA_ASTRO, histtype='step', cumulative=True, density=True, label='Non-neuro')
plt.hist(AREA_NEURO, histtype='step', cumulative=True, density=True, label='Neuro')
plt.xlabel('Cell area(um2)')
plt.ylabel('Cumulative fraction')

plt.figure()
plt.scatter(NeuNIntensity, CELL_AREA)
plt.xlabel('Normalized NeuN')
plt.ylabel('Cell area (um2)')
