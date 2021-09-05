# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 12:00:38 2021

@author: KEVIN-DORGANS
"""
import scipy as sp

FULL_PROMOTER_LIST = []

FULL_eGFP_INTENSITY_normed_ASTRO = []
FULL_eGFP_INTENSITY_ASTRO = []

FULL_eGFP_INTENSITY_normed = []
FULL_eGFP_INTENSITY = []
FULL_AAVPHPS_INTENSITY_normed = []
FULL_AAVPHPS_INTENSITY = []

FULL_AAVPHPS_INTENSITY_normed_ASTRO = []
FULL_AAVPHPS_INTENSITY_ASTRO  = []

FULL_DIVIDED_normed_ASTRO = []
FULL_DIVIDED_ASTRO = []
FULL_DIVIDED_normed = []
FULL_DIVIDED = []

FULL_AREAS = []

plt.figure(figsize=(5,5))
ax1 = plt.subplot(321)
ax2 = plt.subplot(322)
ax3 = plt.subplot(323)
ax4 = plt.subplot(324)
ax5 = plt.subplot(325)
ax6 = plt.subplot(326)

PATHS, IMG_INDEX = load_directory_content_and_sub__()
IMG_INDEX = IMG_INDEX + r'\\IMG_INDEX.csv'
IMG_INDEX = pd.read_csv(IMG_INDEX)

for i in range(len(PATHS)):
    if ('csv' in PATHS[i])==True and ('IMG_INDEX' in PATHS[i])==False:
        PATH = PATHS[i]
        PROMOTER = PATH.split('\\')[-1].split('_')[1]

        CSV_FILE = pd.read_csv(PATH)
        eGFP = CSV_FILE['eGFP']
        tdTomato = CSV_FILE['tdTomato']
        area = CSV_FILE['AREA']
        Divided = np.divide(eGFP, tdTomato)
        zscore_eGFP = sp.stats.zscore(eGFP)
        zscore_tdTomato = sp.stats.zscore(tdTomato)
        zscore_divided = sp.stats.zscore(Divided)
        
        tdTomato_ASTRO = [tdTomato[j] for j in range(len(tdTomato)) if area[j]<100]
        eGFP_ASTRO = [eGFP[j] for j in range(len(tdTomato)) if area[j]<100]
        tdTomato_NEURO = [tdTomato[j] for j in range(len(tdTomato)) if area[j]>100]
        eGFP_NEURO = [eGFP[j] for j in range(len(tdTomato)) if area[j]>100]
        Divided_Astro = [Divided[j] for j in range(len(tdTomato)) if area[j]<100]
        Divided_Neuro = [Divided[j] for j in range(len(tdTomato)) if area[j]>100]

       
        
        tdTomato = tdTomato/np.nanmax(tdTomato)
        eGFP = eGFP/np.nanmax(eGFP)
        Divided_normed = np.divide(eGFP, tdTomato)
        
        tdTomato_ASTRO_normed = [tdTomato[j] for j in range(len(tdTomato)) if area[j]<100]
        eGFP_ASTRO_normed = [eGFP[j] for j in range(len(tdTomato)) if area[j]<100]
        tdTomato_NEURO_normed = [tdTomato[j] for j in range(len(tdTomato)) if area[j]>100]
        eGFP_NEURO_normed = [eGFP[j] for j in range(len(tdTomato)) if area[j]>100]
        Divided_ASTRO_normed = [Divided_normed[j] for j in range(len(tdTomato)) if area[j]<100]
        Divided_NEURO_normed = [Divided_normed[j] for j in range(len(tdTomato)) if area[j]>100]
        
        

        ax1.hist(tdTomato_NEURO_normed, color='black', alpha=0.1, bins=10)
        #ax5.hist(eGFP_NEURO_normed, color='green', alpha=0.1, bins=10)
        ax3.scatter(PROMOTER, np.nanmean(eGFP_NEURO), facecolor='white', edgecolor='purple')
        
        eGFP_ASTRO_NEURO_RATIO = np.nanmean(eGFP_NEURO_normed)/np.nanmean(eGFP_ASTRO_normed)
        tdTomato_ASTRO_NEURO_RATIO = np.nanmean(tdTomato_NEURO_normed)/np.nanmean(tdTomato_ASTRO_normed)

        ax2.scatter(PROMOTER, eGFP_ASTRO_NEURO_RATIO, facecolor='white', edgecolor='purple')
        #ax2.scatter(PROMOTER, tdTomato_ASTRO_NEURO_RATIO, color='red')

        SPARSE_CELLS = [CSV_FILE['eGFP'][j] for j in range(len(zscore_divided)) if zscore_divided[j]>3]
        OUTLINER_CELLS_eGFP = [CSV_FILE['eGFP'][j] for j in range(len(zscore_eGFP)) if zscore_eGFP[j]>3]
        OUTLINER_CELLS_tdTomato = [CSV_FILE['tdTomato'][j] for j in range(len(zscore_tdTomato)) if zscore_tdTomato[j]>3]
        
        if len(OUTLINER_CELLS_eGFP)>0:
            for j in range(len(OUTLINER_CELLS_eGFP)):
                ax4.scatter(PROMOTER, OUTLINER_CELLS_eGFP[j]/np.nanmean(eGFP_NEURO), facecolor='white', edgecolor='purple')
                ax5.scatter(PROMOTER, OUTLINER_CELLS_eGFP[j], facecolor='white', edgecolor='purple')

        else:
            ax4.scatter(PROMOTER, 0, color='black', alpha=0.1)
            ax5.scatter(PROMOTER, 0, color='black', alpha=0.1)
        for j in range(len(IMG_INDEX['FILENAME'])):
            if (IMG_INDEX['FILENAME'][j] in PATH )==True:
                IO_AREA = IMG_INDEX['IOAREA'][j]
                
        #ax5.scatter(PROMOTER, len(OUTLINER_CELLS_eGFP)/IO_AREA, facecolor='white', edgecolor='purple')
        print(PROMOTER+': '+str(len(SPARSE_CELLS))+' : '+str(len(OUTLINER_CELLS_eGFP))+' : '+str(len(OUTLINER_CELLS_tdTomato)))
        a, b = sp.stats.normaltest(eGFP)
        ax6.scatter(PROMOTER, b, facecolor='white', edgecolor='purple')
        
        FULL_PROMOTER_LIST.append(PROMOTER)
        FULL_eGFP_INTENSITY_normed.append(eGFP_NEURO_normed)
        FULL_eGFP_INTENSITY.append(eGFP_NEURO)
        FULL_eGFP_INTENSITY_normed_ASTRO.append(eGFP_ASTRO_normed)
        FULL_eGFP_INTENSITY_ASTRO.append(eGFP_ASTRO)
        
        FULL_AAVPHPS_INTENSITY_normed.append(tdTomato_NEURO_normed)
        FULL_AAVPHPS_INTENSITY.append(tdTomato_NEURO)
        FULL_AAVPHPS_INTENSITY_normed_ASTRO.append(tdTomato_ASTRO_normed)
        FULL_AAVPHPS_INTENSITY_ASTRO.append(tdTomato_ASTRO)
        
        FULL_DIVIDED_normed_ASTRO.append(Divided_ASTRO_normed)
        FULL_DIVIDED_ASTRO.append(Divided_Astro)
        FULL_DIVIDED_normed.append(Divided_NEURO_normed)
        FULL_DIVIDED.append(Divided_Neuro)
        
        FULL_AREAS.append(area)
        
ax1.set_title('AAV9 / AAV.PHP.S')
ax2.set_title('eGFP-Neuro/eGFP-Astro ')
ax3.set_title('eGFP-Neuro ')
ax4.set_title('Sparse expression (%mean)')
ax5.set_title('Strong cells/mm2')
ax6.set_title('Normality test (p-value)')
ax6.set_yscale('log')

LIST_ = [



'Isgf9(1.3)','Isgf9(2.5)','Isgf9(3.7)','5HTr2b(1.0)','5HTr2b(1.8)','5HTr2b(3.0)','5HTr2b(3.7)','PDX1','SUSD4','CAG']

plt.figure(figsize=(1,len(LIST_)))

for j in range(len(LIST_)):
    ax = plt.subplot(len(LIST_), 1, j+1)
    temp_ = []
    for i in range(len(FULL_PROMOTER_LIST)):
        if FULL_PROMOTER_LIST[i]==LIST_[j]:
            temp_.append(FULL_eGFP_INTENSITY_normed[i])
    ax.hist(np.concatenate(temp_), histtype='step', bins=10)

    ax.set_title(LIST_[j])
    
    
plt.figure()
ax = plt.subplot(231)
ax2 = plt.subplot(234, sharex=ax)
ax3 = plt.subplot(235)
ax4 = plt.subplot(236)
ax5 = plt.subplot(232, sharex=ax3)

for j in range(len(LIST_)):
    temp_ = []
    temp_2 = []
    ratio_ = []
    for i in range(len(FULL_PROMOTER_LIST)):
        if FULL_PROMOTER_LIST[i]==LIST_[j]:
            ratio_.append(np.nanmean(FULL_DIVIDED[i]))
            temp_.append(FULL_DIVIDED[i])
    ax.hist(ratio_, histtype='step', density=True)
    MEAN = (np.nanmean(ratio_)-np.nanmean(np.concatenate(FULL_DIVIDED)))*100
    SEM = sp.stats.sem(ratio_)
    try:
        print(LIST_[j]+' '+str(MEAN)+' +/-'+str(SEM)+' n='+str(len(np.concatenate(temp_))))
    except:
        print('No vals. for'+str(LIST_[j]))
    ax2.scatter( MEAN, LIST_[j])
    ax2.plot((MEAN+SEM, MEAN-SEM), (LIST_[j], LIST_[j]))
ax.hist(np.concatenate(FULL_DIVIDED), histtype='step', bins=100, density=True)
ax5.hist([np.nanmean(FULL_eGFP_INTENSITY[j])/np.nanmean(FULL_eGFP_INTENSITY_ASTRO[j]) for j in range(len(FULL_eGFP_INTENSITY))], histtype='step',  density=True)

print('NEURON_SPE')
for j in range(len(LIST_)):
    temp_ = []
    temp_2 = []
    temp_3 = []
    ratio_ = []
    for i in range(len(FULL_PROMOTER_LIST)):
        if FULL_PROMOTER_LIST[i]==LIST_[j]:
            temp_.append(np.nanmean(FULL_eGFP_INTENSITY[i]))
            temp_2.append(np.nanmean(FULL_eGFP_INTENSITY_ASTRO[i]))
            temp_3.append(np.nanmean(FULL_DIVIDED[i]))
            ratio_.append(np.nanmean(FULL_eGFP_INTENSITY[i])/np.nanmean(FULL_eGFP_INTENSITY_ASTRO[i]))
    MEAN = (np.nanmean(ratio_)-1)*100
    SEM = sp.stats.sem(ratio_)
    try:
        print(LIST_[j]+' '+str(MEAN)+' +/-'+str(SEM))
    except:
        print('No vals. for'+str(LIST_[j]))
    ax3.scatter(MEAN, LIST_[j])
    ax3.plot((MEAN+SEM, MEAN-SEM), (LIST_[j], LIST_[j]))
    ax4.scatter(np.nanmean(ratio_), np.nanmean(temp_3))
temp_ = []
temp_2 = []
ratio_ = []
for i in range(len(FULL_AAVPHPS_INTENSITY)):
    temp_.append(np.nanmean(FULL_AAVPHPS_INTENSITY[i]))
    temp_2.append(np.nanmean(FULL_AAVPHPS_INTENSITY_ASTRO[i]))
    ratio_.append(np.nanmean(FULL_AAVPHPS_INTENSITY[i])/np.nanmean(FULL_AAVPHPS_INTENSITY_ASTRO[i]))
MEAN = np.nanmean(ratio_)
SEM = sp.stats.sem(ratio_)
ax3.scatter(MEAN, 'AAV.PHP.S-CAG')
ax3.plot((MEAN+SEM, MEAN-SEM), ('AAV.PHP.S-CAG', 'AAV.PHP.S-CAG'))

ax2.set_xlabel('Normalized expression')
ax3.set_xlabel('Neuron-specificity')
ax4.set_xlabel('Neuron-specificity')
ax4.set_ylabel('Normalized expression')

plt.tight_layout()