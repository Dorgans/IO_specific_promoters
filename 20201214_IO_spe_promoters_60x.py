# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 22:13:44 2020

@author: dorga
"""

#LOADS CSVS

PATHS = load_directory_content__()

ALL_AREAS = []

FULL_tDt_CELLS = []
FULL_tDt_AREA = []
FULL_Igsf9_25_CELLS = []
FULL_Igsf9_25_AREA = []
FULL_Igsf9_13_CELLS = []
FULL_Igsf9_13_AREA = []
FULL_Igsf9_37_CELLS = []
FULL_Igsf9_37_AREA = []
FULL_5HTr2b_CELLS = []
FULL_5HTr2b_AREA = []
FULL_5HTr2b_10_CELLS = []
FULL_5HTr2b_10_AREA = []
FULL_5HTr2b_18_CELLS = []
FULL_5HTr2b_18_AREA = []
FULL_5HTr2b_30_CELLS = []
FULL_5HTr2b_30_AREA = []
FULL_PDX1_CELLS = []
FULL_PDX1_AREA = []
FULL_SUSD4_CELLS = []
FULL_SUSD4_AREA = []

FULL_Igsf9_25_CELLS_tDt = []
FULL_Igsf9_37_CELLS_tDt = []
FULL_Igsf9_13_CELLS_tDt = []
FULL_5HTr2b_CELLS_tDt = []
FULL_5HTr2b_10_CELLS_tDt = []
FULL_5HTr2b_18_CELLS_tDt = []
FULL_5HTr2b_30_CELLS_tDt = []
FULL_PDX1_CELLS_tDt = []
FULL_SUSD4_CELLS_tDt = []

FULL_Igsf9_37_DIVIDED = []
FULL_Igsf9_25_DIVIDED = []
FULL_Igsf9_13_DIVIDED = []
FULL_5HTr2b_10_DIVIDED = []
FULL_5HTr2b_18_DIVIDED = []
FULL_5HTr2b_30_DIVIDED = []
FULL_5HTr2b_DIVIDED = []
FULL_PDX1_DIVIDED = []
FULL_SUSD4_DIVIDED = []
FULL_Igsf9_37_SUBTRACTED = []
FULL_Igsf9_25_SUBTRACTED = []
FULL_Igsf9_13_SUBTRACTED = []
FULL_5HTr2b_10_SUBTRACTED = []
FULL_5HTr2b_18_SUBTRACTED = []
FULL_5HTr2b_30_SUBTRACTED = []
FULL_5HTr2b_SUBTRACTED = []
FULL_PDX1_SUBTRACTED = []
FULL_SUSD4_SUBTRACTED = []

plt.figure(figsize=(14,7))
ax = plt.subplot(461)
ax2 = plt.subplot(462)
ax3 = plt.subplot(463, sharex=ax2, sharey=ax2)
ax4 = plt.subplot(464, sharex=ax2, sharey=ax2)
ax5 = plt.subplot(465, sharex=ax2, sharey=ax2)
ax6 = plt.subplot(466, sharex=ax2, sharey=ax2)

ax7 = plt.subplot(467)
ax8 = plt.subplot(468)
ax9 = plt.subplot(469, sharex=ax8, sharey=ax8)
ax10 = plt.subplot(4, 6, 10, sharex=ax8, sharey=ax8)
ax11 = plt.subplot(4,6 ,11, sharex=ax8, sharey=ax8)
ax12 = plt.subplot(4, 6, 12, sharex=ax8, sharey=ax8)

ax13 = plt.subplot(4, 6, 13)
ax14 = plt.subplot(4, 6, 14)
ax15 = plt.subplot(4, 6, 15, sharex=ax14, sharey=ax14)
ax16 = plt.subplot(4, 6, 16, sharex=ax14, sharey=ax14)
ax17 = plt.subplot(4, 6, 17, sharex=ax14, sharey=ax14)
ax18 = plt.subplot(4, 6, 18, sharex=ax14, sharey=ax14)

ax19 = plt.subplot(4, 6, 19)
ax20 = plt.subplot(4, 6, 20)
ax21 = plt.subplot(4, 6, 21, sharex=ax20, sharey=ax20)
ax22 = plt.subplot(4, 6, 22, sharex=ax20, sharey=ax20)
ax23 = plt.subplot(4, 6, 23, sharex=ax20, sharey=ax20)
ax24 = plt.subplot(4, 6, 24, sharex=ax20, sharey=ax20)

for i in range(len(PATHS[0])):
    file = PATHS[0][i]
    
    if ('KD' in file)==True and ('0x' in file)==True and ('ROI' in file)==False and ('csv' in file)==True:
        print('OK')
        Cell_file = pd.read_csv(file)
        CELL_AREA = Cell_file.values[:,1]
        eGFP_intensity = Cell_file.values[:,3]
        tDt_intensity = Cell_file.values[:,2]
        SUB_CHANNELS = Cell_file.values[:,5]-Cell_file.values[:,4]
        SUB_CHANNELS_2 = Cell_file.values[:,4]
        
        FULL_tDt_CELLS.append(tDt_intensity)
        FULL_tDt_AREA.append(CELL_AREA)
        
        if ('SUSD4' in file)==True:
            FULL_SUSD4_CELLS.append(eGFP_intensity)
            FULL_SUSD4_CELLS_tDt.append(tDt_intensity)
            FULL_SUSD4_AREA.append(CELL_AREA)
            FULL_SUSD4_SUBTRACTED.append(np.subtract(SUB_CHANNELS, SUB_CHANNELS_2))
            FULL_SUSD4_DIVIDED.append(np.divide(eGFP_intensity, tDt_intensity))
            #ax6.scatter(tDt_intensity, eGFP_intensity)
        elif ('PDX1' in file)==True:
            FULL_PDX1_CELLS.append(eGFP_intensity)
            FULL_PDX1_CELLS_tDt.append(tDt_intensity)
            FULL_PDX1_AREA.append(CELL_AREA)
            FULL_PDX1_SUBTRACTED.append(np.subtract(SUB_CHANNELS, SUB_CHANNELS_2))
            FULL_PDX1_DIVIDED.append(np.divide(eGFP_intensity, tDt_intensity))
            #ax5.scatter(tDt_intensity, eGFP_intensity)
        elif ('5HTr2b' in file)==True and ('1.0' in file)==True:
            FULL_5HTr2b_10_CELLS.append(eGFP_intensity)
            FULL_5HTr2b_10_CELLS_tDt.append(tDt_intensity)
            FULL_5HTr2b_10_AREA.append(CELL_AREA)
            FULL_5HTr2b_10_SUBTRACTED.append(np.subtract(SUB_CHANNELS, SUB_CHANNELS_2))
            FULL_5HTr2b_10_DIVIDED.append(np.divide(eGFP_intensity, tDt_intensity))
            #ax4.scatter(tDt_intensity, eGFP_intensity)
        elif ('5HTr2b' in file)==True and ('1.8' in file)==True:
            FULL_5HTr2b_18_CELLS.append(eGFP_intensity)
            FULL_5HTr2b_18_CELLS_tDt.append(tDt_intensity)
            FULL_5HTr2b_18_AREA.append(CELL_AREA)
            FULL_5HTr2b_18_SUBTRACTED.append(np.subtract(SUB_CHANNELS, SUB_CHANNELS_2))
            FULL_5HTr2b_18_DIVIDED.append(np.divide(eGFP_intensity, tDt_intensity))
            #ax4.scatter(tDt_intensity, eGFP_intensity)
        elif ('5HTr2b' in file)==True and ('3.0' in file)==True:
            FULL_5HTr2b_30_CELLS.append(eGFP_intensity)
            FULL_5HTr2b_30_CELLS_tDt.append(tDt_intensity)
            FULL_5HTr2b_30_AREA.append(CELL_AREA)
            FULL_5HTr2b_30_SUBTRACTED.append(np.subtract(SUB_CHANNELS, SUB_CHANNELS_2))
            FULL_5HTr2b_30_DIVIDED.append(np.divide(eGFP_intensity, tDt_intensity))
            #ax4.scatter(tDt_intensity, eGFP_intensity)
        elif ('5HTr2b' in file)==True:# and ('3.7' in file)==True:
            FULL_5HTr2b_CELLS.append(eGFP_intensity)
            FULL_5HTr2b_CELLS_tDt.append(tDt_intensity)
            FULL_5HTr2b_AREA.append(CELL_AREA)
            FULL_5HTr2b_SUBTRACTED.append(np.subtract(SUB_CHANNELS, SUB_CHANNELS_2))
            FULL_5HTr2b_DIVIDED.append(np.divide(eGFP_intensity, tDt_intensity))
            #ax4.scatter(tDt_intensity, eGFP_intensity)
        elif ('Igsf9' in file)==True and ('1.3' in file)==True:
            FULL_Igsf9_13_CELLS.append(eGFP_intensity)
            FULL_Igsf9_13_CELLS_tDt.append(tDt_intensity)
            FULL_Igsf9_13_AREA.append(CELL_AREA)
            FULL_Igsf9_13_SUBTRACTED.append(np.subtract(SUB_CHANNELS, SUB_CHANNELS_2))
            FULL_Igsf9_13_DIVIDED.append(np.divide(eGFP_intensity, tDt_intensity))
            #ax3.scatter(tDt_intensity, eGFP_intensity)
        elif ('Igsf9' in file)==True and ('2.5' in file)==True:
            FULL_Igsf9_25_CELLS.append(eGFP_intensity)
            FULL_Igsf9_25_CELLS_tDt.append(tDt_intensity)
            FULL_Igsf9_25_AREA.append(CELL_AREA)
            FULL_Igsf9_25_SUBTRACTED.append(np.subtract(SUB_CHANNELS, SUB_CHANNELS_2))
            FULL_Igsf9_25_DIVIDED.append(np.divide(eGFP_intensity, tDt_intensity))
            #ax2.scatter(tDt_intensity, eGFP_intensity)
        elif ('Igsf9' in file)==True and ('3.7' in file)==True:
            FULL_Igsf9_37_CELLS.append(eGFP_intensity)
            FULL_Igsf9_37_CELLS_tDt.append(tDt_intensity)
            FULL_Igsf9_37_AREA.append(CELL_AREA)
            FULL_Igsf9_37_SUBTRACTED.append(np.subtract(SUB_CHANNELS, SUB_CHANNELS_2))
            FULL_Igsf9_37_DIVIDED.append(np.divide(eGFP_intensity, tDt_intensity))
            #ax2.scatter(tDt_intensity, eGFP_intensity)
            
            
            
#PLOTS FOR Igsf9 group
            
plt.figure(figsize=(10, 4), num = 'Igsf9 group')
ax1 = plt.subplot(2, 6, 1)
ax2 = plt.subplot(2, 6, 2)
ax3 = plt.subplot(2, 6, 3)
ax4 = plt.subplot(2, 6, 4)
ax5 = plt.subplot(2, 6, 5)
ax6 = plt.subplot(2, 6, 6)
ax7 = plt.subplot(2, 6, 7)
ax8 = plt.subplot(2, 6, 8)
ax9 = plt.subplot(2, 6, 9)
ax10 = plt.subplot(2, 6, 10)
ax11 = plt.subplot(2, 6, 11)
ax12 = plt.subplot(2, 6, 12)
j=0

CELLS = ([FULL_Igsf9_13_CELLS[j][i] for i in range(len(FULL_Igsf9_13_CELLS[j])) if FULL_Igsf9_13_AREA[j][i]>100])
MEAN_Igsf9_13_NEURO = CELLS

CELLS_1 = ([FULL_Igsf9_25_CELLS[j][i] for i in range(len(FULL_Igsf9_25_CELLS[j])) if FULL_Igsf9_25_AREA[j][i]>100])
MEAN_Igsf9_25_NEURO = CELLS_1

CELLS_2 = ([FULL_Igsf9_37_CELLS[j][i] for i in range(len(FULL_Igsf9_37_CELLS[j])) if FULL_Igsf9_37_AREA[j][i]>100])
MEAN_Igsf9_37_NEURO = CELLS_2

ax1.boxplot([CELLS, CELLS_1, CELLS_2])


CELLS = ([FULL_Igsf9_13_CELLS_tDt[j][i] for i in range(len(FULL_Igsf9_13_CELLS[j])) if FULL_Igsf9_13_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
MEAN_Igsf9_13_NEURO_tDt = CELLS

CELLS_1 = ([FULL_Igsf9_25_CELLS_tDt[j][i] for i in range(len(FULL_Igsf9_25_CELLS[j])) if FULL_Igsf9_25_AREA[j][i]>100])
MEAN = np.nanmean(CELLS_1)
MEAN_Igsf9_25_NEURO_tDt = CELLS_1

CELLS_2 = ([FULL_Igsf9_37_CELLS_tDt[j][i] for i in range(len(FULL_Igsf9_37_CELLS[j])) if FULL_Igsf9_37_AREA[j][i]>100])
MEAN = np.nanmean(CELLS_2)
MEAN_Igsf9_37_NEURO_tDt = CELLS_2

ax2.boxplot([CELLS, CELLS_1, CELLS_2])


CELLS = ([FULL_Igsf9_13_SUBTRACTED[j][i] for i in range(len(FULL_Igsf9_13_CELLS[j])) if FULL_Igsf9_13_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
MEAN_Igsf9_13_SUBTRACTED_NEURO = CELLS

CELLS_1 = ([FULL_Igsf9_25_SUBTRACTED[j][i] for i in range(len(FULL_Igsf9_25_CELLS[j])) if FULL_Igsf9_25_AREA[j][i]>100])
MEAN = np.nanmean(CELLS_1)
MEAN_Igsf9_25_SUBTRACTED_NEURO = CELLS_1

CELLS_2 = ([FULL_Igsf9_37_SUBTRACTED[j][i] for i in range(len(FULL_Igsf9_37_CELLS[j])) if FULL_Igsf9_37_AREA[j][i]>100])
MEAN = np.nanmean(CELLS_2)
MEAN_Igsf9_37_SUBTRACTED_NEURO = CELLS_2

ax3.boxplot([CELLS, CELLS_1, CELLS_2])



CELLS = ([FULL_Igsf9_13_DIVIDED[j][i] for i in range(len(FULL_Igsf9_13_CELLS[j])) if FULL_Igsf9_13_AREA[j][i]>100])
MEAN_Igsf9_13_DIVIDED_NEURO = CELLS

CELLS_1 = ([FULL_Igsf9_25_DIVIDED[j][i] for i in range(len(FULL_Igsf9_25_CELLS[j])) if FULL_Igsf9_25_AREA[j][i]>100])
MEAN_Igsf9_25_DIVIDED_NEURO= CELLS_1

CELLS_2 = ([FULL_Igsf9_37_DIVIDED[j][i] for i in range(len(FULL_Igsf9_37_CELLS[j])) if FULL_Igsf9_37_AREA[j][i]>100])
MEAN_Igsf9_37_DIVIDED_NEURO = CELLS_2

ax4.boxplot([CELLS, CELLS_1, CELLS_2])



CELLS = ([FULL_Igsf9_13_CELLS[j][i] for i in range(len(FULL_Igsf9_13_CELLS[j])) if FULL_Igsf9_13_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_Igsf9_13_CELLS_tDt[j][i] for i in range(len(FULL_Igsf9_13_CELLS[j])) if FULL_Igsf9_13_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
ax5.scatter(CELLS, CELLS_tDt)
ax6.hist(CELLS, alpha=0.7)
ax12.hist(CELLS_tDt, alpha=0.7, color='red')

CELLS = ([FULL_Igsf9_25_CELLS[j][i] for i in range(len(FULL_Igsf9_25_CELLS[j])) if FULL_Igsf9_25_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_Igsf9_25_CELLS_tDt[j][i] for i in range(len(FULL_Igsf9_25_CELLS[j])) if FULL_Igsf9_25_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
ax5.scatter(CELLS, CELLS_tDt)
ax6.hist(CELLS, alpha=0.7)
ax12.hist(CELLS_tDt, alpha=0.7, color='red')

CELLS = ([FULL_Igsf9_37_CELLS[j][i] for i in range(len(FULL_Igsf9_37_CELLS[j])) if FULL_Igsf9_37_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_Igsf9_37_CELLS_tDt[j][i] for i in range(len(FULL_Igsf9_37_CELLS[j])) if FULL_Igsf9_37_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
ax5.scatter(CELLS, CELLS_tDt)
ax6.hist(CELLS, alpha=0.7)
ax12.hist(CELLS_tDt, alpha=0.7, color='red')

ax1.set_title('eGFP')
ax2.set_title('tD-Tomato')
ax3.set_title('Equalized Sub.')
ax4.set_title('eGFP/tD-tomato')
ax5.set_xlabel('eGFP')
ax5.set_ylabel('tD-Tomato')


CELLS = ([FULL_Igsf9_13_CELLS[j][i] for i in range(len(FULL_Igsf9_13_CELLS[j])) if FULL_Igsf9_13_AREA[j][i]<100])
MEAN_Igsf9_13_ASTRO = CELLS

CELLS_1 = ([FULL_Igsf9_25_CELLS[j][i] for i in range(len(FULL_Igsf9_25_CELLS[j])) if FULL_Igsf9_25_AREA[j][i]<100])
MEAN_Igsf9_25_ASTRO = CELLS_1

CELLS_2 = ([FULL_Igsf9_37_CELLS[j][i] for i in range(len(FULL_Igsf9_37_CELLS[j])) if FULL_Igsf9_37_AREA[j][i]<100])
MEAN_Igsf9_37_ASTRO= CELLS_2

ax7.boxplot([CELLS, CELLS_1, CELLS_2])


CELLS = ([FULL_Igsf9_13_CELLS_tDt[j][i] for i in range(len(FULL_Igsf9_13_CELLS[j])) if FULL_Igsf9_13_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_1 = ([FULL_Igsf9_25_CELLS_tDt[j][i] for i in range(len(FULL_Igsf9_25_CELLS[j])) if FULL_Igsf9_25_AREA[j][i]<100])
MEAN = np.nanmean(CELLS_1)

CELLS_2 = ([FULL_Igsf9_37_CELLS_tDt[j][i] for i in range(len(FULL_Igsf9_37_CELLS[j])) if FULL_Igsf9_37_AREA[j][i]<100])
MEAN = np.nanmean(CELLS_2)


ax8.boxplot([CELLS, CELLS_1, CELLS_2])


CELLS = ([FULL_Igsf9_13_SUBTRACTED[j][i] for i in range(len(FULL_Igsf9_13_CELLS[j])) if FULL_Igsf9_13_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_1 = ([FULL_Igsf9_25_SUBTRACTED[j][i] for i in range(len(FULL_Igsf9_25_CELLS[j])) if FULL_Igsf9_25_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_2 = ([FULL_Igsf9_37_SUBTRACTED[j][i] for i in range(len(FULL_Igsf9_37_CELLS[j])) if FULL_Igsf9_37_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

ax9.boxplot([CELLS, CELLS_1, CELLS_2])



CELLS = ([FULL_Igsf9_13_DIVIDED[j][i] for i in range(len(FULL_Igsf9_13_CELLS[j])) if FULL_Igsf9_13_AREA[j][i]<100])
MEAN_Igsf9_13_DIVIDED_ASTRO = CELLS

CELLS_1 = ([FULL_Igsf9_25_DIVIDED[j][i] for i in range(len(FULL_Igsf9_25_CELLS[j])) if FULL_Igsf9_25_AREA[j][i]<100])
MEAN_Igsf9_25_DIVIDED_ASTRO = CELLS_1

CELLS_2 = ([FULL_Igsf9_37_DIVIDED[j][i] for i in range(len(FULL_Igsf9_37_CELLS[j])) if FULL_Igsf9_37_AREA[j][i]<100])
MEAN_Igsf9_37_DIVIDED_ASTRO= CELLS_2

ax10.boxplot([CELLS, CELLS_1, CELLS_2])



CELLS = ([FULL_Igsf9_13_CELLS[j][i] for i in range(len(FULL_Igsf9_13_CELLS[j])) if FULL_Igsf9_13_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_Igsf9_13_CELLS_tDt[j][i] for i in range(len(FULL_Igsf9_13_CELLS[j])) if FULL_Igsf9_13_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)
ax11.scatter(CELLS, CELLS_tDt)

CELLS = ([FULL_Igsf9_25_CELLS[j][i] for i in range(len(FULL_Igsf9_25_CELLS[j])) if FULL_Igsf9_25_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_Igsf9_25_CELLS_tDt[j][i] for i in range(len(FULL_Igsf9_25_CELLS[j])) if FULL_Igsf9_25_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)
ax11.scatter(CELLS, CELLS_tDt)

CELLS = ([FULL_Igsf9_37_CELLS[j][i] for i in range(len(FULL_Igsf9_37_CELLS[j])) if FULL_Igsf9_37_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_Igsf9_37_CELLS_tDt[j][i] for i in range(len(FULL_Igsf9_37_CELLS[j])) if FULL_Igsf9_37_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)
ax11.scatter(CELLS, CELLS_tDt)

ax7.set_title('eGFP Astro.')
ax8.set_title('tD-Tomato Astro.')
ax9.set_title('Equalized Sub. Astro.')
ax10.set_title('eGFP/tD-tomato Astro.')
ax11.set_xlabel('eGFP')
ax11.set_ylabel('tD-Tomato')

plt.tight_layout()



#Plots for 5HTr2b group

plt.figure(figsize=(10, 4), num='5Htr2b Group')
ax1 = plt.subplot(2, 6, 1)
ax2 = plt.subplot(2, 6, 2)
ax3 = plt.subplot(2, 6, 3)
ax4 = plt.subplot(2, 6, 4)
ax5 = plt.subplot(2, 6, 5)
ax6 = plt.subplot(2, 6, 6)
ax7 = plt.subplot(2, 6, 7)
ax8 = plt.subplot(2, 6, 8)
ax9 = plt.subplot(2, 6, 9)
ax10 = plt.subplot(2, 6, 10)
ax11 = plt.subplot(2, 6, 11)
ax12 = plt.subplot(2, 6, 12)

j=0

CELLS = ([FULL_5HTr2b_10_CELLS[j][i] for i in range(len(FULL_5HTr2b_10_CELLS[j])) if FULL_5HTr2b_10_AREA[j][i]>100])
MEAN_5HTR2b_10_NEURO = CELLS

CELLS_1 = ([FULL_5HTr2b_18_CELLS[j][i] for i in range(len(FULL_5HTr2b_18_CELLS[j])) if FULL_5HTr2b_18_AREA[j][i]>100])
MEAN_5HTR2b_18_NEURO = CELLS_1

CELLS_2 = ([FULL_5HTr2b_30_CELLS[j][i] for i in range(len(FULL_5HTr2b_30_CELLS[j])) if FULL_5HTr2b_30_AREA[j][i]>100])
MEAN_5HTR2b_30_NEURO = CELLS_2

CELLS_3 = ([FULL_5HTr2b_CELLS[j][i] for i in range(len(FULL_5HTr2b_CELLS[j])) if FULL_5HTr2b_AREA[j][i]>100])
MEAN_5HTR2b_37_NEURO = CELLS_3
ax1.boxplot([CELLS, CELLS_1, CELLS_2, CELLS_3])


CELLS = ([FULL_5HTr2b_10_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_10_CELLS[j])) if FULL_5HTr2b_10_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
MEAN_5HTR2b_10_NEURO_tDt = CELLS

CELLS_1 = ([FULL_5HTr2b_18_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_18_CELLS[j])) if FULL_5HTr2b_18_AREA[j][i]>100])
MEAN = np.nanmean(CELLS_1)
MEAN_5HTR2b_18_NEURO_tDt = CELLS_1

CELLS_2 = ([FULL_5HTr2b_30_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_30_CELLS[j])) if FULL_5HTr2b_30_AREA[j][i]>100])
MEAN = np.nanmean(CELLS_2)
MEAN_5HTR2b_30_NEURO_tDt = CELLS_2

CELLS_3 = ([FULL_5HTr2b_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_CELLS[j])) if FULL_5HTr2b_AREA[j][i]>100])
MEAN = np.nanmean(CELLS_3)
MEAN_5HTR2b_37_NEURO_tDt = CELLS_3

ax2.boxplot([CELLS, CELLS_1, CELLS_2, CELLS_3])


CELLS = ([FULL_5HTr2b_10_SUBTRACTED[j][i] for i in range(len(FULL_5HTr2b_10_CELLS[j])) if FULL_5HTr2b_10_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
MEAN_5HTR2b_10_SUBTRACTED_NEURO = CELLS

CELLS_1 = ([FULL_5HTr2b_18_SUBTRACTED[j][i] for i in range(len(FULL_5HTr2b_18_CELLS[j])) if FULL_5HTr2b_18_AREA[j][i]>100])
MEAN = np.nanmean(CELLS_1)
MEAN_5HTR2b_18_SUBTRACTED_NEURO = CELLS_1

CELLS_2 = ([FULL_5HTr2b_30_SUBTRACTED[j][i] for i in range(len(FULL_5HTr2b_30_CELLS[j])) if FULL_5HTr2b_30_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
MEAN_5HTR2b_30_SUBTRACTED_NEURO = CELLS_2

CELLS_3 = ([FULL_5HTr2b_SUBTRACTED[j][i] for i in range(len(FULL_5HTr2b_CELLS[j])) if FULL_5HTr2b_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
MEAN_5HTR2b_37_SUBTRACTED_NEURO = CELLS_3

ax3.boxplot([CELLS, CELLS_1, CELLS_2, CELLS_3])



CELLS = ([FULL_5HTr2b_10_DIVIDED[j][i] for i in range(len(FULL_5HTr2b_10_CELLS[j])) if FULL_5HTr2b_10_AREA[j][i]>100])
MEAN_5HTR2b_10_DIVIDED_NEURO  = CELLS

CELLS_1 = ([FULL_5HTr2b_18_DIVIDED[j][i] for i in range(len(FULL_5HTr2b_18_CELLS[j])) if FULL_5HTr2b_18_AREA[j][i]>100])
MEAN_5HTR2b_18_DIVIDED_NEURO  = CELLS_1

CELLS_2 = ([FULL_5HTr2b_30_DIVIDED[j][i] for i in range(len(FULL_5HTr2b_30_CELLS[j])) if FULL_5HTr2b_30_AREA[j][i]>100])
MEAN_5HTR2b_30_DIVIDED_NEURO  = CELLS_2

CELLS_3 = ([FULL_5HTr2b_DIVIDED[j][i] for i in range(len(FULL_5HTr2b_CELLS[j])) if FULL_5HTr2b_AREA[j][i]>100])
MEAN_5HTR2b_37_DIVIDED_NEURO  = CELLS_3
ax4.boxplot([CELLS, CELLS_1, CELLS_2, CELLS_3])



CELLS = ([FULL_5HTr2b_10_CELLS[j][i] for i in range(len(FULL_5HTr2b_10_CELLS[j])) if FULL_5HTr2b_10_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_5HTr2b_10_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_10_CELLS[j])) if FULL_5HTr2b_10_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
ax5.scatter(CELLS, CELLS_tDt)
ax6.hist(CELLS, alpha=0.7)
ax12.hist(CELLS_tDt, alpha=0.7, color='red')

CELLS = ([FULL_5HTr2b_18_CELLS[j][i] for i in range(len(FULL_5HTr2b_18_CELLS[j])) if FULL_5HTr2b_18_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_5HTr2b_18_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_18_CELLS[j])) if FULL_5HTr2b_18_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
ax5.scatter(CELLS, CELLS_tDt)
ax6.hist(CELLS, alpha=0.7)
ax12.hist(CELLS_tDt, alpha=0.7, color='red')

CELLS = ([FULL_5HTr2b_30_CELLS[j][i] for i in range(len(FULL_5HTr2b_30_CELLS[j])) if FULL_5HTr2b_30_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_5HTr2b_30_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_30_CELLS[j])) if FULL_5HTr2b_30_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
ax5.scatter(CELLS, CELLS_tDt)
ax6.hist(CELLS, alpha=0.7)
ax12.hist(CELLS_tDt, alpha=0.7, color='red')

CELLS = ([FULL_5HTr2b_CELLS[j][i] for i in range(len(FULL_5HTr2b_CELLS[j])) if FULL_5HTr2b_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_5HTr2b_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_CELLS[j])) if FULL_5HTr2b_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
ax5.scatter(CELLS, CELLS_tDt)
ax6.hist(CELLS, alpha=0.7)
ax12.hist(CELLS_tDt, alpha=0.7, color='red')

ax1.set_title('eGFP')
ax2.set_title('tD-Tomato')
ax3.set_title('Equalized Sub.')
ax4.set_title('eGFP/tD-tomato')
ax5.set_xlabel('eGFP')
ax5.set_ylabel('tD-Tomato')


CELLS = ([FULL_5HTr2b_10_CELLS[j][i] for i in range(len(FULL_5HTr2b_10_CELLS[j])) if FULL_5HTr2b_10_AREA[j][i]<100])
MEAN_5HTR2b_10_ASTRO = CELLS

CELLS_1 = ([FULL_5HTr2b_18_CELLS[j][i] for i in range(len(FULL_5HTr2b_18_CELLS[j])) if FULL_5HTr2b_18_AREA[j][i]<100])
MEAN_5HTR2b_18_ASTRO = CELLS_1

CELLS_2 = ([FULL_5HTr2b_30_CELLS[j][i] for i in range(len(FULL_5HTr2b_30_CELLS[j])) if FULL_5HTr2b_30_AREA[j][i]<100])
MEAN_5HTR2b_30_ASTRO = CELLS_2

CELLS_3 = ([FULL_5HTr2b_CELLS[j][i] for i in range(len(FULL_5HTr2b_CELLS[j])) if FULL_5HTr2b_AREA[j][i]>100])
MEAN_5HTR2b_37_ASTRO = CELLS_3
ax7.boxplot([CELLS, CELLS_1, CELLS_2, CELLS_3])


CELLS = ([FULL_5HTr2b_10_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_10_CELLS[j])) if FULL_5HTr2b_10_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_1 = ([FULL_5HTr2b_18_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_18_CELLS[j])) if FULL_5HTr2b_18_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_2 = ([FULL_5HTr2b_30_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_30_CELLS[j])) if FULL_5HTr2b_30_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_3 = ([FULL_5HTr2b_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_CELLS[j])) if FULL_5HTr2b_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
ax8.boxplot([CELLS, CELLS_1, CELLS_2, CELLS_3])


CELLS = ([FULL_5HTr2b_10_SUBTRACTED[j][i] for i in range(len(FULL_5HTr2b_10_CELLS[j])) if FULL_5HTr2b_10_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_1 = ([FULL_5HTr2b_18_SUBTRACTED[j][i] for i in range(len(FULL_5HTr2b_18_CELLS[j])) if FULL_5HTr2b_18_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_2 = ([FULL_5HTr2b_30_SUBTRACTED[j][i] for i in range(len(FULL_5HTr2b_30_CELLS[j])) if FULL_5HTr2b_30_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_3 = ([FULL_5HTr2b_SUBTRACTED[j][i] for i in range(len(FULL_5HTr2b_CELLS[j])) if FULL_5HTr2b_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)
ax9.boxplot([CELLS, CELLS_1, CELLS_2, CELLS_3])



CELLS = ([FULL_5HTr2b_10_DIVIDED[j][i] for i in range(len(FULL_5HTr2b_10_CELLS[j])) if FULL_5HTr2b_10_AREA[j][i]<100])
MEAN_5HTR2b_10_DIVIDED_ASTRO  = CELLS

CELLS_1 = ([FULL_5HTr2b_18_DIVIDED[j][i] for i in range(len(FULL_5HTr2b_18_CELLS[j])) if FULL_5HTr2b_18_AREA[j][i]<100])
MEAN_5HTR2b_18_DIVIDED_ASTRO  = CELLS_1

CELLS_2 = ([FULL_5HTr2b_30_DIVIDED[j][i] for i in range(len(FULL_5HTr2b_30_CELLS[j])) if FULL_5HTr2b_30_AREA[j][i]<100])
MEAN_5HTR2b_30_DIVIDED_ASTRO = CELLS_2

CELLS_3 = ([FULL_5HTr2b_DIVIDED[j][i] for i in range(len(FULL_5HTr2b_CELLS[j])) if FULL_5HTr2b_AREA[j][i]<100])
MEAN_5HTR2b_37_DIVIDED_ASTRO = CELLS_3
ax10.boxplot([CELLS, CELLS_1, CELLS_2, CELLS_3])



CELLS = ([FULL_5HTr2b_10_CELLS[j][i] for i in range(len(FULL_5HTr2b_10_CELLS[j])) if FULL_5HTr2b_10_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_5HTr2b_10_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_10_CELLS[j])) if FULL_5HTr2b_10_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)
ax11.scatter(CELLS, CELLS_tDt)

CELLS = ([FULL_5HTr2b_18_CELLS[j][i] for i in range(len(FULL_5HTr2b_18_CELLS[j])) if FULL_5HTr2b_18_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_5HTr2b_18_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_18_CELLS[j])) if FULL_5HTr2b_18_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)
ax11.scatter(CELLS, CELLS_tDt)

CELLS = ([FULL_5HTr2b_30_CELLS[j][i] for i in range(len(FULL_5HTr2b_30_CELLS[j])) if FULL_5HTr2b_30_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_5HTr2b_30_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_30_CELLS[j])) if FULL_5HTr2b_30_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)
ax11.scatter(CELLS, CELLS_tDt)

CELLS = ([FULL_5HTr2b_CELLS[j][i] for i in range(len(FULL_5HTr2b_CELLS[j])) if FULL_5HTr2b_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_5HTr2b_CELLS_tDt[j][i] for i in range(len(FULL_5HTr2b_CELLS[j])) if FULL_5HTr2b_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)
ax11.scatter(CELLS, CELLS_tDt)

ax7.set_title('eGFP Astro.')
ax8.set_title('tD-Tomato Astro.')
ax9.set_title('Equalized Sub. Astro.')
ax10.set_title('eGFP/tD-tomato Astro.')
ax11.set_xlabel('eGFP')
ax11.set_ylabel('tD-Tomato')

plt.tight_layout()





 #Plots for 5HTr2b group

plt.figure(figsize=(10, 4), num='Others Group')
ax1 = plt.subplot(2, 6, 1)
ax2 = plt.subplot(2, 6, 2)
ax3 = plt.subplot(2, 6, 3)
ax4 = plt.subplot(2, 6, 4)
ax5 = plt.subplot(2, 6, 5)
ax6 = plt.subplot(2, 6, 6)
ax7 = plt.subplot(2, 6, 7)
ax8 = plt.subplot(2, 6, 8)
ax9 = plt.subplot(2, 6, 9)
ax10 = plt.subplot(2, 6, 10)
ax11 = plt.subplot(2, 6, 11)
ax12 = plt.subplot(2, 6, 12)

j=0

CELLS = ([FULL_SUSD4_CELLS[j][i] for i in range(len(FULL_SUSD4_CELLS[j])) if FULL_SUSD4_AREA[j][i]>100])
MEAN_SUSD4_NEURO = CELLS

CELLS_1 = ([FULL_PDX1_CELLS[j][i] for i in range(len(FULL_PDX1_CELLS[j])) if FULL_PDX1_AREA[j][i]>100])
MEAN_PDX1_NEURO  = CELLS_1


ax1.boxplot([CELLS, CELLS_1])


CELLS = ([FULL_SUSD4_CELLS_tDt[j][i] for i in range(len(FULL_SUSD4_CELLS[j])) if FULL_SUSD4_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
MEAN_SUSD4_NEURO_tDt = CELLS

CELLS_1 = ([FULL_PDX1_CELLS_tDt[j][i] for i in range(len(FULL_PDX1_CELLS[j])) if FULL_PDX1_AREA[j][i]>100])
MEAN = np.nanmean(CELLS_1)
MEAN_PDX1_NEURO_tDt = CELLS_1

ax2.boxplot([CELLS, CELLS_1])


CELLS = ([FULL_SUSD4_SUBTRACTED[j][i] for i in range(len(FULL_SUSD4_CELLS[j])) if FULL_SUSD4_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
MEAN_SUSD4_SUBTRACTED_NEURO = CELLS

CELLS_1 = ([FULL_PDX1_SUBTRACTED[j][i] for i in range(len(FULL_PDX1_CELLS[j])) if FULL_PDX1_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
MEAN_PDX1_SUBTRACTED_NEURO = CELLS_1

ax3.boxplot([CELLS, CELLS_1])



CELLS = ([FULL_SUSD4_DIVIDED[j][i] for i in range(len(FULL_SUSD4_CELLS[j])) if FULL_SUSD4_AREA[j][i]>100])
MEAN_SUSD4_DIVIDED_NEURO = CELLS

CELLS_1 = ([FULL_PDX1_DIVIDED[j][i] for i in range(len(FULL_PDX1_CELLS[j])) if FULL_PDX1_AREA[j][i]>100])
MEAN_PDX1_DIVIDED_NEURO = CELLS_1
ax4.boxplot([CELLS, CELLS_1])



CELLS = ([FULL_SUSD4_CELLS[j][i] for i in range(len(FULL_SUSD4_CELLS[j])) if FULL_SUSD4_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_SUSD4_CELLS_tDt[j][i] for i in range(len(FULL_SUSD4_CELLS[j])) if FULL_SUSD4_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
ax5.scatter(CELLS, CELLS_tDt)
ax6.hist(CELLS, alpha=0.7)
ax12.hist(CELLS_tDt, alpha=0.7, color='red')

CELLS = ([FULL_PDX1_CELLS[j][i] for i in range(len(FULL_PDX1_CELLS[j])) if FULL_PDX1_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_PDX1_CELLS_tDt[j][i] for i in range(len(FULL_PDX1_CELLS[j])) if FULL_PDX1_AREA[j][i]>100])
MEAN = np.nanmean(CELLS)
ax5.scatter(CELLS, CELLS_tDt)
ax6.hist(CELLS, alpha=0.7)
ax12.hist(CELLS_tDt, alpha=0.7, color='red')

ax1.set_title('eGFP')
ax2.set_title('tD-Tomato')
ax3.set_title('Equalized Sub.')
ax4.set_title('eGFP/tD-tomato')
ax5.set_xlabel('eGFP')
ax5.set_ylabel('tD-Tomato')


CELLS = ([FULL_SUSD4_CELLS[j][i] for i in range(len(FULL_SUSD4_CELLS[j])) if FULL_SUSD4_AREA[j][i]<100])
MEAN_SUSD4_ASTRO = CELLS

CELLS_1 = ([FULL_PDX1_CELLS[j][i] for i in range(len(FULL_PDX1_CELLS[j])) if FULL_PDX1_AREA[j][i]<100])
MEAN_PDX1_ASTRO = CELLS


ax7.boxplot([CELLS, CELLS_1])


CELLS = ([FULL_SUSD4_CELLS_tDt[j][i] for i in range(len(FULL_SUSD4_CELLS[j])) if FULL_SUSD4_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_1 = ([FULL_PDX1_CELLS_tDt[j][i] for i in range(len(FULL_PDX1_CELLS[j])) if FULL_PDX1_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

ax8.boxplot([CELLS, CELLS_1])


CELLS = ([FULL_SUSD4_SUBTRACTED[j][i] for i in range(len(FULL_SUSD4_CELLS[j])) if FULL_SUSD4_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_1 = ([FULL_PDX1_SUBTRACTED[j][i] for i in range(len(FULL_PDX1_CELLS[j])) if FULL_PDX1_AREA[j][i]<100])
MEAN = np.nanmean(CELLS_1)

ax9.boxplot([CELLS, CELLS_1])



CELLS = ([FULL_SUSD4_DIVIDED[j][i] for i in range(len(FULL_SUSD4_CELLS[j])) if FULL_SUSD4_AREA[j][i]<100])
MEAN_SUSD4_DIVIDED_ASTRO = CELLS

CELLS_1 = ([FULL_PDX1_DIVIDED[j][i] for i in range(len(FULL_PDX1_CELLS[j])) if FULL_PDX1_AREA[j][i]<100])
MEAN_PDX1_DIVIDED_ASTRO = CELLS_1
ax9.boxplot([CELLS, CELLS_1])



CELLS = ([FULL_SUSD4_CELLS[j][i] for i in range(len(FULL_SUSD4_CELLS[j])) if FULL_SUSD4_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_SUSD4_CELLS_tDt[j][i] for i in range(len(FULL_SUSD4_CELLS[j])) if FULL_SUSD4_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)
ax11.scatter(CELLS, CELLS_tDt)

CELLS = ([FULL_PDX1_CELLS[j][i] for i in range(len(FULL_PDX1_CELLS[j])) if FULL_PDX1_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)

CELLS_tDt = ([FULL_PDX1_CELLS_tDt[j][i] for i in range(len(FULL_PDX1_CELLS[j])) if FULL_PDX1_AREA[j][i]<100])
MEAN = np.nanmean(CELLS)
ax11.scatter(CELLS, CELLS_tDt)

ax7.set_title('eGFP Astro.')
ax8.set_title('tD-Tomato Astro.')
ax9.set_title('Equalized Sub. Astro.')
ax10.set_title('eGFP/tD-tomato Astro.')
ax11.set_xlabel('eGFP')
ax11.set_ylabel('tD-Tomato')

plt.tight_layout()