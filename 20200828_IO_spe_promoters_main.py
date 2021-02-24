# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 15:38:21 2020

@author: dorga
"""
def EqualizeHistogram(a, bins):
	a = np.array(a)
	hist, bins2 = np.histogram(a, bins=bins)

	#Compute CDF from histogram
	cdf = np.cumsum(hist, dtype=np.float64)
	cdf = np.hstack(([0], cdf))
	cdf = cdf / cdf[-1]

	#Do equalization
	binnum = np.digitize(a, bins, True)-1
	neg = np.where(binnum < 0)
	binnum[neg] = 0

	aeq = cdf[binnum] * bins[-1]

	return aeq

#IMAGE ANALYSIS PLOT
import matplotlib.backends.backend_pdf
import seaborn as sns
from PIL import Image

pdf = matplotlib.backends.backend_pdf.PdfPages(r'C:\Users\dorga\Desktop\output.pdf')
DIR__, DIRECTORY = load_directory_content__()

HALF_LABEL_LIST = pd.read_csv(DIRECTORY+r'\20200828_HALF_IO_DILUTION_LABELS.csv')

FULL_RAW_HIST_eGFP = []
FULL_RAW_HIST_tDt = []

FULL_tDt_IO_OUT_RATIO = []
FULL_eGFP_IO_OUT_RATIO = []
FULL_tDt_IO_IN = []
FULL_eGFP_IO_IN = []
FULL_tDt_IO_OUT = []
FULL_eGFP_IO_OUT = []
FULL_tDt_IO_IN_UP = []
FULL_eGFP_IO_IN_UP = []
FULL_tDt_IO_IN_DOWN = []
FULL_eGFP_IO_IN_DOWN = []

FULL_eGFP_SPREAD_IN = []
FULL_eGFP_SPREAD_OUT = []
FULL_tDt_SPREAD_IN = []
FULL_tDt_SPREAD_OUT = []

PixNumeGFP = []
PixNumTdT = []

FULL_PROMOTER_NAMES = []
eGFP_HALF_IO = []
eGFP_HALF_IO_OUT = []
tDt_HALF_IO = []
tDt_HALF_IO_OUT = []
eGFP_HALF_IO_DILUTION = []
eGFP_HALF_IO_LABEL = []

for PATH in DIR__:   
    if ('.tif' in PATH)==True and ('KD' in PATH)==True:
        if True: 
            #PATH = load_specific_file__()
            
            w = 1024
            IM = Image.open(PATH)
            SPLIT_IMAGE = Image.Image.split(IM.convert('RGB'))
            
            fig = plt.figure(num=PATH, figsize=(13,3))
            ax= plt.subplot(151)
            ax2= plt.subplot(152)
            ax3= plt.subplot(153)
            ax4= plt.subplot(154)
            ax5= plt.subplot(155)
            
            ax.imshow(np.array(SPLIT_IMAGE[0]), cmap=cm.Reds)
            ax2.imshow(np.array(SPLIT_IMAGE[1]), cmap=cm.Greens)
            ax3.imshow(np.array(SPLIT_IMAGE[2]), cmap=cm.Blues)
            
            tdTomato = np.concatenate(np.array(SPLIT_IMAGE[0]))
            eGFP = np.concatenate(np.array(SPLIT_IMAGE[1]))
            
            """
            I = EqualizeHistogram(eGFP, np.linspace(0, 255, 257))
            I_equalized = []
            for k in range(len(I)):
                if I[k]>240:
                    I_equalized.append(I[k])
                else:
                    I_equalized.append(np.nan)
            eGFP = I_equalized
            
            I = EqualizeHistogram(tdTomato, np.linspace(0, 255, 257))
            I_equalized = []
            for k in range(len(I)):
                if I[k]>240:
                    I_equalized.append(I[k])
                else:
                    I_equalized.append(np.nan)
            tdTomato = I_equalized
            """
            
            FULL_RAW_HIST_eGFP.append(eGFP)
            FULL_RAW_HIST_tDt.append(tdTomato)
            
            SUB = np.subtract(np.nan_to_num(eGFP), np.nan_to_num(tdTomato))
            
            ax4.imshow(np.reshape(SUB, (-1, 1024)), cmap=cm.magma,interpolation='hamming')
            
            
            
            IM_ARRAY = np.array(np.reshape(eGFP, (-1, w)))
            Z_SCORE = np.nanmean(IM_ARRAY) + 3*np.nanstd(IM_ARRAY)
            ROI_crop = pd.read_csv(PATH.split('tif')[0]+'csv')
            eGFP_IO = []
            eGFP_IO_sup = []
            eGFP_IO_inf = []
            eGFP_OUT_sup = []
            eGFP_OUT_inf = []
            
            for k in range(len(ROI_crop.X)):
                eGFP_IO.append(IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]])
                if ROI_crop.Y[k]<w*0.5:
                    eGFP_IO_sup.append(IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]])
                elif ROI_crop.Y[k]>w*0.5:
                    eGFP_IO_inf.append(IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]])
                IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]] = 0
            
            eGFP_ARRAY_CONC = np.concatenate(IM_ARRAY)
            #Z_SCORE = np.nanmedian(np.concatenate(np.transpose(IM_ARRAY)[0:int(w/4)]))+3*np.std(np.concatenate(np.transpose(IM_ARRAY)[0:int(w/2)]))
            
            eGFP_ARRAY_CONC = [num for num in eGFP_ARRAY_CONC if num>Z_SCORE]
            eGFP_IO = [num for num in eGFP_IO if num>Z_SCORE]
            
            eGFP_OUT_sup = [eGFP_ARRAY_CONC[k] for k in range(len(eGFP_ARRAY_CONC)) if eGFP_ARRAY_CONC[k]>Z_SCORE and  k<int(len(eGFP_ARRAY_CONC)*0.4)]
            eGFP_OUT_inf = [eGFP_ARRAY_CONC[k] for k in range(len(eGFP_ARRAY_CONC)) if eGFP_ARRAY_CONC[k]>Z_SCORE and  k>int(len(eGFP_ARRAY_CONC)*0.6)]
            
            PixNumeGFP.append([len(eGFP_OUT_sup), len(eGFP_OUT_inf)])
            
            IM_ARRAY = np.array(np.reshape(tdTomato, (-1, w)))
            Z_SCORE = np.nanmean(IM_ARRAY) + 3*np.nanstd(IM_ARRAY)
            ROI_crop = pd.read_csv(PATH.split('tif')[0]+'csv')
            tDTomato_IO = []
            tDTomato_IO_sup = []
            tDTomato_IO_inf = []
            tDTomato_OUT_sup = []
            tDTomato_OUT_inf = []
            for k in range(len(ROI_crop.X)):
                tDTomato_IO.append(IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]])
                if ROI_crop.Y[k]<w*0.5:
                    tDTomato_IO_sup.append(IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]])
                elif ROI_crop.Y[k]>w*0.5:
                    tDTomato_IO_inf.append(IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]])
                IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]] = 0
            
            tDt_CONC = np.concatenate(IM_ARRAY)
            #Z_SCORE = np.nanmedian(np.concatenate(np.transpose(IM_ARRAY)[0:int(w/4)]))+3*np.std(np.concatenate(np.transpose(IM_ARRAY)[0:int(w/2)]))
            
            tDt_CONC = [num for num in tDt_CONC if num>Z_SCORE]
            tDTomato_IO = [num for num in tDTomato_IO if num>Z_SCORE] 
            
            tDTomato_OUT_sup = [tDt_CONC[k] for k in range(len(tDt_CONC)) if tDt_CONC[k]>Z_SCORE and  k<int(len(tDt_CONC)*0.4)]
            tDTomato_OUT_inf = [tDt_CONC[k] for k in range(len(tDt_CONC)) if tDt_CONC[k]>Z_SCORE and  k>int(len(tDt_CONC)*0.6)]
            
            PixNumTdT.append([len(tDTomato_OUT_sup), len(tDTomato_OUT_inf)])
            
            SUB = np.multiply(SUB, 1)
            ax5.boxplot([tDt_CONC, eGFP_ARRAY_CONC, tDTomato_IO, eGFP_IO], showfliers=False)
            
            promoter_name_ = PATH.split('\\')[-1].split('-')[0].split('_')[-1]
            
            filename = PATH.split('\\')[-1].split('.tif')[0]
            
            for k in range(len(HALF_LABEL_LIST.NAME.tolist())):
                if HALF_LABEL_LIST.NAME[k] == filename:
                    position__ = HALF_LABEL_LIST.POSITION[k]
                    dilution__ = HALF_LABEL_LIST.DILUTION[k]
                    if position__=='UP':
                        eGFP_HALF_IO.append(np.nanmean(eGFP_IO_sup))
                        eGFP_HALF_IO_OUT.append(np.nanmean(eGFP_OUT_sup))
                        tDt_HALF_IO.append(np.nanmean(tDTomato_IO_sup))
                        tDt_HALF_IO_OUT.append(np.nanmean(tDTomato_OUT_sup))
                    elif position__=='DOWN':
                        eGFP_HALF_IO.append(np.nanmean(eGFP_IO_inf))
                        eGFP_HALF_IO_OUT.append(np.nanmean(eGFP_OUT_inf))
                        tDt_HALF_IO.append(np.nanmean(tDTomato_IO_inf))
                        tDt_HALF_IO_OUT.append(np.nanmean(tDTomato_OUT_inf))
                    eGFP_HALF_IO_LABEL.append(promoter_name_)
                    eGFP_HALF_IO_DILUTION.append(dilution__)
            
            FULL_PROMOTER_NAMES.append(promoter_name_)
            FULL_eGFP_IO_OUT_RATIO.append(np.nanmean(eGFP_IO) / np.nanmean(eGFP_ARRAY_CONC))
            FULL_tDt_IO_OUT_RATIO.append(np.nanmean(tDTomato_IO) / np.nanmean(tDt_CONC))
            
            FULL_tDt_IO_IN.append(np.nanmean(tDTomato_IO))
            FULL_eGFP_IO_IN.append(np.nanmean(eGFP_IO))
            FULL_tDt_IO_OUT.append(np.nanmean(tDt_CONC))
            FULL_eGFP_IO_OUT.append(np.nanmean(eGFP_ARRAY_CONC))
            
            
            FULL_eGFP_SPREAD_IN.append(len(eGFP_IO))
            FULL_eGFP_SPREAD_OUT.append(len(eGFP_ARRAY_CONC))
            FULL_tDt_SPREAD_IN.append(len(tDTomato_IO))
            FULL_tDt_SPREAD_OUT.append(len(tDt_CONC))

            
            FULL_tDt_IO_IN_UP.append(np.nanmean(tDTomato_IO_sup))
            FULL_eGFP_IO_IN_UP.append(np.nanmean(eGFP_IO_sup))
            FULL_tDt_IO_IN_DOWN.append(np.nanmean(tDTomato_IO_inf))
            FULL_eGFP_IO_IN_DOWN.append(np.nanmean(eGFP_IO_inf))
            
            ax.set_title('CAG-tdTomato')
            ax2.set_title('p-eGFP')
            ax3.set_title('DAPI')
            ax4.set_title('Difference image')
            ax5.set_title('Difference image')
            ax.set_axis_off()
            ax2.set_axis_off()
            ax3.set_axis_off()
            ax4.set_axis_off()
            plt.title(PATH)
            plt.tight_layout()
            pdf.savefig(fig)
        #except:
            #pass

pdf.close()


'''
#EXPRESSION PLOT
plt.figure(figsize=(15,2))
ax = plt.subplot(1, 10, 1)
ax2 = plt.subplot(1, 10, 2, sharex=ax, sharey=ax)
ax3 = plt.subplot(1, 10, 3, sharex=ax, sharey=ax)
ax4 = plt.subplot(1, 10, 4, sharex=ax, sharey=ax)
ax5 = plt.subplot(1, 10, 5, sharex=ax, sharey=ax)
ax6 = plt.subplot(1, 10, 6)
ax7 = plt.subplot(1, 10, 7)
ax8 = plt.subplot(1, 10, 8)
ax9 = plt.subplot(1, 10, 9)
ax10 = plt.subplot(1, 10, 10)

for i in range(len(FULL_PROMOTER_NAMES)):
    if ('2.5' in FULL_PROMOTER_NAMES[i])==True:
        ax.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='blue')
        #ax5.scatter('Igsf', FULL_eGFP_IO_OUT_RATIO[i], color='blue')
    elif ('5HTr' in FULL_PROMOTER_NAMES[i])==True and ('1.0' in FULL_PROMOTER_NAMES[i])==False and ('1.8' in FULL_PROMOTER_NAMES[i])==False:
        ax3.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='red')
        #ax5.scatter('5HTr', FULL_eGFP_IO_OUT_RATIO[i], color='red')
    elif ('PDX' in FULL_PROMOTER_NAMES[i])==True:
        ax4.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='purple')
        #ax5.scatter('PDX1', FULL_eGFP_IO_OUT_RATIO[i], color='purple')
    elif ('SUS' in FULL_PROMOTER_NAMES[i])==True:
        ax5.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='orange')
        #ax5.scatter('SUSD4', FULL_eGFP_IO_OUT_RATIO[i], color='orange')
    elif ('1.3' in FULL_PROMOTER_NAMES[i])==True:
        ax2.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='orange')
        #ax5.scatter('SUSD4', FULL_eGFP_IO_OUT_RATIO[i], color='orange')



ax.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)
ax2.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)
ax3.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)
ax4.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)
ax5.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)

a = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='3.7kb']
d = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']
ax6.boxplot([a, b, c, d, e, f, g, h], labels=['Igsf9(2.5)', 'Igsf9(1.3)', 'Igsf9(3.7)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4'])


a = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='3.7kb']
d = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_eGFP_IO_IN[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']
ax7.boxplot([a, b, c, d, e, f, g, h], labels=['Igsf9(2.5)', 'Igsf9(1.3)', 'Igsf9(3.7)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4'])


a = [FULL_eGFP_IO_OUT[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_eGFP_IO_OUT[k]for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_eGFP_IO_OUT[k]for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
d = [FULL_eGFP_IO_OUT[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_eGFP_IO_OUT[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_eGFP_IO_OUT[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_eGFP_IO_OUT[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_eGFP_IO_OUT[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']
ax8.boxplot([a, b, c, d, e, f, g, h], labels=['Igsf9(2.5)', 'Igsf9(1.3)', 'Igsf9(3.7)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4'])



a = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
d = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']
ax9.boxplot([a, b, c, d, e, f, g, h], labels=['Igsf9(2.5)', 'Igsf9(1.3)', 'Igsf9(3.7)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4'])


a = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='3.7kb']
d = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']
ax10.boxplot([a, b, c, d, e, f, g, h], labels=['Igsf9(2.5)', 'Igsf9(1.3)', 'Igsf9(3.7)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4'])

ax.set_ylabel('eGFP IO / eGFP OUT')
ax2.set_ylabel('eGFP IO / eGFP OUT')
ax3.set_ylabel('eGFP IO / eGFP OUT')
ax4.set_ylabel('eGFP IO / eGFP OUT')
ax5.set_ylabel('eGFP IO / eGFP OUT')

ax.set_xlabel('tDt-IO / tDt OUT')
ax2.set_xlabel('tDt-IO / tDt OUT')
ax3.set_xlabel('tDt-IO / tDt OUT')
ax4.set_xlabel('tDt-IO / tDt OUT')

ax.set_title('Igsf9(2.5kb)-eGFP')
ax2.set_title('Igsf9(1.3kb)-eGFP')
ax3.set_title('5HTR2b-eGFP')
ax4.set_title('PDX1-eGFP')
ax5.set_title('SUSD4-eGFP')
ax6.set_title('eGFP(IO)/tDt(IO)')
ax7.set_title('eGFP(IO)')
ax8.set_title('eGFP(OUT)')
plt.tight_layout()
'''

#CONCENTRATION PLOT
plt.figure()
ax = plt.subplot(5, 2, 1)
ax2 = plt.subplot(5, 2, 2)
ax3 = plt.subplot(5, 2, 3)
ax4 = plt.subplot(5, 2, 4)
ax5 = plt.subplot(5, 2, 5)
ax6 = plt.subplot(5, 2, 6)
ax7 = plt.subplot(5, 2, 7)
ax8 = plt.subplot(5, 2, 8)
ax9 = plt.subplot(5, 2, 9)
ax10 = plt.subplot(5, 2, 10)

promoter = '2.5kb'
a = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax.set_title(promoter)

a = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax2.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax.set_title('Igsf9(2.5kb)-eGFP')
ax2.set_title('CAG-tDt')
ax.set_xlabel('Dilution')
ax2.set_xlabel('Dilution')
ax.set_ylabel('IO-spe./non-spe.')

promoter = '1.3kb'
a = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax3.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax3.set_title(promoter)

a = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax4.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax3.set_title('Igsf9(1.3kb)-eGFP')
ax4.set_title('CAG-tDt')
ax3.set_xlabel('Dilution')
ax4.set_xlabel('Dilution')
ax3.set_ylabel('IO-spe./non-spe.')

promoter = '5HT'
a = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax5.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax5.set_title(promoter)

a = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax6.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax5.set_title(promoter+'-eGFP')
ax6.set_title('CAG-tDt')
ax5.set_xlabel('Dilution')
ax6.set_xlabel('Dilution')
ax5.set_ylabel('IO-spe./non-spe.')

promoter = 'PDX'
a = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax7.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax7.set_title(promoter)

a = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax8.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax7.set_title(promoter+'-eGFP')
ax8.set_title('CAG-tDt')
ax7.set_xlabel('Dilution')
ax8.set_xlabel('Dilution')
ax7.set_ylabel('IO-spe./non-spe.')

promoter = 'SUSD4'
a = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax9.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax9.set_title(promoter)

a = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax10.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax9.set_title(promoter+'-eGFP')
ax10.set_title('CAG-tDt')
ax9.set_xlabel('Dilution')
ax10.set_xlabel('Dilution')
ax9.set_ylabel('IO-spe./non-spe.')



#EXPRESSION PLOT
plt.figure(figsize=(15,2))
ax = plt.subplot(1, 12, 1)
ax2 = plt.subplot(1, 12, 2, sharex=ax, sharey=ax)
ax3 = plt.subplot(1, 12, 3, sharex=ax, sharey=ax)
ax4 = plt.subplot(1, 12, 4, sharex=ax, sharey=ax)
ax5 = plt.subplot(1, 12, 5, sharex=ax, sharey=ax)
ax6 = plt.subplot(1, 12, 6)
ax7 = plt.subplot(1, 12, 7)
ax8 = plt.subplot(1, 12, 8)
ax9 = plt.subplot(1, 12, 9)
ax10 = plt.subplot(1, 12, 10)
ax11 = plt.subplot(1, 12, 11)
ax12 = plt.subplot(1, 12, 12)

for i in range(len(FULL_PROMOTER_NAMES)):
    if ('2.5' in FULL_PROMOTER_NAMES[i])==True:
        ax.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='blue')
        #ax5.scatter('Igsf', FULL_eGFP_IO_OUT_RATIO[i], color='blue')
    elif ('5HTr' in FULL_PROMOTER_NAMES[i])==True:
        ax3.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='red')
        #ax5.scatter('5HTr', FULL_eGFP_IO_OUT_RATIO[i], color='red')
    elif ('PDX' in FULL_PROMOTER_NAMES[i])==True:
        ax4.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='purple')
        #ax5.scatter('PDX1', FULL_eGFP_IO_OUT_RATIO[i], color='purple')
    elif ('SUS' in FULL_PROMOTER_NAMES[i])==True:
        ax5.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='orange')
        #ax5.scatter('SUSD4', FULL_eGFP_IO_OUT_RATIO[i], color='orange')
    elif ('1.3' in FULL_PROMOTER_NAMES[i])==True:
        ax2.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='orange')
        #ax5.scatter('SUSD4', FULL_eGFP_IO_OUT_RATIO[i], color='orange')



ax.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)
ax2.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)
ax3.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)
ax4.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)
ax5.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)

a = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='Igsf9(3.7kb)']
d = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']

FULL_eGFP_IO_IN_ALL = [np.nanmean(a), np.nanmean(b), np.nanmean(c), np.nanmean(d), np.nanmean(e), np.nanmean(f), np.nanmean(g), np.nanmean(h)]
ax6.boxplot([a, b, c, d, e, f, g, h], labels=['Igsf9(2.5)', 'Igsf9(1.3)', 'Igsf9(3.7)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4'])


a = [FULL_tDt_IO_IN[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_tDt_IO_IN[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_tDt_IO_IN[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='Igsf9(3.7kb)']
d = [FULL_tDt_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_tDt_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_tDt_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_tDt_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_tDt_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']
FULL_tDt_IO_IN_ALL = [np.nanmean(a), np.nanmean(b), np.nanmean(c), np.nanmean(d), np.nanmean(e), np.nanmean(f), np.nanmean(g), np.nanmean(h)]
ax7.boxplot([a, b, c, d, e, f, g, h], labels=['Igsf9(2.5)', 'Igsf9(1.3)', 'Igsf9(3.7)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4'])


a = [FULL_eGFP_IO_OUT_RATIO[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_eGFP_IO_OUT_RATIO[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_eGFP_IO_OUT_RATIO[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='Igsf9(3.7kb)']
d = [FULL_eGFP_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_eGFP_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_eGFP_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_eGFP_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_eGFP_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']
FULL_eGFP_IO_OUT_RATIO_ALL = [np.nanmean(a), np.nanmean(b), np.nanmean(c), np.nanmean(d), np.nanmean(e), np.nanmean(f), np.nanmean(g), np.nanmean(h)]
ax8.boxplot([a, b, c, d, e, f, g, h], labels=['Igsf9(2.5)', 'Igsf9(1.3)', 'Igsf9(3.7)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4'])

a = [FULL_tDt_IO_OUT[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_tDt_IO_OUT[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_tDt_IO_OUT[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='Igsf9(3.7kb)']
d = [FULL_tDt_IO_OUT[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_tDt_IO_OUT[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_tDt_IO_OUT[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_tDt_IO_OUT[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_tDt_IO_OUT[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']
FULL_tDt_IO_OUT_ALL = [np.nanmean(a), np.nanmean(b), np.nanmean(c), np.nanmean(d), np.nanmean(e), np.nanmean(f), np.nanmean(g), np.nanmean(h)]

a = [FULL_eGFP_IO_OUT[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_eGFP_IO_OUT[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_eGFP_IO_OUT[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='Igsf9(3.7kb)']
d = [FULL_eGFP_IO_OUT[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_eGFP_IO_OUT[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_eGFP_IO_OUT[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_eGFP_IO_OUT[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_eGFP_IO_OUT[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']
FULL_eGFP_IO_OUT_ALL = [np.nanmean(a), np.nanmean(b), np.nanmean(c), np.nanmean(d), np.nanmean(e), np.nanmean(f), np.nanmean(g), np.nanmean(h)]

a = [FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='Igsf9(3.7kb)']
d = [FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']
FULL_tDt_IO_OUT_RATIO_ALL = [np.nanmean(a), np.nanmean(b), np.nanmean(c), np.nanmean(d), np.nanmean(e), np.nanmean(f), np.nanmean(g), np.nanmean(h)]
ax9.boxplot([a, b, c, d, e, f, g, h], labels=['Igsf9(2.5)', 'Igsf9(1.3)', 'Igsf9(3.7)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4'])


a = [FULL_eGFP_IO_OUT_RATIO[k]/FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_eGFP_IO_OUT_RATIO[k]/FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_eGFP_IO_OUT_RATIO[k]/FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='Igsf9(3.7kb)']
d = [FULL_eGFP_IO_OUT_RATIO[k]/FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_eGFP_IO_OUT_RATIO[k]/FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_eGFP_IO_OUT_RATIO[k]/FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_eGFP_IO_OUT_RATIO[k]/FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_eGFP_IO_OUT_RATIO[k]/FULL_tDt_IO_OUT_RATIO[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']
FULL_eGFP_IO_OUT_RATIO_div_tdT_OUT_RATIO_ALL = [np.nanmean(a), np.nanmean(b), np.nanmean(c), np.nanmean(d), np.nanmean(e), np.nanmean(f), np.nanmean(g), np.nanmean(h)]
ax10.boxplot([a, b, c, d, e, f, g, h], labels=['Igsf9(2.5)', 'Igsf9(1.3)', 'Igsf9(3.7)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4'])

a = [FULL_eGFP_IO_IN[k]/ FULL_tDt_IO_IN[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_eGFP_IO_IN[k]/ FULL_tDt_IO_IN[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_eGFP_IO_IN[k]/ FULL_tDt_IO_IN[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='Igsf9(3.7kb)']
d = [FULL_eGFP_IO_IN[k]/ FULL_tDt_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_eGFP_IO_IN[k]/ FULL_tDt_IO_IN[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_eGFP_IO_IN[k]/ FULL_tDt_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_eGFP_IO_IN[k]/ FULL_tDt_IO_IN[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_eGFP_IO_IN[k]/ FULL_tDt_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']
FULL_eGFP_IO_IN_RATIO_div_tdT_IN_RATIO_ALL = [np.nanmean(a), np.nanmean(b), np.nanmean(c), np.nanmean(d), np.nanmean(e), np.nanmean(f), np.nanmean(g), np.nanmean(h)]
ax11.boxplot([a, b, c, d, e, f, g, h], labels=['Igsf9(2.5)', 'Igsf9(1.3)', 'Igsf9(3.7)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4'])

a = [FULL_eGFP_IO_OUT[k]/ FULL_tDt_IO_OUT[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='2.5kb']
b = [FULL_eGFP_IO_OUT[k]/ FULL_tDt_IO_OUT[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='1.3kb']
c = [FULL_eGFP_IO_OUT[k]/ FULL_tDt_IO_OUT[k]  for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='Igsf9(3.7kb)']
d = [FULL_eGFP_IO_OUT[k]/ FULL_tDt_IO_OUT[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b']
e = [FULL_eGFP_IO_OUT[k]/ FULL_tDt_IO_OUT[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.8kb)']
f = [FULL_eGFP_IO_OUT[k]/ FULL_tDt_IO_OUT[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b(1.0kb)']
g = [FULL_eGFP_IO_OUT[k]/ FULL_tDt_IO_OUT[k]  for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1']
h = [FULL_eGFP_IO_OUT[k]/ FULL_tDt_IO_OUT[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4']
FULL_eGFP_IO_OUT_div_tdT_OUT_ALL = [np.nanmean(a), np.nanmean(b), np.nanmean(c), np.nanmean(d), np.nanmean(e), np.nanmean(f), np.nanmean(g), np.nanmean(h)]
ax12.boxplot([a, b, c, d, e, f, g, h], labels=['Igsf9(2.5)', 'Igsf9(1.3)', 'Igsf9(3.7)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4'])

ax.set_ylabel('eGFP IO / eGFP OUT')
ax2.set_ylabel('eGFP IO / eGFP OUT')
ax3.set_ylabel('eGFP IO / eGFP OUT')
ax4.set_ylabel('eGFP IO / eGFP OUT')
ax5.set_ylabel('eGFP IO / eGFP OUT')

ax.set_xlabel('tDt-IO / tDt OUT')
ax2.set_xlabel('tDt-IO / tDt OUT')
ax3.set_xlabel('tDt-IO / tDt OUT')
ax4.set_xlabel('tDt-IO / tDt OUT')

ax.set_title('Igsf9(2.5kb)-eGFP')
ax2.set_title('Igsf9(1.3kb)-eGFP')
ax3.set_title('5HTR2b-eGFP')
ax4.set_title('PDX1-eGFP')
ax5.set_title('SUSD4-eGFP')

ax6.set_title('eGFP(IO)')
ax7.set_title('tdTomato(IO)')
ax8.set_title('eGFP(IO)/eGFP(OUT)')
ax9.set_title('tdTomato(IO)/tdTomato(OUT)')
ax10.set_title('eGFP(RATIO)/tdTomato(RATIO)')
ax11.set_title('eGFP(IO)/tdTomato(IO)')
ax12.set_title('eGFP(OUT)/tdTomato(OUT)')
plt.tight_layout()


#AVERAGE RECAP MAP
plt.figure(figsize=(6,3))
s = ['Igsf9(2.5)', 'Igsf9(1.3)', 'Igsf9(3.7)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4']
ax1 = plt.subplot(141)
ax2 = plt.subplot(142)
ax3 = plt.subplot(143)
ax4 = plt.subplot(144)

ax1.scatter(FULL_eGFP_IO_IN_ALL, FULL_tDt_IO_IN_ALL)
for i in range(len(s)):
    ax1.text(FULL_eGFP_IO_IN_ALL[i], FULL_tDt_IO_IN_ALL[i], s[i]) 
ax1.set_xlabel('eGFP IO')
ax1.set_ylabel('tdT IO')

ax2.scatter(FULL_eGFP_IO_OUT_ALL, FULL_tDt_IO_OUT_ALL)
for i in range(len(s)):
    ax2.text(FULL_eGFP_IO_OUT_ALL[i], FULL_tDt_IO_OUT_ALL[i], s[i]) 
ax2.set_xlabel('eGFP OUT')
ax2.set_ylabel('tdT OUT')

ax3.scatter(FULL_eGFP_IO_OUT_RATIO_ALL, FULL_tDt_IO_OUT_RATIO_ALL)
for i in range(len(s)):
    ax3.text(FULL_eGFP_IO_OUT_RATIO_ALL[i], FULL_tDt_IO_OUT_RATIO_ALL[i], s[i]) 
ax3.set_xlabel('eGFP (IN/OUT)')
ax3.set_ylabel('tdT (IN/OUT)')

ax4.scatter(FULL_eGFP_IO_IN_RATIO_div_tdT_IN_RATIO_ALL, FULL_eGFP_IO_OUT_div_tdT_OUT_ALL)
for i in range(len(s)):
    ax4.text(FULL_eGFP_IO_IN_RATIO_div_tdT_IN_RATIO_ALL[i], FULL_eGFP_IO_OUT_div_tdT_OUT_ALL[i], s[i]) 
ax4.set_xlabel('Inside IO (eGFP/tDt)')
ax4.set_ylabel('Around IO (eGFP/tDt)')




#CONCENTRATION PLOT
plt.figure()
ax = plt.subplot(5, 2, 1)
ax2 = plt.subplot(5, 2, 2)
ax3 = plt.subplot(5, 2, 3)
ax4 = plt.subplot(5, 2, 4)
ax5 = plt.subplot(5, 2, 5)
ax6 = plt.subplot(5, 2, 6)
ax7 = plt.subplot(5, 2, 7)
ax8 = plt.subplot(5, 2, 8)
ax9 = plt.subplot(5, 2, 9)
ax10 = plt.subplot(5, 2, 10)

promoter = '2.5kb'
a = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax.set_title(promoter)

a = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax2.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax.set_title('Igsf9(2.5kb)-eGFP')
ax2.set_title('CAG-tDt')
ax.set_xlabel('Dilution')
ax2.set_xlabel('Dilution')
ax.set_ylabel('IO-spe./non-spe.')

promoter = '1.3kb'
a = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax3.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax3.set_title(promoter)

a = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax4.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax3.set_title('Igsf9(1.3kb)-eGFP')
ax4.set_title('CAG-tDt')
ax3.set_xlabel('Dilution')
ax4.set_xlabel('Dilution')
ax3.set_ylabel('IO-spe./non-spe.')

promoter = '5HT'
a = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax5.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax5.set_title(promoter)

a = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax6.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax5.set_title(promoter+'-eGFP')
ax6.set_title('CAG-tDt')
ax5.set_xlabel('Dilution')
ax6.set_xlabel('Dilution')
ax5.set_ylabel('IO-spe./non-spe.')

promoter = 'PDX'
a = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax7.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax7.set_title(promoter)

a = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax8.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax7.set_title(promoter+'-eGFP')
ax8.set_title('CAG-tDt')
ax7.set_xlabel('Dilution')
ax8.set_xlabel('Dilution')
ax7.set_ylabel('IO-spe./non-spe.')

promoter = 'SUSD4'
a = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [eGFP_HALF_IO[k]/eGFP_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax9.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax9.set_title(promoter)

a = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==1 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
b = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==5 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
c = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==20 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]
d = [tDt_HALF_IO[k]/tDt_HALF_IO_OUT[k] for k in range(len(eGFP_HALF_IO_LABEL)) if eGFP_HALF_IO_DILUTION[k]==200 and(promoter in  eGFP_HALF_IO_LABEL[k])==True]

ax10.boxplot([a,b,c,d], labels=['1', '5', '20', '200'])
ax9.set_title(promoter+'-eGFP')
ax10.set_title('CAG-tDt')
ax9.set_xlabel('Dilution')
ax10.set_xlabel('Dilution')
ax9.set_ylabel('IO-spe./non-spe.')

#5x MAPPING INSIDE / OUTSIDE IO

plt.figure(figsize=(6,3))
s = ['Igsf9(2.5)', 'Igsf9(1.3)','5HTr2b(3.6)','5HTr2b(1.8)','5HTr2b(1.0)','PDX1','SUSD4']

ax1 = plt.subplot(141)
ax2 = plt.subplot(142)
ax3 = plt.subplot(143)
ax4 = plt.subplot(144)

ax1.scatter(FULL_eGFP_IO_IN_ALL, FULL_tDt_IO_IN_ALL)
for i in range(len(s)):
    ax1.text(FULL_eGFP_IO_IN_ALL[i], FULL_tDt_IO_IN_ALL[i], s[i]) 
ax1.set_xlabel('eGFP IO')
ax1.set_ylabel('tdT IO')

ax2.scatter(FULL_eGFP_IO_OUT_ALL, FULL_tDt_IO_OUT_ALL)
for i in range(len(s)):
    ax2.text(FULL_eGFP_IO_OUT_ALL[i], FULL_tDt_IO_OUT_ALL[i], s[i]) 
ax2.set_xlabel('eGFP OUT')
ax2.set_ylabel('tdT OUT')

ax3.scatter(FULL_eGFP_IO_OUT_RATIO_ALL, FULL_tDt_IO_OUT_RATIO_ALL)
for i in range(len(s)):
    ax3.text(FULL_eGFP_IO_OUT_RATIO_ALL[i], FULL_tDt_IO_OUT_RATIO_ALL[i], s[i]) 
ax3.set_xlabel('eGFP (IN/OUT)')
ax3.set_ylabel('tdT (IN/OUT)')

ax4.scatter(FULL_eGFP_IO_IN_RATIO_div_tdT_IN_RATIO_ALL, FULL_eGFP_IO_OUT_div_tdT_OUT_ALL)
for i in range(len(s)):
    ax4.text(FULL_eGFP_IO_IN_RATIO_div_tdT_IN_RATIO_ALL[i], FULL_eGFP_IO_OUT_div_tdT_OUT_ALL[i], s[i]) 
ax4.set_xlabel('Inside IO (eGFP/tDt)')
ax4.set_ylabel('Around IO (eGFP/tDt)')



#60x images _ ASTROCYTE VERSUS NEURON
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
FULL_PDX1_CELLS_tDt = []
FULL_SUSD4_CELLS_tDt = []

FULL_Igsf9_37_DIVIDED = []
FULL_Igsf9_25_DIVIDED = []
FULL_Igsf9_13_DIVIDED = []
FULL_5HTr2b_10_DIVIDED = []
FULL_5HTr2b_18_DIVIDED = []
FULL_5HTr2b_DIVIDED = []
FULL_PDX1_DIVIDED = []
FULL_SUSD4_DIVIDED = []
FULL_Igsf9_37_SUBTRACTED = []
FULL_Igsf9_25_SUBTRACTED = []
FULL_Igsf9_13_SUBTRACTED = []
FULL_5HTr2b_10_SUBTRACTED = []
FULL_5HTr2b_18_SUBTRACTED = []
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

'''
FULL_SUSD4_CELLS = np.concatenate([FULL_SUSD4_CELLS[k]/np.nanmax(FULL_SUSD4_CELLS[k]) for k in range(len(FULL_SUSD4_CELLS))])
FULL_PDX1_CELLS = np.concatenate([FULL_PDX1_CELLS[k]/np.nanmax(FULL_PDX1_CELLS[k]) for k in range(len(FULL_PDX1_CELLS))])
FULL_Igsf9_25_CELLS = np.concatenate([FULL_Igsf9_25_CELLS[k]/np.nanmax(FULL_Igsf9_25_CELLS[k]) for k in range(len(FULL_Igsf9_25_CELLS))])
FULL_Igsf9_13_CELLS = np.concatenate([FULL_Igsf9_13_CELLS[k]/np.nanmax(FULL_Igsf9_13_CELLS[k]) for k in range(len(FULL_Igsf9_13_CELLS))])
FULL_5HTr2b_10_CELLS = np.concatenate([FULL_5HTr2b_10_CELLS[k]/np.nanmax(FULL_5HTr2b_10_CELLS[k]) for k in range(len(FULL_5HTr2b_10_CELLS))])
FULL_5HTr2b_18_CELLS = np.concatenate([FULL_5HTr2b_18_CELLS[k]/np.nanmax(FULL_5HTr2b_18_CELLS[k]) for k in range(len(FULL_5HTr2b_18_CELLS))])
FULL_5HTr2b_CELLS = np.concatenate([FULL_5HTr2b_CELLS[k]/np.nanmax(FULL_5HTr2b_CELLS[k]) for k in range(len(FULL_5HTr2b_CELLS))])
FULL_tDt_CELLS = np.concatenate([FULL_tDt_CELLS[k]/np.nanmax(FULL_tDt_CELLS[k]) for k in range(len(FULL_tDt_CELLS))])

FULL_SUSD4_CELLS_tDt = np.concatenate([FULL_SUSD4_CELLS_tDt[k]/np.nanmax(FULL_SUSD4_CELLS_tDt[k]) for k in range(len(FULL_SUSD4_CELLS_tDt))])
FULL_PDX1_CELLS_tDt = np.concatenate([FULL_PDX1_CELLS_tDt[k]/np.nanmax(FULL_PDX1_CELLS_tDt[k]) for k in range(len(FULL_PDX1_CELLS_tDt))])
FULL_Igsf9_25_CELLS_tDt = np.concatenate([FULL_Igsf9_25_CELLS_tDt[k]/np.nanmax(FULL_Igsf9_25_CELLS_tDt[k]) for k in range(len(FULL_Igsf9_25_CELLS_tDt))])
FULL_Igsf9_13_CELLS_tDt = np.concatenate([FULL_Igsf9_13_CELLS_tDt[k]/np.nanmax(FULL_Igsf9_13_CELLS_tDt[k]) for k in range(len(FULL_Igsf9_13_CELLS_tDt))])
FULL_5HTr2b_10_CELLS_tDt = np.concatenate([FULL_5HTr2b_10_CELLS_tDt[k]/np.nanmax(FULL_5HTr2b_10_CELLS_tDt[k]) for k in range(len(FULL_5HTr2b_10_CELLS_tDt))])
FULL_5HTr2b_18_CELLS_tDt = np.concatenate([FULL_5HTr2b_18_CELLS_tDt[k]/np.nanmax(FULL_5HTr2b_18_CELLS_tDt[k]) for k in range(len(FULL_5HTr2b_18_CELLS_tDt))])
FULL_5HTr2b_CELLS_tDt = np.concatenate([FULL_5HTr2b_CELLS_tDt[k]/np.nanmax(FULL_5HTr2b_CELLS_tDt[k]) for k in range(len(FULL_5HTr2b_CELLS_tDt))])
'''
"""
FULL_SUSD4_CELLS = np.concatenate(FULL_SUSD4_CELLS)
FULL_PDX1_CELLS = np.concatenate(FULL_PDX1_CELLS)
FULL_Igsf9_25_CELLS = np.concatenate(FULL_Igsf9_25_CELLS)
FULL_Igsf9_13_CELLS = np.concatenate(FULL_Igsf9_13_CELLS)
FULL_5HTr2b_10_CELLS = np.concatenate(FULL_5HTr2b_10_CELLS)
FULL_5HTr2b_18_CELLS = np.concatenate(FULL_5HTr2b_18_CELLS)
FULL_5HTr2b_CELLS = np.concatenate(FULL_5HTr2b_CELLS)
FULL_tDt_CELLS = np.concatenate(FULL_tDt_CELLS)

FULL_SUSD4_CELLS_tDt = np.concatenate(FULL_SUSD4_CELLS_tDt)
FULL_PDX1_CELLS_tDt = np.concatenate(FULL_PDX1_CELLS_tDt)
FULL_Igsf9_25_CELLS_tDt = np.concatenate(FULL_Igsf9_25_CELLS_tDt)
FULL_Igsf9_13_CELLS_tDt = np.concatenate(FULL_Igsf9_13_CELLS_tDt)
FULL_5HTr2b_10_CELLS_tDt = np.concatenate(FULL_5HTr2b_10_CELLS_tDt)
FULL_5HTr2b_18_CELLS_tDt = np.concatenate(FULL_5HTr2b_18_CELLS_tDt)
FULL_5HTr2b_CELLS_tDt = np.concatenate(FULL_5HTr2b_CELLS_tDt)
"""
ax.scatter(np.concatenate(FULL_tDt_AREA), np.concatenate(FULL_tDt_CELLS))


ax2.scatter(np.concatenate(FULL_Igsf9_25_AREA), np.concatenate(FULL_Igsf9_25_DIVIDED))
ax3.scatter(np.concatenate(FULL_Igsf9_13_AREA), np.concatenate(FULL_Igsf9_13_DIVIDED))
ax4.scatter(np.concatenate(FULL_5HTr2b_AREA), np.concatenate(FULL_5HTr2b_DIVIDED))
ax5.scatter(np.concatenate(FULL_PDX1_AREA), np.concatenate(FULL_PDX1_DIVIDED))
ax6.scatter(np.concatenate(FULL_SUSD4_AREA), np.concatenate(FULL_SUSD4_DIVIDED))

"""
ax8.scatter(np.concatenate(FULL_Igsf9_25_AREA), np.concatenate(FULL_Igsf9_25_SUBTRACTED))
ax9.scatter(np.concatenate(FULL_Igsf9_13_AREA), np.concatenate(FULL_Igsf9_13_SUBTRACTED))
ax10.scatter(np.concatenate(FULL_5HTr2b_AREA), np.concatenate(FULL_5HTr2b_SUBTRACTED))
ax11.scatter(np.concatenate(FULL_PDX1_AREA), np.concatenate(FULL_PDX1_SUBTRACTED))
ax12.scatter(np.concatenate(FULL_SUSD4_AREA), np.concatenate(FULL_SUSD4_SUBTRACTED))
"""

#ax2.scatter(FULL_Igsf9_25_CELLS,FULL_tDt_CELLS, alpha=0.5, normed=True )
#ax3.hist(FULL_Igsf9_13_CELLS, FULL_tDt_CELLS, alpha=0.5, normed=True )
#ax4.hist(FULL_5HTr2b_CELLS, FULL_tDt_CELLS, alpha=0.5, normed=True )
#ax5.hist(FULL_PDX1_CELLS, FULL_tDt_CELLS, alpha=0.5, normed=True )
#ax6.hist(FULL_SUSD4_CELLS, FULL_tDt_CELLS, alpha=0.5, normed=True )

#ax13.hist(FULL_Igsf9_25_CELLS, histtype='step', alpha=0.5, normed=True, cumulative=True, bins=100 )
#ax13.hist(FULL_Igsf9_13_CELLS, histtype='step', alpha=0.5, normed=True, cumulative=True, bins=100 )
#ax13.hist(FULL_5HTr2b_CELLS, histtype='step', alpha=0.5, normed=True, cumulative=True, bins=100 )
#ax13.hist(FULL_PDX1_CELLS, histtype='step', alpha=0.5, normed=True, cumulative=True, bins=100 )
#ax13.hist(FULL_SUSD4_CELLS, histtype='step', alpha=0.5, normed=True, cumulative=True, bins=100 )
#ax13.hist(FULL_tDt_CELLS, histtype='step', alpha=0.5, color='red', normed=True, cumulative=True, bins=100 )

ASTRO_LABEL = []
NEURON_LABEL = []
ASTRO_VAR = []
NEURON_VAR = []
EXPE_ID = []

ASTROS = [np.concatenate(FULL_Igsf9_25_DIVIDED)[k] for k in range(len(np.concatenate(FULL_Igsf9_25_DIVIDED))) if np.concatenate(FULL_Igsf9_25_AREA)[k]<100]
NEURONS = [np.concatenate(FULL_Igsf9_25_DIVIDED)[k] for k in range(len(np.concatenate(FULL_Igsf9_25_DIVIDED))) if np.concatenate(FULL_Igsf9_25_AREA)[k]>=100]
Igsf9_25_NEURONS = NEURONS
Igsf9_25_ASTROS = ASTROS
Igsf9_25_BOOST = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]>=1.25]
Igsf9_25_LAB = [NEURONS[k] for k in range(len(NEURONS)) if 0.75<NEURONS[k]<1.25]
Igsf9_25_NON_LAB = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]<=0.75]
Igsf9_25_BOOST_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if ASTROS[k]>=1.25]
Igsf9_25_LAB_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if 0.75<ASTROS[k]<1.25]
Igsf9_25_NON_LAB_ASTROS = [NEURONS[k] for k in range(len(ASTROS)) if ASTROS[k]<=0.75]
ASTRO_LABEL.append(np.nanmean(ASTROS))
NEURON_LABEL.append(np.nanmean(NEURONS))

ASTRO_VAR.append(np.nanvar(ASTROS)/np.nanmean(ASTROS))
NEURON_VAR.append(np.nanvar(NEURONS)/np.nanmean(NEURONS))
EXPE_ID.append('Igsf9(2.5)')
ax8.boxplot([ASTROS, NEURONS], labels=['Astro.', 'Neu.'], vert=False)
ax14.hist(NEURONS, histtype='stepfilled', alpha=0.5, normed=True )
ax20.hist(ASTROS, histtype='stepfilled', alpha=0.5, normed=True )

ASTROS = [np.concatenate(FULL_Igsf9_13_DIVIDED)[k] for k in range(len(np.concatenate(FULL_Igsf9_13_DIVIDED))) if np.concatenate(FULL_Igsf9_13_AREA)[k]<100]
NEURONS = [np.concatenate(FULL_Igsf9_13_DIVIDED)[k] for k in range(len(np.concatenate(FULL_Igsf9_13_DIVIDED))) if np.concatenate(FULL_Igsf9_13_AREA)[k]>=100]
Igsf9_13_NEURONS = NEURONS
Igsf9_13_ASTROS = ASTROS
Igsf9_13_BOOST = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]>=1.25]
Igsf9_13_LAB = [NEURONS[k] for k in range(len(NEURONS)) if 0.75<NEURONS[k]<1.25]
Igsf9_13_NON_LAB = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]<=0.75]
Igsf9_13_BOOST_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if ASTROS[k]>=1.25]
Igsf9_13_LAB_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if 0.75<ASTROS[k]<1.25]
Igsf9_13_NON_LAB_ASTROS = [NEURONS[k] for k in range(len(ASTROS)) if ASTROS[k]<=0.75]
ASTRO_LABEL.append(np.nanmean(ASTROS))
NEURON_LABEL.append(np.nanmean(NEURONS))

ASTRO_VAR.append(np.nanvar(ASTROS)/np.nanmean(ASTROS))
NEURON_VAR.append(np.nanvar(NEURONS)/np.nanmean(NEURONS))
EXPE_ID.append('Igsf9(1.3)')
ax9.boxplot([ASTROS, NEURONS], labels=['Astro.', 'Neu.'], vert=False)
ax15.hist(NEURONS, histtype='stepfilled', alpha=0.5, normed=True )
ax21.hist(ASTROS, histtype='stepfilled', alpha=0.5, normed=True )

ASTROS = [np.concatenate(FULL_5HTr2b_10_DIVIDED)[k] for k in range(len(np.concatenate(FULL_5HTr2b_10_DIVIDED))) if np.concatenate(FULL_5HTr2b_10_AREA)[k]<100]
NEURONS = [np.concatenate(FULL_5HTr2b_10_DIVIDED)[k] for k in range(len(np.concatenate(FULL_5HTr2b_10_DIVIDED))) if np.concatenate(FULL_5HTr2b_10_AREA)[k]>=100]
HTR2b_10_NEURONS = NEURONS
HTR2b_10_ASTROS = ASTROS
HTR2b_10_BOOST = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]>=1.25]
HTR2b_10_LAB = [NEURONS[k] for k in range(len(NEURONS)) if 0.75<NEURONS[k]<1.25]
HTR2b_10_NON_LAB = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]<=0.75]
HTR2b_10_BOOST_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if ASTROS[k]>=1.25]
HTR2b_10_LAB_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if 0.75<ASTROS[k]<1.25]
HTR2b_10_NON_LAB_ASTROS = [NEURONS[k] for k in range(len(ASTROS)) if ASTROS[k]<=0.75]
ASTRO_LABEL.append(np.nanmean(ASTROS))
NEURON_LABEL.append(np.nanmean(NEURONS))

ASTRO_VAR.append(np.nanvar(ASTROS)/np.nanmean(ASTROS))
NEURON_VAR.append(np.nanvar(NEURONS)/np.nanmean(NEURONS))

EXPE_ID.append('5Htr2b_10')
#ax10.boxplot([ASTROS, NEURONS], labels=['Astro.', 'Neu.'], vert=False)
#ax16.hist(NEURONS, histtype='stepfilled', alpha=0.5, normed=True )
#ax22.hist(ASTROS, histtype='stepfilled', alpha=0.5, normed=True )

ASTROS = [np.concatenate(FULL_5HTr2b_18_DIVIDED)[k] for k in range(len(np.concatenate(FULL_5HTr2b_18_DIVIDED))) if np.concatenate(FULL_5HTr2b_18_AREA)[k]<100]
NEURONS = [np.concatenate(FULL_5HTr2b_18_DIVIDED)[k] for k in range(len(np.concatenate(FULL_5HTr2b_18_DIVIDED))) if np.concatenate(FULL_5HTr2b_18_AREA)[k]>=100]
HTR2b_18_NEURONS = NEURONS
HTR2b_18_ASTROS = ASTROS
HTR2b_18_BOOST = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]>=1.25]
HTR2b_18_LAB = [NEURONS[k] for k in range(len(NEURONS)) if 0.75<NEURONS[k]<1.25]
HTR2b_18_NON_LAB = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]<=0.75]
HTR2b_18_BOOST_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if ASTROS[k]>=1.25]
HTR2b_18_LAB_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if 0.75<ASTROS[k]<1.25]
HTR2b_18_NON_LAB_ASTROS = [NEURONS[k] for k in range(len(ASTROS)) if ASTROS[k]<=0.75]
ASTRO_LABEL.append(np.nanmean(ASTROS))
NEURON_LABEL.append(np.nanmean(NEURONS))

ASTRO_VAR.append(np.nanvar(ASTROS)/np.nanmean(ASTROS))
NEURON_VAR.append(np.nanvar(NEURONS)/np.nanmean(NEURONS))

EXPE_ID.append('5Htr2b_18')
#ax10.boxplot([ASTROS, NEURONS], labels=['Astro.', 'Neu.'], vert=False)
#ax16.hist(NEURONS, histtype='stepfilled', alpha=0.5, normed=True )
#ax22.hist(ASTROS, histtype='stepfilled', alpha=0.5, normed=True )

ASTROS = [np.concatenate(FULL_5HTr2b_DIVIDED)[k] for k in range(len(np.concatenate(FULL_5HTr2b_DIVIDED))) if np.concatenate(FULL_5HTr2b_AREA)[k]<100]
NEURONS = [np.concatenate(FULL_5HTr2b_DIVIDED)[k] for k in range(len(np.concatenate(FULL_5HTr2b_DIVIDED))) if np.concatenate(FULL_5HTr2b_AREA)[k]>=100]
HTR2b_NEURONS = NEURONS
HTR2b_ASTROS = ASTROS
HTR2b_BOOST = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]>=1.25]
HTR2b_LAB = [NEURONS[k] for k in range(len(NEURONS)) if 0.75<NEURONS[k]<1.25]
HTR2b_NON_LAB = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]<=0.75]
HTR2b_BOOST_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if ASTROS[k]>=1.25]
HTR2b_LAB_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if 0.75<ASTROS[k]<1.25]
HTR2b_NON_LAB_ASTROS = [NEURONS[k] for k in range(len(ASTROS)) if ASTROS[k]<=0.75]
ASTRO_LABEL.append(np.nanmean(ASTROS))
NEURON_LABEL.append(np.nanmean(NEURONS))

ASTRO_VAR.append(np.nanvar(ASTROS)/np.nanmean(ASTROS))
NEURON_VAR.append(np.nanvar(NEURONS)/np.nanmean(NEURONS))

EXPE_ID.append('5Htr2b')
ax10.boxplot([ASTROS, NEURONS], labels=['Astro.', 'Neu.'], vert=False)
ax16.hist(NEURONS, histtype='stepfilled', alpha=0.5, normed=True )
ax22.hist(ASTROS, histtype='stepfilled', alpha=0.5, normed=True )

ASTROS = [np.concatenate(FULL_PDX1_DIVIDED)[k] for k in range(len(np.concatenate(FULL_PDX1_DIVIDED))) if np.concatenate(FULL_PDX1_AREA)[k]<100]
NEURONS = [np.concatenate(FULL_PDX1_DIVIDED)[k] for k in range(len(np.concatenate(FULL_PDX1_DIVIDED))) if np.concatenate(FULL_PDX1_AREA)[k]>=100]
PDX1_NEURONS = NEURONS
PDX1_ASTROS = ASTROS
PDX1_BOOST = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]>=1.25]
PDX1_LAB = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]<1.25]
PDX1_NON_LAB = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]<=0.75]
PDX1_BOOST_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if ASTROS[k]>=1.25]
PDX1_LAB_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if 0.75<ASTROS[k]<1.25]
PDX1_NON_LAB_ASTROS = [NEURONS[k] for k in range(len(ASTROS)) if ASTROS[k]<=0.75]
ASTRO_LABEL.append(np.nanmean(ASTROS))
NEURON_LABEL.append(np.nanmean(NEURONS))

ASTRO_VAR.append(np.nanvar(ASTROS)/np.nanmean(ASTROS))
NEURON_VAR.append(np.nanvar(NEURONS)/np.nanmean(NEURONS))
EXPE_ID.append('PDX1)')
ax11.boxplot([ASTROS, NEURONS], labels=['Astro.', 'Neu.'], vert=False)
ax17.hist(NEURONS, histtype='stepfilled', alpha=0.5, normed=True )
ax23.hist(ASTROS, histtype='stepfilled', alpha=0.5, normed=True )

ASTROS = [np.concatenate(FULL_SUSD4_DIVIDED)[k] for k in range(len(np.concatenate(FULL_SUSD4_DIVIDED))) if np.concatenate(FULL_SUSD4_AREA)[k]<100]
NEURONS = [np.concatenate(FULL_SUSD4_DIVIDED)[k] for k in range(len(np.concatenate(FULL_SUSD4_DIVIDED))) if np.concatenate(FULL_SUSD4_AREA)[k]>=100]
SUSD4_NEURONS = NEURONS
SUSD4_ASTROS = ASTROS
SUSD4_BOOST = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]>=1.25]
SUSD4_LAB = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]<1.25]
SUSD4_NON_LAB = [NEURONS[k] for k in range(len(NEURONS)) if NEURONS[k]<=0.75]
SUSD4_BOOST_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if ASTROS[k]>=1.25]
SUSD4_LAB_ASTROS = [ASTROS[k] for k in range(len(ASTROS)) if 0.75<ASTROS[k]<1.25]
SUSD4_NON_LAB_ASTROS = [NEURONS[k] for k in range(len(ASTROS)) if ASTROS[k]<=0.75]
ASTRO_LABEL.append(np.nanmean(ASTROS))
NEURON_LABEL.append(np.nanmean(NEURONS))

ASTRO_VAR.append(np.nanvar(ASTROS)/np.nanmean(ASTROS))
NEURON_VAR.append(np.nanvar(NEURONS)/np.nanmean(NEURONS))

EXPE_ID.append('SUSD4')
ax12.boxplot([ASTROS, NEURONS], labels=['Astro.', 'Neu.'], vert=False)
ax18.hist(NEURONS, histtype='stepfilled', alpha=0.5, normed=True )
ax24.hist(ASTROS, histtype='stepfilled', alpha=0.5, normed=True )

x, y = np.histogram(np.nan_to_num(np.concatenate(FULL_tDt_AREA)), bins=40)

ax7.plot(y[1::], x)

ax13.scatter(ASTRO_VAR, NEURON_VAR)
for i in range(len(ASTRO_VAR)):
    ax13.text(ASTRO_VAR[i], NEURON_VAR[i], EXPE_ID[i])

ax19.scatter(ASTRO_LABEL, NEURON_LABEL)
for i in range(len(ASTRO_LABEL)):
    ax19.text(ASTRO_LABEL[i], NEURON_LABEL[i], EXPE_ID[i])


#ax14.plot((1, 1), (0, 6), color='black', ls='--')
#ax15.plot((1, 1), (0, 6), color='black', ls='--')
#ax16.plot((1, 1), (0, 6), color='black', ls='--')
#ax17.plot((1, 1), (0, 6), color='black', ls='--')
#ax18.plot((1, 1), (0, 6), color='black', ls='--')
#ax20.plot((1, 1), (0, 6), color='black', ls='--')
#ax21.plot((1, 1), (0, 6), color='black', ls='--')
#ax22.plot((1, 1), (0, 6), color='black', ls='--')
#ax23.plot((1, 1), (0, 6), color='black', ls='--')
#ax24.plot((1, 1), (0, 6), color='black', ls='--')

ax.set_title('CAG-tDt')
ax2.set_title('Igsf9(2.5kb)-eGFP')
ax3.set_title('Igsf9(1.3kb)-eGFP')
ax4.set_title('5HTr2b-eGFP')
ax5.set_title('PDX1-eGFP')
ax6.set_title('SUSD4-eGFP')
ax13.set_xlabel('Astrocytes (VAR/MEAN)')
ax13.set_ylabel('Neurons (VAR/MEAN)')
ax19.set_xlabel('Astrocytes (SUB)')
ax19.set_ylabel('Neurons (SUB)')



non_labeled_SUSD4_NEURONS = len(SUSD4_NON_LAB)/(len(SUSD4_BOOST)+len(SUSD4_LAB)+len(SUSD4_NON_LAB))
non_labeled_PDX1_NEURONS =len(PDX1_NON_LAB)/(len(PDX1_BOOST)+len(PDX1_LAB)+len(PDX1_NON_LAB))
non_labeled_5HTr2b_NEURONS =len(HTR2b_NON_LAB)/(len(HTR2b_BOOST)+len(HTR2b_LAB)+len(HTR2b_NON_LAB))
non_labeled_Igsf9_25_NEURONS = len(Igsf9_25_NON_LAB)/(len(Igsf9_25_BOOST)+len(Igsf9_25_LAB)+len(Igsf9_25_NON_LAB))
non_labeled_Igsf9_13_NEURONS =len(Igsf9_13_NON_LAB)/(len(Igsf9_13_BOOST)+len(Igsf9_13_LAB)+len(Igsf9_13_NON_LAB))


BOOST_SUSD4_NEURONS= len(SUSD4_BOOST)/(len(SUSD4_BOOST)+len(SUSD4_LAB)+len(SUSD4_NON_LAB))
BOOST_PDX1_NEURONS =len(PDX1_BOOST)/(len(PDX1_BOOST)+len(PDX1_LAB)+len(PDX1_NON_LAB))
BOOST_5HTr2b_NEURONS =len(HTR2b_BOOST)/(len(HTR2b_BOOST)+len(HTR2b_LAB)+len(HTR2b_NON_LAB))
BOOST_Igsf9_25_NEURONS = len(Igsf9_25_BOOST)/(len(Igsf9_25_BOOST)+len(Igsf9_25_LAB)+len(Igsf9_25_NON_LAB))
BOOST_Igsf9_13_NEURONS =len(Igsf9_13_BOOST)/(len(Igsf9_13_BOOST)+len(Igsf9_13_LAB)+len(Igsf9_13_NON_LAB))

non_labeled_SUSD4_ASTROS = len(SUSD4_NON_LAB_ASTROS)/(len(SUSD4_BOOST_ASTROS)+len(SUSD4_LAB_ASTROS)+len(SUSD4_NON_LAB_ASTROS))
non_labeled_PDX1_ASTROS =len(PDX1_NON_LAB_ASTROS)/(len(PDX1_BOOST_ASTROS)+len(PDX1_LAB_ASTROS)+len(PDX1_NON_LAB_ASTROS))
non_labeled_5HTr2b_ASTROS =len(HTR2b_NON_LAB_ASTROS)/(len(HTR2b_BOOST_ASTROS)+len(HTR2b_LAB_ASTROS)+len(HTR2b_NON_LAB_ASTROS))
non_labeled_Igsf9_25_ASTROS = len(Igsf9_25_NON_LAB_ASTROS)/(len(Igsf9_25_BOOST_ASTROS)+len(Igsf9_25_LAB_ASTROS)+len(Igsf9_25_NON_LAB_ASTROS))
non_labeled_Igsf9_13_ASTROS =len(Igsf9_13_NON_LAB_ASTROS)/(len(Igsf9_13_BOOST_ASTROS)+len(Igsf9_13_LAB_ASTROS)+len(Igsf9_13_NON_LAB_ASTROS))


BOOST_SUSD4_ASTROS = len(SUSD4_BOOST_ASTROS)/(len(SUSD4_BOOST_ASTROS)+len(SUSD4_LAB_ASTROS)+len(SUSD4_NON_LAB_ASTROS))
BOOST_PDX1_ASTROS =len(PDX1_BOOST_ASTROS)/(len(PDX1_BOOST_ASTROS)+len(PDX1_LAB_ASTROS)+len(PDX1_NON_LAB_ASTROS))
BOOST_5HTr2b_ASTROS =len(HTR2b_BOOST_ASTROS)/(len(HTR2b_BOOST_ASTROS)+len(HTR2b_LAB)+len(HTR2b_NON_LAB_ASTROS))
BOOST_Igsf9_25_ASTROS = len(Igsf9_25_BOOST_ASTROS)/(len(Igsf9_25_BOOST_ASTROS)+len(Igsf9_25_LAB_ASTROS)+len(Igsf9_25_NON_LAB_ASTROS))
BOOST_Igsf9_13_ASTROS =len(Igsf9_13_BOOST_ASTROS)/(len(Igsf9_13_BOOST_ASTROS)+len(Igsf9_13_LAB_ASTROS)+len(Igsf9_13_NON_LAB_ASTROS))



plt.figure(figsize=(4,4))
ax = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)
ax.bar(height = [non_labeled_Igsf9_25_ASTROS, non_labeled_Igsf9_13_ASTROS, non_labeled_5HTr2b_ASTROS, non_labeled_PDX1_ASTROS, non_labeled_SUSD4_ASTROS], x=['Igsf9(2.5)', 'Igsf9(1.3)','5HTr2b','PDX1','SUSD4'])
ax2.bar(height =[non_labeled_Igsf9_25_NEURONS, non_labeled_Igsf9_13_NEURONS, non_labeled_5HTr2b_NEURONS, non_labeled_PDX1_NEURONS, non_labeled_SUSD4_NEURONS], x=['Igsf9(2.5)', 'Igsf9(1.3)','5HTr2b','PDX1','SUSD4'])
ax3.bar(height =[BOOST_Igsf9_25_ASTROS, BOOST_Igsf9_13_ASTROS, BOOST_5HTr2b_ASTROS, BOOST_PDX1_ASTROS, BOOST_SUSD4_ASTROS], x=['Igsf9(2.5)', 'Igsf9(1.3)','5HTr2b','PDX1','SUSD4'])
ax4.bar(height =[BOOST_Igsf9_25_NEURONS, BOOST_Igsf9_13_NEURONS, BOOST_5HTr2b_NEURONS, BOOST_PDX1_NEURONS, BOOST_SUSD4_NEURONS], x=['Igsf9(2.5)', 'Igsf9(1.3)','5HTr2b','PDX1','SUSD4'])

ax.set_title('ASTROS')
ax2.set_title('NEURONS')
ax.set_ylabel('No-label')
ax3.set_ylabel('Boost')
plt.tight_layout()


plt.figure(figsize=(10,3))
ax = plt.subplot(121)
DATA = pd.DataFrame()
DATA['Igsf9_25']=pd.Series(Igsf9_25_NEURONS)
DATA['Igsf9_13']=pd.Series(Igsf9_13_NEURONS)
DATA['5HTr2b']=pd.Series(HTR2b_NEURONS)
DATA['5HTr2b_10']=pd.Series(HTR2b_10_NEURONS)
DATA['5HTr2b_18']=pd.Series(HTR2b_18_NEURONS)
DATA['PDX1']=pd.Series(PDX1_NEURONS)
DATA['SUSD4']=pd.Series(SUSD4_NEURONS)
sns.violinplot(data = [DATA['Igsf9_25'].values, DATA['Igsf9_13'].values, DATA['5HTr2b_10'].values, DATA['5HTr2b_18'].values, DATA['5HTr2b'].values, DATA['PDX1'].values, DATA['SUSD4'].values], palette='Set3')

ax2 = plt.subplot(122)
DATA_ASTROS = pd.DataFrame()
DATA_ASTROS['Igsf9_25']=pd.Series(Igsf9_25_ASTROS)
DATA_ASTROS['Igsf9_13']=pd.Series(Igsf9_13_ASTROS)
DATA_ASTROS['5HTr2b_10']=pd.Series(HTR2b_10_ASTROS)
DATA_ASTROS['5HTr2b_18']=pd.Series(HTR2b_18_ASTROS)
DATA_ASTROS['5HTr2b']=pd.Series(HTR2b_ASTROS)
DATA_ASTROS['PDX1']=pd.Series(PDX1_ASTROS)
DATA_ASTROS['SUSD4']=pd.Series(SUSD4_ASTROS)
sns.violinplot(data = [DATA_ASTROS['Igsf9_25'].values, DATA_ASTROS['Igsf9_13'].values, DATA_ASTROS['5HTr2b_10'].values, DATA_ASTROS['5HTr2b_18'].values, DATA_ASTROS['5HTr2b'].values, DATA_ASTROS['PDX1'].values, DATA_ASTROS['SUSD4'].values], palette='Set3')

ax.plot((-1, 5), (1, 1), color='black', ls='--')
ax2.plot((-1, 5), (1, 1), color='black', ls='--')
ax.set_ylabel('Expression level (CAG-tDt ratio)')
ax.set_title('Neurons')
ax2.set_title('Astrocytes')




"""
202012
60x DATA VISU 2

plt.figure()
ax1 = plt.subplot(381)
ax2 = plt.subplot(382, sharex=ax1, sharey=ax1)
ax3 = plt.subplot(383, sharex=ax1, sharey=ax1)
ax4 = plt.subplot(384, sharex=ax1, sharey=ax1)
ax5 = plt.subplot(385, sharex=ax1, sharey=ax1)
ax6 = plt.subplot(386, sharex=ax1, sharey=ax1)
ax7 = plt.subplot(387, sharex=ax1, sharey=ax1)
ax8 = plt.subplot(388, sharex=ax1, sharey=ax1)


ax9 = plt.subplot(3, 8, 9)
ax10 = plt.subplot(3, 8, 10, sharex= ax9)
ax11 = plt.subplot(3, 8, 11, sharex= ax9)
ax12 = plt.subplot(3, 8, 12, sharex= ax9)
ax13 = plt.subplot(3, 8, 13, sharex= ax9)
ax14 = plt.subplot(3, 8, 14, sharex= ax9)
ax15 = plt.subplot(3, 8, 15, sharex= ax9)
ax16 = plt.subplot(3, 8, 16, sharex= ax9)

ax17 = plt.subplot(3, 8, 17)
ax18 = plt.subplot(3, 8, 18, sharex= ax17)
ax19 = plt.subplot(3, 8, 19, sharex= ax17)
ax20 = plt.subplot(3, 8, 20, sharex= ax17)
ax21 = plt.subplot(3, 8, 21, sharex= ax17)
ax22 = plt.subplot(3, 8, 22, sharex= ax17)
ax23 = plt.subplot(3, 8, 23, sharex= ax17)
ax24 = plt.subplot(3, 8, 24, sharex= ax17)


all_x = []
all_y = []

for j in range(len(FULL_PDX1_CELLS)):
    to_plot_x = [FULL_PDX1_SUBTRACTED[j][i] for i in range(len(FULL_PDX1_SUBTRACTED[j])) if FULL_PDX1_AREA[j][i]>100]
    to_plot_y = [FULL_PDX1_CELLS[j][i] for i in range(len(FULL_PDX1_CELLS[j])) if FULL_PDX1_AREA[j][i]>100]
    ax1.scatter(to_plot_x, to_plot_y, facecolors='none', edgecolors='black')
    r_score = sp.stats.pearsonr(to_plot_x, to_plot_y)
    ax1.text(0, 0, r_score[0])
    ax9.hist(to_plot_x, density=True, stacked=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_x)
    ax9.scatter(mu, 0, facecolors='none', edgecolors='black')  
    ax9.plot((mu + sig, mu-sig), (0, 0), color='black')
    all_x.append(to_plot_x)
    all_y.append(to_plot_y)
    
    ax17.hist(to_plot_y, density=True, stacked=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_y)
    ax17.scatter(mu, 0, facecolors='none', edgecolors='black')  
    ax17.plot((mu + sig, mu-sig), (0, 0), color='black')
    

for j in range(len(FULL_SUSD4_CELLS)):
    to_plot_x = [FULL_SUSD4_SUBTRACTED[j][i] for i in range(len(FULL_SUSD4_SUBTRACTED[j])) if FULL_SUSD4_AREA[j][i]>100]
    to_plot_y = [FULL_SUSD4_CELLS[j][i] for i in range(len(FULL_SUSD4_CELLS[j])) if FULL_SUSD4_AREA[j][i]>100]
    ax2.scatter(to_plot_x, to_plot_y, facecolors='none', edgecolors='black')
    r_score = sp.stats.pearsonr(to_plot_x, to_plot_y)
    ax2.text(0, 0, r_score[0])
    ax10.hist(to_plot_x, density=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_x)
    ax10.scatter(mu, 0, facecolors='none', edgecolors='black')
    ax10.plot((mu + sig, mu-sig), (0, 0), color='black')
    all_x.append(to_plot_x)
    all_y.append(to_plot_y)
    
    ax18.hist(to_plot_y, density=True, stacked=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_y)
    ax18.scatter(mu, 0, facecolors='none', edgecolors='black')  
    ax18.plot((mu + sig, mu-sig), (0, 0), color='black')    

for j in range(len(FULL_5HTr2b_10_CELLS)):
    to_plot_x = [FULL_5HTr2b_10_SUBTRACTED[j][i] for i in range(len(FULL_5HTr2b_10_SUBTRACTED[j])) if FULL_5HTr2b_10_AREA[j][i]>100]
    to_plot_y = [FULL_5HTr2b_10_CELLS[j][i] for i in range(len(FULL_5HTr2b_10_CELLS[j])) if FULL_5HTr2b_10_AREA[j][i]>100]
    ax3.scatter(to_plot_x, to_plot_y, facecolors='none', edgecolors='black')
    r_score = sp.stats.pearsonr(to_plot_x, to_plot_y)
    ax3.text(0, 0, r_score[0])
    ax11.hist(to_plot_x, density=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_x)
    ax11.scatter(mu, 0, facecolors='none', edgecolors='black')
    ax11.plot((mu + sig, mu-sig), (0, 0), color='black')
    all_x.append(to_plot_x)
    all_y.append(to_plot_y)
    
    ax19.hist(to_plot_y, density=True, stacked=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_y)
    ax19.scatter(mu, 0, facecolors='none', edgecolors='black')  
    ax19.plot((mu + sig, mu-sig), (0, 0), color='black')

for j in range(len(FULL_5HTr2b_18_CELLS)):
    to_plot_x = [FULL_5HTr2b_18_SUBTRACTED[j][i] for i in range(len(FULL_5HTr2b_18_SUBTRACTED[j])) if FULL_5HTr2b_18_AREA[j][i]>100]
    to_plot_y = [FULL_5HTr2b_18_CELLS[j][i] for i in range(len(FULL_5HTr2b_18_CELLS[j])) if FULL_5HTr2b_18_AREA[j][i]>100]
    ax4.scatter(to_plot_x, to_plot_y, facecolors='none', edgecolors='black')
    r_score = sp.stats.pearsonr(to_plot_x, to_plot_y)
    ax4.text(0, 0, r_score[0])
    ax12.hist(to_plot_x, density=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_x)
    ax12.scatter(mu, 0, facecolors='none', edgecolors='black')
    ax12.plot((mu + sig, mu-sig), (0, 0), color='black')
    all_x.append(to_plot_x)
    all_y.append(to_plot_y)
    
    ax20.hist(to_plot_y, density=True, stacked=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_y)
    ax20.scatter(mu, 0, facecolors='none', edgecolors='black')  
    ax20.plot((mu + sig, mu-sig), (0, 0), color='black')


for j in range(len(FULL_5HTr2b_CELLS)):
    to_plot_x = [FULL_5HTr2b_SUBTRACTED[j][i] for i in range(len(FULL_5HTr2b_SUBTRACTED[j])) if FULL_5HTr2b_AREA[j][i]>100]
    to_plot_y = [FULL_5HTr2b_CELLS[j][i] for i in range(len(FULL_5HTr2b_CELLS[j])) if FULL_5HTr2b_AREA[j][i]>100]
    ax5.scatter(to_plot_x, to_plot_y, facecolors='none', edgecolors='black')
    r_score = sp.stats.pearsonr(to_plot_x, to_plot_y)
    ax5.text(0, 0, r_score[0])
    ax13.hist(to_plot_x, density=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_x)
    ax13.scatter(mu, 0, facecolors='none', edgecolors='black')
    ax13.plot((mu + sig, mu-sig), (0, 0), color='black')
    all_x.append(to_plot_x)
    all_y.append(to_plot_y)
    
    ax21.hist(to_plot_y, density=True, stacked=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_y)
    ax21.scatter(mu, 0, facecolors='none', edgecolors='black')  
    ax21.plot((mu + sig, mu-sig), (0, 0), color='black')

for j in range(len(FULL_Igsf9_13_CELLS)):
    to_plot_x = [FULL_Igsf9_13_SUBTRACTED[j][i] for i in range(len(FULL_Igsf9_13_SUBTRACTED[j])) if FULL_Igsf9_13_AREA[j][i]>100]
    to_plot_y = [FULL_Igsf9_13_CELLS[j][i] for i in range(len(FULL_Igsf9_13_CELLS[j])) if FULL_Igsf9_13_AREA[j][i]>100]
    ax6.scatter(to_plot_x, to_plot_y, facecolors='none', edgecolors='black')
    r_score = sp.stats.pearsonr(to_plot_x, to_plot_y)
    ax6.text(0, 0, r_score[0])
    ax14.hist(to_plot_x, density=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_x)
    ax14.scatter(mu, 0, facecolors='none', edgecolors='black')
    ax14.plot((mu + sig, mu-sig), (0, 0), color='black')
    all_x.append(to_plot_x)
    all_y.append(to_plot_y)

    ax22.hist(to_plot_y, density=True, stacked=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_y)
    ax22.scatter(mu, 0, facecolors='none', edgecolors='black')  
    ax22.plot((mu + sig, mu-sig), (0, 0), color='black')

for j in range(len(FULL_Igsf9_25_CELLS)):   
    to_plot_x = [FULL_Igsf9_25_SUBTRACTED[j][i] for i in range(len(FULL_Igsf9_25_SUBTRACTED[j])) if FULL_Igsf9_25_AREA[j][i]>100]
    to_plot_y = [FULL_Igsf9_25_CELLS[j][i] for i in range(len(FULL_Igsf9_25_CELLS[j])) if FULL_Igsf9_25_AREA[j][i]>100]
    ax7.scatter(to_plot_x, to_plot_y, facecolors='none', edgecolors='black')
    r_score = sp.stats.pearsonr(to_plot_x, to_plot_y)
    ax7.text(0, 0, r_score[0])
    ax15.hist(to_plot_x, density=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_x)
    ax15.scatter(mu, 0, facecolors='none', edgecolors='black')
    ax15.plot((mu + sig, mu-sig), (0, 0), color='black')
    all_x.append(to_plot_x)
    all_y.append(to_plot_y)
    
    ax23.hist(to_plot_y, density=True, stacked=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_y)
    ax23.scatter(mu, 0, facecolors='none', edgecolors='black')  
    ax23.plot((mu + sig, mu-sig), (0, 0), color='black')

for j in range(len(FULL_Igsf9_37_CELLS)):   
    to_plot_x = [FULL_Igsf9_37_SUBTRACTED[j][i] for i in range(len(FULL_Igsf9_37_SUBTRACTED[j])) if FULL_Igsf9_37_AREA[j][i]>100]
    to_plot_y = [FULL_Igsf9_37_CELLS[j][i] for i in range(len(FULL_Igsf9_37_CELLS[j])) if FULL_Igsf9_37_AREA[j][i]>100]
    ax8.scatter(to_plot_x, to_plot_y, facecolors='none', edgecolors='black')
    r_score = sp.stats.pearsonr(to_plot_x, to_plot_y)
    ax8.text(0, 0, r_score[0])
    ax16.hist(to_plot_x, density=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_x)
    ax16.scatter(mu, 0, facecolors='none', edgecolors='black')
    ax16.plot((mu + sig, mu-sig), (0, 0), color='black')
    all_x.append(to_plot_x)
    all_y.append(to_plot_y)
    
    ax24.hist(to_plot_y, density=True, stacked=True, histtype='step')
    mu, sig = sp.stats.norm.fit(to_plot_y)
    ax24.scatter(mu, 0, facecolors='none', edgecolors='black')  
    ax24.plot((mu + sig, mu-sig), (0, 0), color='black')

ax1.set_title('PDX1-eGFP')
ax2.set_title('SUSD4-eGFP')
ax3.set_title('5HTr2b(1.0kb)-eGFP')
ax4.set_title('5HTr2b(1.8kb)-eGFP')
ax5.set_title('5HTr2b(3.7kb)-eGFP')
ax6.set_title('Igsf9(1.3kb)-eGFP')
ax7.set_title('Igsf9(2.5kb)-eGFP')
ax8.set_title('Igsf9(3.7kb)-eGFP')

ax1.set_ylabel('Lab. Homogeneity')



"""