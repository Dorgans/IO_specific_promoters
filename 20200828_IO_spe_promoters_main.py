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
import file_functions
import matplotlib.backends.backend_pdf
import seaborn as sns
from PIL import image

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
            
            FULL_RAW_HIST_eGFP.append(eGFP)
            FULL_RAW_HIST_tDt.append(tdTomato)
            
            SUB = np.subtract(np.nan_to_num(eGFP), np.nan_to_num(tdTomato))
            
            ax4.imshow(np.reshape(SUB, (-1, 1024)), cmap=cm.magma,interpolation='hamming')
            
            
            
            IM_ARRAY = np.array(np.reshape(eGFP, (-1, w)))
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
            Z_SCORE = 0#np.nanmean(IM_ARRAY) + 3*np.nanstd(IM_ARRAY)
            eGFP_ARRAY_CONC = [num for num in eGFP_ARRAY_CONC if num>Z_SCORE]
            eGFP_IO = [num for num in eGFP_IO if num>Z_SCORE]
            
            eGFP_OUT_sup = [eGFP_ARRAY_CONC[k] for k in range(len(eGFP_ARRAY_CONC)) if eGFP_ARRAY_CONC[k]>Z_SCORE and  k<int(len(eGFP_ARRAY_CONC)*0.4)]
            eGFP_OUT_inf = [eGFP_ARRAY_CONC[k] for k in range(len(eGFP_ARRAY_CONC)) if eGFP_ARRAY_CONC[k]>Z_SCORE and  k>int(len(eGFP_ARRAY_CONC)*0.6)]
            
            PixNumeGFP.append([len(eGFP_OUT_sup), len(eGFP_OUT_inf)])
            
            IM_ARRAY = np.array(np.reshape(tdTomato, (-1, w)))
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
            Z_SCORE = 0#np.nanmean(IM_ARRAY) + 3*np.nanstd(IM_ARRAY)
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


#60x images _ ASTROCYTE VERSUS NEURON
PATHS = load_directory_content__()

ALL_AREAS = []

FULL_tDt_CELLS = []
FULL_tDt_AREA = []
FULL_Igsf9_25_CELLS = []
FULL_Igsf9_25_AREA = []
FULL_Igsf9_13_CELLS = []
FULL_Igsf9_13_AREA = []
FULL_5HTr2b_CELLS = []
FULL_5HTr2b_AREA = []
FULL_PDX1_CELLS = []
FULL_PDX1_AREA = []
FULL_SUSD4_CELLS = []
FULL_SUSD4_AREA = []

FULL_Igsf9_25_SUBTRACTED = []
FULL_Igsf9_13_SUBTRACTED = []
FULL_5HTr2b_SUBTRACTED = []
FULL_PDX1_SUBTRACTED = []
FULL_SUSD4_SUBTRACTED = []

for i in range(len(PATHS[0])):
    file = PATHS[0][i]
    
    if ('KD' in file)==True and ('60x' in file)==True and ('ROI' in file)==False:
        print('OK')
        Cell_file = pd.read_csv(file)
        CELL_AREA = Cell_file.values[:,1]
        eGFP_intensity = Cell_file.values[:,2]
        tDt_intensity = Cell_file.values[:,3]
        SUB_CHANNELS = Cell_file.values[:,5]
        

        FULL_tDt_CELLS.append(tDt_intensity)
        FULL_tDt_AREA.append(CELL_AREA)
        
        if ('SUSD4' in file)==True:
            FULL_SUSD4_CELLS.append(eGFP_intensity)
            FULL_SUSD4_AREA.append(CELL_AREA)
            FULL_SUSD4_SUBTRACTED.append(SUB_CHANNELS)
        elif ('PDX1' in file)==True:
            FULL_PDX1_CELLS.append(eGFP_intensity)
            FULL_PDX1_AREA.append(CELL_AREA)
            FULL_PDX1_SUBTRACTED.append(SUB_CHANNELS)
        elif ('5HTr2b' in file)==True:
            FULL_5HTr2b_CELLS.append(eGFP_intensity)
            FULL_5HTr2b_AREA.append(CELL_AREA)
            FULL_5HTr2b_SUBTRACTED.append(SUB_CHANNELS)
        elif ('Igsf9_1.3' in file)==True:
            FULL_Igsf9_13_CELLS.append(eGFP_intensity)
            FULL_Igsf9_13_AREA.append(CELL_AREA)
            FULL_Igsf9_13_SUBTRACTED.append(SUB_CHANNELS)
        elif ('Igsf9_2.5' in file)==True and ('2.5' in file)==True:
            FULL_Igsf9_25_CELLS.append(eGFP_intensity)
            FULL_Igsf9_25_AREA.append(CELL_AREA)
            FULL_Igsf9_25_SUBTRACTED.append(SUB_CHANNELS)
            
plt.figure(figsize=(14,7))
ax = plt.subplot(461)
ax2 = plt.subplot(462, sharex=ax, sharey=ax)
ax3 = plt.subplot(463, sharex=ax, sharey=ax)
ax4 = plt.subplot(464, sharex=ax, sharey=ax)
ax5 = plt.subplot(465, sharex=ax, sharey=ax)
ax6 = plt.subplot(466, sharex=ax, sharey=ax)

ax7 = plt.subplot(467)
ax8 = plt.subplot(468, sharex=ax)
ax9 = plt.subplot(469, sharex=ax, sharey=ax8)
ax10 = plt.subplot(4, 6, 10, sharex=ax, sharey=ax8)
ax11 = plt.subplot(4,6 ,11, sharex=ax, sharey=ax8)
ax12 = plt.subplot(4, 6, 12, sharex=ax, sharey=ax8)

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

ax.scatter(np.concatenate(FULL_tDt_AREA), np.concatenate(FULL_tDt_CELLS))
ax2.scatter(np.concatenate(FULL_Igsf9_25_AREA), np.concatenate(FULL_Igsf9_25_CELLS))
ax3.scatter(np.concatenate(FULL_Igsf9_13_AREA), np.concatenate(FULL_Igsf9_13_CELLS))
ax4.scatter(np.concatenate(FULL_5HTr2b_AREA), np.concatenate(FULL_5HTr2b_CELLS))
ax5.scatter(np.concatenate(FULL_PDX1_AREA), np.concatenate(FULL_PDX1_CELLS))
ax6.scatter(np.concatenate(FULL_SUSD4_AREA), np.concatenate(FULL_SUSD4_CELLS))

ax8.scatter(np.concatenate(FULL_Igsf9_25_AREA), np.concatenate(FULL_Igsf9_25_SUBTRACTED))
ax9.scatter(np.concatenate(FULL_Igsf9_13_AREA), np.concatenate(FULL_Igsf9_13_SUBTRACTED))
ax10.scatter(np.concatenate(FULL_5HTr2b_AREA), np.concatenate(FULL_5HTr2b_SUBTRACTED))
ax11.scatter(np.concatenate(FULL_PDX1_AREA), np.concatenate(FULL_PDX1_SUBTRACTED))
ax12.scatter(np.concatenate(FULL_SUSD4_AREA), np.concatenate(FULL_SUSD4_SUBTRACTED))

ASTROS = [np.concatenate(FULL_Igsf9_25_SUBTRACTED)[k] for k in range(len(np.concatenate(FULL_Igsf9_25_SUBTRACTED))) if np.concatenate(FULL_Igsf9_25_AREA)[k]<100]
NEURONS = [np.concatenate(FULL_Igsf9_25_SUBTRACTED)[k] for k in range(len(np.concatenate(FULL_Igsf9_25_SUBTRACTED))) if np.concatenate(FULL_Igsf9_25_AREA)[k]>=100]
ax14.boxplot([ASTROS, NEURONS], labels=['Astro.', 'Neu.'], vert=False)
ax20.hist([ASTROS, NEURONS], histtype='stepfilled', alpha=0.5 )

ASTROS = [np.concatenate(FULL_Igsf9_13_SUBTRACTED)[k] for k in range(len(np.concatenate(FULL_Igsf9_13_SUBTRACTED))) if np.concatenate(FULL_Igsf9_13_AREA)[k]<100]
NEURONS = [np.concatenate(FULL_Igsf9_25_SUBTRACTED)[k] for k in range(len(np.concatenate(FULL_Igsf9_13_SUBTRACTED))) if np.concatenate(FULL_Igsf9_13_AREA)[k]>=100]
ax15.boxplot([ASTROS, NEURONS], labels=['Astro.', 'Neu.'], vert=False)
ax21.hist([ASTROS, NEURONS], histtype='stepfilled', alpha=0.5 )

ASTROS = [np.concatenate(FULL_5HTr2b_SUBTRACTED)[k] for k in range(len(np.concatenate(FULL_5HTr2b_SUBTRACTED))) if np.concatenate(FULL_5HTr2b_AREA)[k]<100]
NEURONS = [np.concatenate(FULL_5HTr2b_SUBTRACTED)[k] for k in range(len(np.concatenate(FULL_5HTr2b_SUBTRACTED))) if np.concatenate(FULL_5HTr2b_AREA)[k]>=100]
ax16.boxplot([ASTROS, NEURONS], labels=['Astro.', 'Neu.'], vert=False)
ax22.hist([ASTROS, NEURONS], histtype='stepfilled', alpha=0.5 )

ASTROS = [np.concatenate(FULL_PDX1_SUBTRACTED)[k] for k in range(len(np.concatenate(FULL_PDX1_SUBTRACTED))) if np.concatenate(FULL_PDX1_AREA)[k]<100]
NEURONS = [np.concatenate(FULL_PDX1_SUBTRACTED)[k] for k in range(len(np.concatenate(FULL_PDX1_SUBTRACTED))) if np.concatenate(FULL_PDX1_AREA)[k]>=100]
ax17.boxplot([ASTROS, NEURONS], labels=['Astro.', 'Neu.'], vert=False)
ax23.hist([ASTROS, NEURONS], histtype='stepfilled', alpha=0.5 )

ASTROS = [np.concatenate(FULL_SUSD4_SUBTRACTED)[k] for k in range(len(np.concatenate(FULL_SUSD4_SUBTRACTED))) if np.concatenate(FULL_SUSD4_AREA)[k]<100]
NEURONS = [np.concatenate(FULL_SUSD4_SUBTRACTED)[k] for k in range(len(np.concatenate(FULL_SUSD4_SUBTRACTED))) if np.concatenate(FULL_SUSD4_AREA)[k]>=100]
ax18.boxplot([ASTROS, NEURONS], labels=['Astro.', 'Neu.'], vert=False)
ax24.hist([ASTROS, NEURONS], histtype='stepfilled', alpha=0.5 )


x, y = np.histogram(np.nan_to_num(np.concatenate(FULL_tDt_AREA)), bins=40)
ax7.plot(y[1::], x)

ax.set_title('CAG-tDt')
ax2.set_title('Igsf9(2.5kb)-eGFP')
ax3.set_title('Igsf9(1.3kb)-eGFP')
ax4.set_title('5HTr2b-eGFP')
ax5.set_title('PDX1-eGFP')
ax6.set_title('SUSD4-eGFP')