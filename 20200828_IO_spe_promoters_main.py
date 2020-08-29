# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 15:38:21 2020

@author: dorga
"""

#IMAGE ANALYSIS PLOT
import matplotlib.backends.backend_pdf

pdf = matplotlib.backends.backend_pdf.PdfPages(r'C:\Users\dorga\Desktop\output.pdf')
DIR__, DIRECTORY = load_directory_content__()

HALF_LABEL_LIST = pd.read_csv(DIRECTORY+r'\20200828_HALF_IO_DILUTION_LABELS.csv')

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
            tdTomato = tdTomato/(np.nanmax(tdTomato)*1)
            eGFP = eGFP/(np.nanmax(eGFP)*1)
            tdTomato = tdTomato-np.nanmin(tdTomato)
            eGFP = eGFP-np.nanmin(eGFP)
            
            
            SUB = np.subtract(eGFP, tdTomato)
            SUB = SUB/(np.nanmax(SUB)*1)
            SUB = SUB-(np.nanmin(SUB)*1)
            
            SUB = np.multiply(SUB, eGFP)
            ax4.imshow(np.reshape(SUB, (-1, 1024)), vmin=0, vmax=0.5, cmap=cm.magma,interpolation='hamming')
            
           
            
            IM_ARRAY = np.array(SPLIT_IMAGE[1])
            ROI_crop = pd.read_csv(PATH.split('tif')[0]+'csv')
            eGFP_IO = []
            eGFP_IO_sup = []
            eGFP_IO_inf = []
            eGFP_OUT_sup = []
            eGFP_OUT_inf = []
            for k in range(len(ROI_crop.X)):
                eGFP_IO.append(IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]])
                if ROI_crop.Y[k]<w*0.4:
                    eGFP_IO_sup.append(IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]])
                elif ROI_crop.Y[k]>w*0.6:
                    eGFP_IO_inf.append(IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]])
                IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]] = 0
            eGFP_OUT_sup = IM_ARRAY[int(len(IM_ARRAY)*0.4):]
            eGFP_OUT_inf =IM_ARRAY[int(len(IM_ARRAY)*0.6):]
            eGFP_ARRAY_CONC = np.concatenate(IM_ARRAY)
            Z_SCORE = np.nanmin(eGFP_ARRAY_CONC)+2*np.std(eGFP_ARRAY_CONC)
            eGFP_ARRAY_CONC = [num for num in eGFP_ARRAY_CONC if num>Z_SCORE]
            eGFP_IO = [num for num in eGFP_IO if num>Z_SCORE]
            
            IM_ARRAY = np.array(SPLIT_IMAGE[0])
            ROI_crop = pd.read_csv(PATH.split('tif')[0]+'csv')
            tDTomato_IO = []
            tDTomato_IO_sup = []
            tDTomato_IO_inf = []
            tDTomato_OUT_sup = []
            tDTomato_OUT_inf = []
            for k in range(len(ROI_crop.X)):
                tDTomato_IO.append(IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]])
                if ROI_crop.Y[k]<w*0.4:
                    tDTomato_IO_sup.append(IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]])
                elif ROI_crop.Y[k]>w*0.6:
                    tDTomato_IO_inf.append(IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]])
                IM_ARRAY[ROI_crop.Y[k]][ROI_crop.X[k]] = 0
            tDTomato_OUT_sup = IM_ARRAY[int(len(IM_ARRAY)*0.4):]
            tDTomato_OUT_inf =IM_ARRAY[int(len(IM_ARRAY)*0.6):] 
            tDt_CONC = np.concatenate(IM_ARRAY)
            Z_SCORE = np.nanmin(tDt_CONC)+2*np.std(tDt_CONC)
            tDt_CONC = [num for num in tDt_CONC if num>Z_SCORE]
            tDTomato_IO = [num for num in tDTomato_IO if num>Z_SCORE] 

            SUB = np.multiply(SUB, 1)
            ax5.boxplot([tDt_CONC, eGFP_ARRAY_CONC, tDTomato_IO, eGFP_IO], showfliers=False)
            
            promoter_name_ = PATH.split('_')[-8]

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

            FULL_tDt_IO_IN_UP.append(np.nanmean(tDTomato_sup))
            FULL_eGFP_IO_IN_UP.append(np.nanmean(eGFP_IO_sup))
            FULL_tDt_IO_IN_DOWN.append(np.nanmean(tDTomato_inf))
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

#EXPRESSION PLOT
plt.figure(figsize=(15,2))
ax = plt.subplot(171)
ax2 = plt.subplot(172, sharex=ax, sharey=ax)
ax3 = plt.subplot(173, sharex=ax, sharey=ax)
ax4 = plt.subplot(174, sharex=ax, sharey=ax)
ax5 = plt.subplot(175)
ax6 = plt.subplot(176)
ax7 = plt.subplot(177)

for i in range(len(FULL_PROMOTER_NAMES)):
    if ('Igsf' in FULL_PROMOTER_NAMES[i])==True:
        ax.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='blue')
        #ax5.scatter('Igsf', FULL_eGFP_IO_OUT_RATIO[i], color='blue')
    elif ('5HTr' in FULL_PROMOTER_NAMES[i])==True:
        ax2.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='red')
        #ax5.scatter('5HTr', FULL_eGFP_IO_OUT_RATIO[i], color='red')
    elif ('PDX' in FULL_PROMOTER_NAMES[i])==True:
        ax3.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='purple')
        #ax5.scatter('PDX1', FULL_eGFP_IO_OUT_RATIO[i], color='purple')
    elif ('SUS' in FULL_PROMOTER_NAMES[i])==True:
        ax4.scatter(FULL_tDt_IO_OUT_RATIO[i], FULL_eGFP_IO_OUT_RATIO[i], color='orange')
        #ax5.scatter('SUSD4', FULL_eGFP_IO_OUT_RATIO[i], color='orange')
        
ax.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)
ax2.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)
ax3.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)
ax4.scatter(FULL_tDt_IO_OUT_RATIO, FULL_eGFP_IO_OUT_RATIO, color='black', alpha=0.1)


a = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_eGFP_IO_OUT_RATIO)) if FULL_PROMOTER_NAMES[k]=='Igsf9']
b = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b-eGFP(G)']
c = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1-eGFP(G)']
d = [FULL_eGFP_IO_OUT_RATIO[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4-eGFP(G)']
ax5.boxplot([a, b, c, d], labels=['Igsf9','5HTr2b','PDX1','SUSD4'])


a = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='Igsf9']
b = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b-eGFP(G)']
c = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1-eGFP(G)']
d = [FULL_eGFP_IO_IN[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4-eGFP(G)']
ax6.boxplot([a, b, c, d], labels=['Igsf9','5HTr2b','PDX1','SUSD4'])

a = [FULL_eGFP_IO_OUT[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='Igsf9']
b = [FULL_eGFP_IO_OUT[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='5HTr2b-eGFP(G)']
c = [FULL_eGFP_IO_OUT[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='PDX1-eGFP(G)']
d = [FULL_eGFP_IO_OUT[k] for k in range(len(FULL_PROMOTER_NAMES)) if FULL_PROMOTER_NAMES[k]=='SUSD4-eGFP(G)']
ax7.boxplot([a, b, c, d], labels=['Igsf9','5HTr2b','PDX1','SUSD4'])

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
ax2.set_title('5HTR2b-eGFP')
ax3.set_title('PDX1-eGFP')
ax4.set_title('SUSD4-eGFP')
ax5.set_title('all')
plt.tight_layout()


#CONCENTRATION PLOT
plt.figure()
ax = plt.subplot(121)
ax2 = plt.subplot(122)
promoter = '5HT'
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
ax.set_title(promoter)