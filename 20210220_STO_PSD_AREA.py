# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 12:13:58 2021

@author: KEVIN-DORGANS
"""
LIST_ = ['AAV9-5HTr2b(1.8)-GCamp6s',
 'AAV9-5HTr2b(3.0)-GCamp6s',
 'AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s',
 'AAV.PHP.S-5HTr2b(3.7)-tTA_AAV.PHP.S-TRE-GCamp6s',
 'AAV9-CAG-GCamp6s',
 'AAV9-Igsf9(1.3)-GCamp6s',
 'AAV9-Igsf9(2.5)-GCamp6s',
 'AAV9-SUSD4(2.4)-GCamp6s',]


options__= ['calculate intensity epicenter',
 'exclude low intensity pixels']

ACTIVE_PIXELS_ = []
ACTIVE_ROI_ = []
PIXEL_PSD_SIGNAL_FULL = []
PIXEL_PSD_BACKGROUND_FULL = []
PIXEL_PSD_ = []
PIXEL_PSD_INTENSITY_RATIO_ = []
PIXEL_INTENSITY_FULL_ = []
PIXEL_INTENSITY_TAIL_ = []
PIXEL_INTENSITY_NON_STO_FULL_ = []
PIXEL_INTENSITY_STO_FULL_ = []
LED_INTENSITY_ = []
PROMOTER_NAME_ = []
FILE_NAME = []
NOISE_THRESHOLD_ = []
SIGNAL_ACORR_ = []
DOMINANT_FREQUENCIES_ABOVE_FULL_ = []
DOMINANT_FREQUENCIES_BELOW_FULL_ = []
DOMINANT_FREQUENCIES_ABOVE_freqs_FULL_ = []
DOMINANT_FREQUENCIES_BELOW_freqs_FULL_ = []
OPTICS_CLUSTERS_FULL_ = []

IM_ARRAY_FULL_ = []
IM_DEC_FULL = []

PATHS = load_directory_content_and_sub__()[0]
for i in range(len(PATHS)):
    if ('tif' in PATHS[i])==True and (('60fps' in PATHS[i])==True or ('30fps' in PATHS[i])==True)==True and ('IRED' in PATHS[i])==False:
        PATH = PATHS[i]
        try:
            
            
            IM_ARRAY_, w, h,  AVERAGE_IMAGE, SAMPLING_FREQUENCY, image_peak_coordinates = preprocess_image(PATHS[i], denoise_method__)
            #IM, AVERAGE_FreqDec_Pixels, CLUST_CENTERS, CLUST, p = Image_Decomposition(np.array(IM_ARRAY_), w, h, AVERAGE_IMAGE,5, 12, SAMPLING_FREQUENCY, IMG_NUM, 10, bin_size, method__, denoise_method__, clustering_method_, residuals__, options__    )
            IM = get_frequency_power_from_pixels(IM_ARRAY_, 3.7, 10, w, h, SAMPLING_FREQUENCY)

            #clr = cm.tab20c(np.linspace(0, 1, np.nanmax(CLUST)+1)).tolist()
            
            plt.figure(figsize=(9,6), num=PATH)
            ax1 = plt.subplot(241)
            ax2 = plt.subplot(242, sharex=ax1, sharey=ax1)
            ax3 = plt.subplot(243, sharex=ax1, sharey=ax1)
            ax4 = plt.subplot(244, sharex=ax1, sharey=ax1)
            ax5 = plt.subplot(245)
            ax6 = plt.subplot(246)
            ax7 = plt.subplot(247, sharex=ax1, sharey=ax1)
            ax8 = plt.subplot(248)
            
            ax1.imshow(AVERAGE_IMAGE)
            ax1.set_title('Fluo.Intensity')
            ax2.imshow(IM)
            ax2.set_title('PSD IM.')
            
            FluoIntensity_FULL = np.concatenate(AVERAGE_IMAGE)
            STO_PSD_FULL = STO_PSD = np.concatenate(IM)
            
            
            DIVIDED = np.divide(np.log10(STO_PSD_FULL) , np.log10(FluoIntensity_FULL))
            ZSCORE_DIVIDED = sp.stats.zscore(DIVIDED)
            ZSCORE_INTENSITY = sp.stats.zscore(np.log10(FluoIntensity_FULL))
            
            DIVIDED_THRESHOLDED = []
            PSD_THRESHOLDED = []
            PSD_NOISE = []
            THRESHOLDED_ARRAY_ABOVE = []
            THRESHOLDED_ARRAY_BELOW = []
            ZSCORE_THRESHOLDED = []
            FluoIntensity_NON_STO = []
            FluoIntensity_STO = []
            
            
            for j in range(len(DIVIDED)):
                if ZSCORE_DIVIDED[j]>3:
                    DIVIDED_THRESHOLDED.append(DIVIDED[j])
                    ZSCORE_THRESHOLDED.append(ZSCORE_DIVIDED[j])
                    THRESHOLDED_ARRAY_ABOVE.append(IM_ARRAY_[j])
                    PSD_THRESHOLDED.append(STO_PSD_FULL[j])
                    FluoIntensity_STO.append(FluoIntensity_FULL[j])                    
                else:
                    DIVIDED_THRESHOLDED.append(np.nan)
                    ZSCORE_THRESHOLDED.append(np.nan)
                    
                if True:
                    THRESHOLDED_ARRAY_BELOW.append(IM_ARRAY_[j])
                    FluoIntensity_NON_STO.append(FluoIntensity_FULL[j])
                    PSD_NOISE.append(STO_PSD_FULL[j])




            Z_THRESHOLD = np.nanmin(DIVIDED_THRESHOLDED)
            
            if len(DIVIDED_THRESHOLDED)==0:
                #RandomPixels = [np.random.randint(len(DIVIDED)) for l in range(10)]
                #DIVIDED_THRESHOLDED = [DIVIDED[RandomPixels[l]] for l in range(len(RandomPixels))]
                #THRESHOLDED_ARRAY = [IM_ARRAY_[RandomPixels[l]] for l in range(len(RandomPixels))]
                #PSD_THRESHOLDED= [STO_PSD_FULL[RandomPixels[l]] for l in range(len(RandomPixels))]
                ACTIVE_ROI_.append(0)
            else:
                ACTIVE_ROI_.append(1)
                
            ax3.imshow(np.array(DIVIDED).reshape(-1, w))
            ax3.set_title('Fluo.I/PSD.IM')
            ax4.imshow(np.array(DIVIDED_THRESHOLDED).reshape(-1, w))
            ax4.set_title('Fluo.I/PSD.IM (Thresh.)')
            ax7.imshow(np.array(ZSCORE_THRESHOLDED).reshape(-1, w))
            ax7.set_title('Fluo.I/PSD.IM (Z-SCORE)')
            
            SCALED_TO_MIN = [THRESHOLDED_ARRAY_BELOW[i]-np.nanmin(THRESHOLDED_ARRAY_BELOW[i]) for i in range(len(THRESHOLDED_ARRAY_BELOW))]
            MEAN = np.nanmean(SCALED_TO_MIN, axis=0)
            SEM = sp.stats.sem(SCALED_TO_MIN, axis=0, nan_policy='omit')
            ax5.plot(np.linspace(0, 10, len(MEAN)), MEAN, color='black', lw=0.5)
            ax5.fill_between(np.linspace(0, 10, len(MEAN)), MEAN+SEM, MEAN-SEM, color='black', alpha=0.1)

            try:
                MED_FILT = sp.signal.medfilt(MEAN, np.int(0.5*SAMPLING_FREQUENCY))
            except:
                MED_FILT = sp.signal.medfilt(MEAN, np.int(0.5*SAMPLING_FREQUENCY)-1)
            
            df = pd.DataFrame()
            df['1'] = MEAN - MED_FILT
            df['2'] = MEAN - MED_FILT
            rs = [crosscorr(df['1'], df['1'], lag) for lag in range(-int(SAMPLING_FREQUENCY),int(SAMPLING_FREQUENCY+1))]
            MEAN = MEAN - MED_FILT
            
            ax6.plot(np.linspace(-1, 1, len(rs)), rs, color='black', lw=0.5)
            f, Pxx_spec = signal.welch(MEAN , SAMPLING_FREQUENCY , window='blackman', noverlap = SAMPLING_FREQUENCY/1.5 , nperseg=SAMPLING_FREQUENCY*5 , nfft=SAMPLING_FREQUENCY*10, average='mean')
            ax8.plot(f, Pxx_spec, color='black', lw=0.5)
            DOMINANT_FREQUENCIES_BELOW_freqs_FULL_.append(f)
            DOMINANT_FREQUENCIES_BELOW_FULL_.append(Pxx_spec)

            try:
                SCALED_TO_MIN = [THRESHOLDED_ARRAY_ABOVE[i]-np.nanmin(THRESHOLDED_ARRAY_ABOVE[i]) for i in range(len(THRESHOLDED_ARRAY_ABOVE))]
                MEAN = np.nanmean(SCALED_TO_MIN, axis=0)
                SEM = sp.stats.sem(SCALED_TO_MIN, axis=0)
                ax5.plot(np.linspace(0, 10, len(MEAN)), MEAN, color='orange', lw=0.5)
                ax5.fill_between(np.linspace(0, 10, len(MEAN)), MEAN+SEM, MEAN-SEM, color='orange', alpha=0.1)
                
                """
                try:
                    MED_FILT = sp.signal.medfilt(MEAN, np.int(0.5*SAMPLING_FREQUENCY))
                except:
                    MED_FILT = sp.signal.medfilt(MEAN, np.int(0.5*SAMPLING_FREQUENCY)-1)
                """
                
                df = pd.DataFrame()
                df['1'] = MEAN# - MED_FILT
                df['2'] = MEAN# - MED_FILT
                rs = [crosscorr(df['1'], df['1'], lag) for lag in range(-int(SAMPLING_FREQUENCY),int(SAMPLING_FREQUENCY+1))]
                MEAN = MEAN - MED_FILT
                
                ax6.plot(np.linspace(-1, 1, len(rs)),rs, color='orange', lw=0.5)
                f, Pxx_spec = signal.welch(MEAN , SAMPLING_FREQUENCY , window='blackman', noverlap = SAMPLING_FREQUENCY/0.5 , nperseg=SAMPLING_FREQUENCY*20 , nfft=SAMPLING_FREQUENCY*20, average='mean', scaling='spectrum')
                ax8.plot(f, Pxx_spec, color='orange', lw=0.5)
                DOMINANT_FREQUENCIES_ABOVE_FULL_.append(Pxx_spec)
                DOMINANT_FREQUENCIES_ABOVE_freqs_FULL_.append(f)
            except:
                DOMINANT_FREQUENCIES_ABOVE_FULL_.append(np.nan)
                DOMINANT_FREQUENCIES_ABOVE_freqs_FULL_.append(np.nan)
                
            SIGNAL_ACORR_.append(rs)
            try:
                LED_INTENSITY_.append(np.int(PATH.split('LED')[-1].split('_')[0]))
            except:
                LED_INTENSITY_.append(np.nan)
            
                
            PROMOTER_NAME_.append(PATH.split('EXPORT\\')[1].split('\\')[0])
            FILE_NAME.append(PATH)
            PIXEL_PSD_INTENSITY_RATIO_.append(np.nanmean(DIVIDED_THRESHOLDED))
            PIXEL_INTENSITY_FULL_.append(np.nanmean(np.nanmean(IM_ARRAY_, axis=0)[0:int(SAMPLING_FREQUENCY)]))
            PIXEL_INTENSITY_TAIL_.append(np.nanmean(np.nanmean(IM_ARRAY_, axis=0)[-int(SAMPLING_FREQUENCY):]))
            PIXEL_INTENSITY_NON_STO_FULL_.append(np.nanmean(FluoIntensity_NON_STO[0:int(SAMPLING_FREQUENCY)]))
            PIXEL_INTENSITY_STO_FULL_.append(np.nanmean(FluoIntensity_STO[0:int(SAMPLING_FREQUENCY)]))
            PIXEL_PSD_.append(PSD_THRESHOLDED)
            PIXEL_PSD_BACKGROUND_FULL.append(PSD_NOISE)
            ACTIVE_PIXELS_.append(np.count_nonzero(np.nan_to_num(DIVIDED_THRESHOLDED)))
            NOISE_THRESHOLD_.append(Z_THRESHOLD)
            #OPTICS_CLUSTERS_FULL_.append(AVERAGE_FreqDec_Pixels)
            IM_ARRAY_FULL_.append(FluoIntensity_FULL)
            IM_DEC_FULL.append(STO_PSD_FULL)
        except:
            pass
        
        
        


LIST_ = [ 'WHITE_NOISE', 
         'AAV.PHP.S-5HTr2b(3.7)-tTA_AAV.PHP.S-TRE-GCamp6s',
         'AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s',
         'AAV.PHP.eB-5HTr2b(3.7)-tTA_AAV.PHP.eB-TRE-GCamp6s']

LIST_ = [
 'AAV9-5HTr2b(1.8)-GCamp6s',
 'AAV9-5HTr2b(3.0)-GCamp6s',
 'AAV9-Igsf9(1.3)-GCamp6s',
 'AAV9-Igsf9(2.5)-GCamp6s',
 'AAV9-SUSD4(2.4)-GCamp6s',
 'AAV.PHP.S-5HTr2b(3.7)-tTA_AAV.PHP.S-TRE-GCamp6s',
 'AAV.PHP.eB-5HTr2b(3.7)-tTA_AAV.PHP.eB-TRE-GCamp6s',
  'AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s',
 'WHITE_NOISE']


LIST_ = [
 'AAV9-5HTr2b(3.0)-GCamp6s',
 'AAV9-5HTr2b(1.8)-GCamp6s',
 'AAV9-CAG-GCamp6s',
 'AAV9-Igsf9(1.3)-GCamp6s',
 'AAV9-Igsf9(2.5)-GCamp6s',
 'AAV9-SUSD4(2.4)-GCamp6s',
 'AAV.PHP.S-5HTr2b(3.7)-tTA_AAV.PHP.S-TRE-GCamp6s',
 'AAV.PHP.S-5HTr2b(3.7)-tTA_AAV.PHP.S-TRE-GCamp6s_IVInj',
 'AAV.PHP.S-5HTr2b(3.7)-tTA_AAV.PHP.S-TRE-GCamp6s_ROIVInj',
 'AAV.PHP.eB-5HTr2b(3.7)-tTA_AAV.PHP.eB-TRE-GCamp6s',
 'AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-ASAP3kv',
 'AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s_PATCH',
 'AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s_6MA',
  'AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s',
 'AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s_InVivo',
 'WHITE_NOISE']

LIST_ = [ 'AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s','AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s_PATCH','AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s_6MA']

LIST_ = ['AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s',
 'AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s_InVivo']
       
LIST_ = ['WHITE_NOISE',
 'AAV9-5HTr2b(1.8)-GCamp6s',
 'AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s',
 'AAV9-Igsf9(1.3)-GCamp6s',
 'AAV9-Igsf9(2.5)-GCamp6s',
 'AAV9-SUSD4(2.4)-GCamp6s']

LIST_ = ['AAV.PHP.S-5HTr2b(3.7)-tTA_AAV.PHP.S-TRE-GCamp6s',
 'AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s']
       
LIST_ = ['AAV9-5HTr2b(1.8)-GCamp6s',
 'AAV9-5HTr2b(3.0)-GCamp6s',
 'AAV9-Igsf9(1.3)-GCamp6s',
 'AAV9-Igsf9(2.5)-GCamp6s',
 'AAV9-SUSD4(2.4)-GCamp6s',
 'WHITE_NOISE','AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s']

plt.figure()
ax1 = plt.subplot(241)
ax2 = plt.subplot(242)
ax3 = plt.subplot(243)
ax4 = plt.subplot(244)
ax5 = plt.subplot(245)
ax6 = plt.subplot(246)
ax7 = plt.subplot(247)
ax8 = plt.subplot(248)

ACTIVE_PIXEL_PARAM_PER_PROMOTER = []
PSD_PARAM_PER_PROMOTER = []
AverageFluo_PER_PROMOTER = []
PixelPowerPerProm_PER_PROMOTER = []
AverageFluorescenceBleach_PER_PROMOTER = []
AverageFluoCells_PER_PROMOTER = []
AverageSTOPSDfromSpec_PER_PROMOTER = []
AverageSTOPSDfromSpec_Frequency_PER_PROMOTER = []
AverageSPIKEPSDfromSpec_PER_PROMOTER = []
AverageSPIKEPSDfromSpec_Frequency_PER_PROMOTER = []
LED_INTENSITY_PER_PROMOTER = []
FILE_NAME_PER_PROM = []
SAMPLING_FREQUENCY_PER_PROM =  []

for PromNum in LIST_:
    X = [LED_INTENSITY_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    LED_INTENSITY_PER_PROMOTER.append(X)
    CellActivityPerProm = [np.nanmean(PIXEL_PSD_[j])/np.nanmean(PIXEL_PSD_BACKGROUND_FULL[j]) for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    PixelPowerPerProm = [np.nanmean(PIXEL_PSD_[j]) for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    ActivePixPerProm = [np.nanmean(ACTIVE_PIXELS_[j]) for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    SignalBleachPerProm = [PIXEL_INTENSITY_TAIL_[j] /PIXEL_INTENSITY_FULL_[j]  for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    AverageFluo = [PIXEL_INTENSITY_FULL_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    AverageFLuoCells = [PIXEL_INTENSITY_STO_FULL_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    AverageFluorescenceActivityNormalized = [PIXEL_PSD_INTENSITY_RATIO_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    FileNamePerProm =  [FILE_NAME[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    SamplingFreqPerProm = [np.int(FILE_NAME[j].split('fps')[0].split('_')[-1]) for j in range(len(FILE_NAME)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]

    AverageSTOPSDfromSpec = []
    AverageSTOPSDfromSpec_Frequency = []
    AverageSPIKEPSDfromSpec = []
    AverageSPIKEPSDfromSpec_Frequency = []
                    
    for j in range(len(LED_INTENSITY_)):
        if PROMOTER_NAME_[j]==PromNum and LED_INTENSITY_[j]>=40:
            SF = np.int(FILE_NAME[j].split('fps')[0].split('_')[-1])
            #Isolation of putative STO signal band-  welch method
            try:
                AVG_PSD_FREQBAND = np.nanmax([DOMINANT_FREQUENCIES_BELOW_FULL_[j][k] for k in range(len(DOMINANT_FREQUENCIES_BELOW_FULL_[j])) if 15>DOMINANT_FREQUENCIES_BELOW_freqs_FULL_[j][k]>1.5])
                DOMINANT_STO_FREQ = DOMINANT_FREQUENCIES_BELOW_freqs_FULL_[j][np.array(DOMINANT_FREQUENCIES_BELOW_FULL_[j]).tolist().index(AVG_PSD_FREQBAND)]
                if DOMINANT_STO_FREQ >2.5 and AVG_PSD_FREQBAND>0.1: #This is 3xSTD+AVG 1mW WHITE NOISE THRESHOLD
                    AverageSTOPSDfromSpec.append(AVG_PSD_FREQBAND*0.0014665418*SF) #bit to nW conversion ... fW.s/(1/SAMPLING in sec) Knowing Micam03WellDepth(e-), eV_to_W.s_conversion
                    AverageSTOPSDfromSpec_Frequency.append(DOMINANT_STO_FREQ )
                else:
                    AverageSTOPSDfromSpec_Frequency.append(np.nan)
                    AverageSTOPSDfromSpec.append(np.nan)
            except:
                AverageSTOPSDfromSpec.append(np.nan)
                AverageSTOPSDfromSpec_Frequency.append(np.nan)

            #Isolation of putative spike signal band - welch method
            try:
                AVG_PSD_FREQBAND = np.nanmax([DOMINANT_FREQUENCIES_ABOVE_FULL_[j][k] for k in range(len(DOMINANT_FREQUENCIES_ABOVE_FULL_[j])) if 2>DOMINANT_FREQUENCIES_ABOVE_freqs_FULL_[j][k]>0.3])
                DOMINANT_STO_FREQ = DOMINANT_FREQUENCIES_ABOVE_freqs_FULL_[j][np.array(DOMINANT_FREQUENCIES_ABOVE_FULL_[j]).tolist().index(AVG_PSD_FREQBAND)]
                if 2000>AVG_PSD_FREQBAND>10 : #This is 3xSTD+AVG 1mW WHITE NOISE THRESHOLD + above 800b is a aberrant value (eye checked on 100pct trials)
                    AVG_PSD_FREQBAND = np.nanmax([DOMINANT_FREQUENCIES_ABOVE_FULL_[j][k] for k in range(len(DOMINANT_FREQUENCIES_ABOVE_FULL_[j])) if 2>DOMINANT_FREQUENCIES_ABOVE_freqs_FULL_[j][k]>0.33])
                    AverageSPIKEPSDfromSpec.append(AVG_PSD_FREQBAND*0.0014665418*SF) 
                    AverageSPIKEPSDfromSpec_Frequency.append(DOMINANT_STO_FREQ )
                else:
                    AverageSPIKEPSDfromSpec.append(np.nan)
                    AverageSPIKEPSDfromSpec_Frequency.append(np.nan)  
            except:
                AverageSPIKEPSDfromSpec.append(np.nan)
                AverageSPIKEPSDfromSpec_Frequency.append(np.nan)

    PSD_PARAM_PER_PROMOTER.append([CellActivityPerProm[j] for j in range(len(CellActivityPerProm))])
    ACTIVE_PIXEL_PARAM_PER_PROMOTER.append([PixelPowerPerProm[j] for j in range(len(AverageFLuoCells))])
    AverageFluo_PER_PROMOTER.append([AverageFluo[j] for j in range(len(AverageFLuoCells))])
    PixelPowerPerProm_PER_PROMOTER.append([PixelPowerPerProm[j] for j in range(len(AverageFLuoCells))])
    AverageFluorescenceBleach_PER_PROMOTER.append([SignalBleachPerProm[j] for j in range(len(AverageFLuoCells))])
    AverageFluoCells_PER_PROMOTER.append([AverageFLuoCells[j] for j in range(len(AverageFLuoCells))])
    AverageSTOPSDfromSpec_PER_PROMOTER.append(AverageSTOPSDfromSpec)
    AverageSTOPSDfromSpec_Frequency_PER_PROMOTER.append(AverageSTOPSDfromSpec_Frequency)
    AverageSPIKEPSDfromSpec_PER_PROMOTER.append(AverageSPIKEPSDfromSpec)
    AverageSPIKEPSDfromSpec_Frequency_PER_PROMOTER.append(AverageSPIKEPSDfromSpec_Frequency)
    FILE_NAME_PER_PROM.append(FileNamePerProm)
    SAMPLING_FREQUENCY_PER_PROM.append(SamplingFreqPerProm)
    
for i in range(len(LIST_)):
    Displayed_ = []
    Displayed_SEM_ = []

    PromNum = LIST_[i]
    
    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([PSD_PARAM_PER_PROMOTER[i][j] for j in range(len(PSD_PARAM_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([PSD_PARAM_PER_PROMOTER[i][j] for j in range(len(PSD_PARAM_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)], nan_policy='omit'))

    Displayed_ = np.array(Displayed_) 
    Displayed_SEM_ = np.array(Displayed_SEM_)  
    ax1.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax1.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)

    ax1.text(100, Displayed_[-1], PromNum)
    
    Displayed_ = []
    Displayed_SEM_ = []

    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([PixelPowerPerProm_PER_PROMOTER[i][j]*SAMPLING_FREQUENCY_PER_PROM[i][j] for j in range(len(PixelPowerPerProm_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([PixelPowerPerProm_PER_PROMOTER[i][j]*SAMPLING_FREQUENCY_PER_PROM[i][j] for j in range(len(PixelPowerPerProm_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)], nan_policy='omit'))


    Displayed_ = np.array(Displayed_) *0.0014665418 #bit to nW conversion ... fW.s/(1/SAMPLING in sec) Knowing Micam03WellDepth(e-), eV_to_W.s_conversion
    Displayed_SEM_ = np.array(Displayed_SEM_) *0.0014665418 #bit to nW conversion ... fW.s/(1/SAMPLING in sec) Knowing Micam03WellDepth(e-), eV_to_W.s_conversion
    ax2.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax2.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    #print(PromNum+str(np.nanmean(SignalThresholdPerProm )))
    ax2.text(100, Displayed_[-1], PromNum)
    
    Displayed_ = []
    Displayed_SEM_ = []

    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([ACTIVE_PIXEL_PARAM_PER_PROMOTER[i][j] for j in range(len(ACTIVE_PIXEL_PARAM_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([ACTIVE_PIXEL_PARAM_PER_PROMOTER[i][j] for j in range(len(ACTIVE_PIXEL_PARAM_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)], nan_policy='omit'))

    Displayed_ = np.array(Displayed_) 
    Displayed_SEM_ = np.array(Displayed_SEM_) 
    ax3.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax3.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    #print(PromNum+str(np.nanmean(SignalThresholdPerProm )))
    ax3.text(100, Displayed_[-1], PromNum)
    
    Displayed_ = []
    Displayed_SEM_ = []

    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([AverageSTOPSDfromSpec_Frequency_PER_PROMOTER[i][j] for j in range(len(AverageSTOPSDfromSpec_Frequency_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([AverageSTOPSDfromSpec_Frequency_PER_PROMOTER[i][j] for j in range(len(AverageSTOPSDfromSpec_Frequency_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)], nan_policy='omit'))

    Displayed_ = np.array(Displayed_) 
    Displayed_SEM_ = np.array(Displayed_SEM_)
    ax4.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax4.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    #print(PromNum+str(np.nanmean(SignalThresholdPerProm )))
    ax4.text(100, Displayed_[-1], PromNum)

    Displayed_ = []
    Displayed_SEM_ = []

    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([AverageFluo_PER_PROMOTER[i][j]*SAMPLING_FREQUENCY_PER_PROM[i][j] for j in range(len(AverageFluo_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([AverageFluo_PER_PROMOTER[i][j]*SAMPLING_FREQUENCY_PER_PROM[i][j] for j in range(len(AverageFluo_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)], nan_policy='omit'))

    Displayed_ = np.array(Displayed_)*0.0014665418 #bit to nW conversion ... fW.s/(1/SAMPLING in sec) Knowing Micam03WellDepth(e-), eV_to_W.s_conversion
    Displayed_SEM_ = np.array(Displayed_SEM_)*0.0014665418 #bit to nW conversion ... fW.s/(1/SAMPLING in sec) Knowing Micam03WellDepth(e-), eV_to_W.s_conversion
    ax5.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax5.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    #print(PromNum+str(np.nanmean(SignalThresholdPerProm )))
    ax5.text(100, Displayed_[-1], PromNum)
    
    Displayed_ = []
    Displayed_SEM_ = []
    Displayed_ref = []
    Displayed_SEM_ref = []
    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([AverageFluoCells_PER_PROMOTER[i][j]*SAMPLING_FREQUENCY_PER_PROM[i][j]  for j in range(len(AverageFluoCells_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([AverageFluoCells_PER_PROMOTER[i][j]*SAMPLING_FREQUENCY_PER_PROM[i][j]  for j in range(len(AverageFluoCells_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)], nan_policy='omit'))
        Displayed_ref.append(np.nanmean([AverageFluoCells_PER_PROMOTER[0][j]  for j in range(len(AverageFluoCells_PER_PROMOTER[0])) if np.int(LED_INTENSITY_PER_PROMOTER[0][j])==np.int(k)]))
        Displayed_SEM_ref.append(sp.stats.sem([AverageFluoCells_PER_PROMOTER[0][j]  for j in range(len(AverageFluoCells_PER_PROMOTER[0])) if np.int(LED_INTENSITY_PER_PROMOTER[0][j])==np.int(k)], nan_policy='omit'))

    Displayed_ = np.array(Displayed_)*0.0014665418 #bit to nW conversion ... fW.s/(1/SAMPLING in sec) Knowing Micam03WellDepth(e-), eV_to_W.s_conversion
    Displayed_SEM_ = np.array(Displayed_SEM_)*0.0014665418 #bit to nW conversion ... fW.s/(1/SAMPLING in sec) Knowing Micam03WellDepth(e-), eV_to_W.s_conversion
    ax6.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax6.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    #print(PromNum+str(np.nanmean(SignalThresholdPerProm )))
    ax6.text(100, Displayed_[-1], PromNum)
    
    Displayed_ = []
    Displayed_SEM_ = []

    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([AverageFluorescenceBleach_PER_PROMOTER[i][j] for j in range(len(AverageFluorescenceBleach_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([AverageFluorescenceBleach_PER_PROMOTER[i][j] for j in range(len(AverageFluorescenceBleach_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)], nan_policy='omit'))

    Displayed_ = np.array(Displayed_) 
    Displayed_SEM_ = np.array(Displayed_SEM_) 
    ax7.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax7.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    #print(PromNum+str(np.nanmean(SignalThresholdPerProm )))
    ax7.text(100, Displayed_[-1], PromNum)

    Displayed_ = []
    Displayed_SEM_ = []

    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([AverageSTOPSDfromSpec_PER_PROMOTER[i][j] for j in range(len(AverageSTOPSDfromSpec_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([AverageSTOPSDfromSpec_PER_PROMOTER[i][j] for j in range(len(AverageSTOPSDfromSpec_PER_PROMOTER[i])) if np.int(LED_INTENSITY_PER_PROMOTER[i][j])==np.int(k)], nan_policy='omit'))

    Displayed_ = np.array(Displayed_) 
    Displayed_SEM_ = np.array(Displayed_SEM_)
    ax8.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax8.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    #print(PromNum+str(np.nanmean(SignalThresholdPerProm )))
    ax8.text(100, Displayed_[-1], PromNum)
      
    
Z = [NOISE_THRESHOLD_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]=='SUSD4(2.4)-GCamp6s']
#plt.plot(X, np.array(Z))
ax1.set_title('CellActivityPerProm')
ax2.set_title('PixelPowerPerProm')
ax3.set_title('ActivePixPerPro ')
ax4.set_title('STO dominant Freq')
ax5.set_title('AverageIntensity')
ax6.set_title('AverageIntensity(activePix)')
ax7.set_title('Average PSD_IM/F')
ax8.set_title('SpecFreqSTO-PSD')


plt.figure()
ax = plt.subplot(811)
ax2 = plt.subplot(812)
ax3 = plt.subplot(813)
ax4 = plt.subplot(814)
ax5 = plt.subplot(815)
ax6 = plt.subplot(816)
ax7 = plt.subplot(817)
ax8 = plt.subplot(818)

ax.boxplot(ACTIVE_PIXEL_PARAM_PER_PROMOTER, labels=LIST_, vert=False)
ax2.boxplot(PSD_PARAM_PER_PROMOTER, labels=LIST_, vert=False)
ax3.boxplot(AverageFluo_PER_PROMOTER, labels=LIST_, vert=False)
ax4.boxplot(PixelPowerPerProm_PER_PROMOTER, labels=LIST_, vert=False)
ax5.boxplot(AverageFluorescenceBleach_PER_PROMOTER, labels=LIST_, vert=False)
ax6.boxplot(AverageFluoCells_PER_PROMOTER, labels=LIST_, vert=False)
ax7.boxplot(AverageSTOPSDfromSpec_PER_PROMOTER, labels=LIST_, vert=False)
ax8.boxplot(AverageSTOPSDfromSpec_Frequency_PER_PROMOTER, labels=LIST_, vert=False)


ax.set_xlabel('Active pixel #')
ax2.set_xlabel('STO-PSD (ActivePixel/Average)')
ax3.set_xlabel('Avg fluo-intensity (FullImage)')
ax4.set_xlabel('STO-PSD (ActivePix)')
ax5.set_xlabel('PSD-IM/F (ActivePix)')
ax6.set_xlabel('Avg fluo-intensity (ActivePix)')
ax7.set_xlabel('Calcium Oscillation power (3-12Hz)')
ax8.set_xlabel('Calcium Oscillation Freq (3-12Hz)')
#for j in range(len(LIST_)):
    #ax3.scatter([LIST_[j] for i in range(len(PSD_PARAM_PER_PROMOTER[i]))], PSD_PARAM_PER_PROMOTER[i], alpha=0.2)

plt.figure(figsize=(3,3))
for i in range(len(LIST_)):
    idx  = np.where(np.array(LED_INTENSITY_PER_PROMOTER[i])>40)[0]
    X = np.nanmean([np.nanmean(ACTIVE_PIXEL_PARAM_PER_PROMOTER[i][j]) for j in range(len(LED_INTENSITY_PER_PROMOTER[i])) if (j in idx)==True])
    Y = np.nanmean([np.nanmean(PSD_PARAM_PER_PROMOTER[i][j]) for j in range(len(LED_INTENSITY_PER_PROMOTER[i])) if (j in idx)==True])
    plt.scatter(X, Y)
    plt.text(X, Y, LIST_[i])

plt.ylabel('STO-PSD SNR (ActivePixel/Average)')
plt.xlabel('F-RMS (4-12Hz)')
plt.tight_layout()

plt.figure(figsize=(3,3))
Freq_of_recorded_STO_IO_FIELD_InVitro = []
AVG_STO_IO_FIELD_InVitro = []
AVG_SPIKE_IO_FIELD_InVitro = []
AVG_Fluo_IO_FIELD_per_recorded_STO = []
AVG_Fluo_IO_FIELD_per_recorded_SPIKE = []

ALL_FLUO_NZ_OBS = []
ALL_SPIKE_NZ_OBS = []

for i in range(len(AverageFluoCells_PER_PROMOTER)):
    plt.scatter(AverageSTOPSDfromSpec_Frequency_PER_PROMOTER[i], AverageSTOPSDfromSpec_PER_PROMOTER[i], label=LIST_[i])
    NON_ZERO_OBSERVATIONS = [AverageSTOPSDfromSpec_PER_PROMOTER[i][j] for j in range(len(AverageSTOPSDfromSpec_PER_PROMOTER[i])) if AverageSTOPSDfromSpec_PER_PROMOTER[i][j]>0]
    NON_ZERO_OBSERVATIONS_Fluo = [ AverageFluoCells_PER_PROMOTER[i][j]*SAMPLING_FREQUENCY_PER_PROM[i][j]*0.0014665418  for j in range(len(AverageSTOPSDfromSpec_PER_PROMOTER[i])) if AverageSTOPSDfromSpec_PER_PROMOTER[i][j]>0]
    NON_ZERO_OBSERVATIONS_Spike = [AverageSPIKEPSDfromSpec_PER_PROMOTER[i][j] for j in range(len(AverageSPIKEPSDfromSpec_PER_PROMOTER[i])) if AverageSPIKEPSDfromSpec_PER_PROMOTER[i][j]>0]
    NON_ZERO_OBSERVATIONS_Spike_Fluo = [AverageFluoCells_PER_PROMOTER[i][j]*SAMPLING_FREQUENCY_PER_PROM[i][j]*0.0014665418 for j in range(len(AverageSPIKEPSDfromSpec_PER_PROMOTER[i])) if AverageSPIKEPSDfromSpec_PER_PROMOTER[i][j]>0]

    ALL_FLUO_NZ_OBS.append(NON_ZERO_OBSERVATIONS_Spike_Fluo)
    ALL_SPIKE_NZ_OBS.append(NON_ZERO_OBSERVATIONS_Spike)
    
    AVG_STO_IO_FIELD_InVitro.append(np.nanmean(NON_ZERO_OBSERVATIONS))
    AVG_SPIKE_IO_FIELD_InVitro.append(np.nanmean(NON_ZERO_OBSERVATIONS_Spike))
    AVG_Fluo_IO_FIELD_per_recorded_STO.append(np.nanmean(NON_ZERO_OBSERVATIONS_Fluo))
    AVG_Fluo_IO_FIELD_per_recorded_SPIKE.append(np.nanmean(NON_ZERO_OBSERVATIONS_Spike_Fluo))

    Freq_of_recorded_STO_IO_FIELD_InVitro.append(np.count_nonzero(np.nan_to_num(AverageSTOPSDfromSpec_PER_PROMOTER[i]))/len(AverageSTOPSDfromSpec_PER_PROMOTER[i]))
plt.legend(loc="upper left")

plt.figure(figsize=(3,3))
for i in range(len(LIST_)):
    plt.scatter(Freq_of_recorded_STO_IO_FIELD_InVitro[i], AVG_STO_IO_FIELD_InVitro[i], label=LIST_[i])
    plt.text(Freq_of_recorded_STO_IO_FIELD_InVitro[i], AVG_STO_IO_FIELD_InVitro[i],LIST_[i])
plt.show()