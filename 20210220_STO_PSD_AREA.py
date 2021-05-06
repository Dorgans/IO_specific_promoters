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
PIXEL_INTENSITY_NON_STO_FULL_ = []
PIXEL_INTENSITY_STO_FULL_ = []
LED_INTENSITY_ = []
PROMOTER_NAME_ = []
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
    if ('tif' in PATHS[i])==True and ('30fps' in PATHS[i])==True and ('IRED' in PATHS[i])==False:
        PATH = PATHS[i]
        try:
            
            
            IM_ARRAY_, w, h,  AVERAGE_IMAGE, SAMPLING_FREQUENCY, image_peak_coordinates = preprocess_image(PATHS[i], denoise_method__)
            #IM, AVERAGE_FreqDec_Pixels, CLUST_CENTERS, CLUST, p = Image_Decomposition(np.array(IM_ARRAY_), w, h, AVERAGE_IMAGE,5, 12, SAMPLING_FREQUENCY, IMG_NUM, 10, bin_size, method__, denoise_method__, clustering_method_, residuals__, options__    )
            IM = get_frequency_power_from_pixels(IM_ARRAY_, 4, 12, w, h, SAMPLING_FREQUENCY)

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
            
            
            DIVIDED = np.divide(np.log10(STO_PSD_FULL), np.log10(FluoIntensity_FULL))
            ZSCORE_DIVIDED = sp.stats.zscore(DIVIDED)
            ZSCORE_INTENSITY = sp.stats.zscore(np.log10(FluoIntensity_FULL))
            
            DIVIDED_THRESHOLDED = []
            PSD_THRESHOLDED = []
            THRESHOLDED_ARRAY_ABOVE = []
            THRESHOLDED_ARRAY_BELOW = []
            ZSCORE_THRESHOLDED = []
            FluoIntensity_NON_STO = []
            FluoIntensity_STO = []
            
            
            for j in range(len(DIVIDED)):
                if ZSCORE_DIVIDED[j]>2:
                    DIVIDED_THRESHOLDED.append(DIVIDED[j])
                    ZSCORE_THRESHOLDED.append(ZSCORE_DIVIDED[j])
                    THRESHOLDED_ARRAY_ABOVE.append(IM_ARRAY_[j])
                    PSD_THRESHOLDED.append(STO_PSD_FULL[j])
                    FluoIntensity_STO.append(FluoIntensity_FULL[j])
                else:
                    DIVIDED_THRESHOLDED.append(np.nan)
                    THRESHOLDED_ARRAY_BELOW.append(IM_ARRAY_[j])
                    FluoIntensity_NON_STO.append(FluoIntensity_FULL[j])
                    ZSCORE_THRESHOLDED.append(np.nan)
                    
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
            
            MEAN = np.nanmean(THRESHOLDED_ARRAY_BELOW, axis=0)
            SEM = sp.stats.sem(THRESHOLDED_ARRAY_BELOW, axis=0)
            ax5.plot(np.linspace(0, 10, len(MEAN)), MEAN, color='black', lw=0.5)
            ax5.fill_between(np.linspace(0, 10, len(MEAN)), MEAN+SEM, MEAN-SEM, color='black', alpha=0.1)
            MED_FILT = sp.signal.medfilt(MEAN, np.int(0.5*SAMPLING_FREQUENCY))
            df = pd.DataFrame()
            df['1'] = MEAN-MED_FILT
            df['2'] = MEAN-MED_FILT
            rs = [crosscorr(df['1'], df['1'], lag) for lag in range(-int(SAMPLING_FREQUENCY),int(SAMPLING_FREQUENCY+1))]
            ax6.plot(np.linspace(-1, 1, len(rs)), rs, color='black', lw=0.5)
            f, Pxx_spec = signal.welch(rs, SAMPLING_FREQUENCY , window = sp.signal.hann(len(rs), False), nfft = SAMPLING_FREQUENCY*200)
            ax8.plot(f, Pxx_spec, color='black', lw=0.5)
            DOMINANT_FREQUENCIES_BELOW_freqs_FULL_.append(f)
            DOMINANT_FREQUENCIES_BELOW_FULL_.append(Pxx_spec)

            try:
                MEAN = np.nanmean(THRESHOLDED_ARRAY_ABOVE, axis=0)
                SEM = sp.stats.sem(THRESHOLDED_ARRAY_ABOVE, axis=0)
                ax5.plot(np.linspace(0, 10, len(MEAN)), MEAN, color='orange', lw=0.5)
                ax5.fill_between(np.linspace(0, 10, len(MEAN)), MEAN+SEM, MEAN-SEM, color='orange', alpha=0.1)
                MED_FILT = sp.signal.medfilt(MEAN, np.int(0.5*SAMPLING_FREQUENCY))
                df = pd.DataFrame()
                df['1'] = MEAN-MED_FILT
                df['2'] = MEAN-MED_FILT
                rs = [crosscorr(df['1'], df['1'], lag) for lag in range(-int(SAMPLING_FREQUENCY),int(SAMPLING_FREQUENCY+1))]
                ax6.plot(np.linspace(-1, 1, len(rs)),rs, color='orange', lw=0.5)
                f, Pxx_spec = signal.welch(rs, SAMPLING_FREQUENCY , window = sp.signal.hann(len(rs), False), nfft = SAMPLING_FREQUENCY*200)
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
            PIXEL_PSD_INTENSITY_RATIO_.append(np.nanmean(DIVIDED_THRESHOLDED))
            PIXEL_INTENSITY_FULL_.append(np.nanmean(FluoIntensity_FULL))
            PIXEL_INTENSITY_NON_STO_FULL_.append(np.nanmean(FluoIntensity_NON_STO))
            PIXEL_INTENSITY_STO_FULL_.append(np.nanmean(FluoIntensity_STO))
            PIXEL_PSD_.append(PSD_THRESHOLDED)
            ACTIVE_PIXELS_.append(np.count_nonzero(np.nan_to_num(DIVIDED_THRESHOLDED)))
            NOISE_THRESHOLD_.append(Z_THRESHOLD)
            #OPTICS_CLUSTERS_FULL_.append(AVERAGE_FreqDec_Pixels)
            IM_ARRAY_FULL_.append(FluoIntensity_FULL)
            IM_DEC_FULL.append(STO_PSD_FULL)
        except:
            pass
        

LIST_ = [ 'AAV.PHP.S-5HTr2b(3.7)-tTA_AAV.PHP.S-TRE-GCamp6s',
         'AAV9-5HTr2b(3.7)-tTA_AAV9-TRE-GCamp6s',
         'AAV.PHP.S-5HTr2b(3.7)-tTA_AAV.PHP.S-TRE-GCamp6s_ROIVInj',
         'AAV.PHP.S-5HTr2b(3.7)-tTA_AAV.PHP.S-TRE-GCamp6s_IVInj',
         'AAV9-SUSD4(2.4)-GCamp6s',
          'AAV9-Igsf9(2.5)-GCamp6s',
         'AAV9-5HTr2b(3.0)-GCamp6s']

         
         
plt.figure()
ax1 = plt.subplot(231)
ax2 = plt.subplot(232)
ax3 = plt.subplot(233)
ax4 = plt.subplot(234)
ax5 = plt.subplot(235)
ax6 = plt.subplot(236)

ACTIVE_PIXEL_PARAM_PER_PROMOTER = []
PSD_PARAM_PER_PROMOTER = []

for PromNum in LIST_:
    X = [LED_INTENSITY_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    CellActivityPerProm = [np.nanmean(PIXEL_PSD_[j])/np.nanmean(IM_DEC_FULL[j]) for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    PixelPowerPerProm = [np.nanmean(PIXEL_PSD_[j]) for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    ActivePixPerProm = [np.nanmean(ACTIVE_PIXELS_[j]) for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    SignalThresholdPerProm = [NOISE_THRESHOLD_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    AverageFluo = [PIXEL_INTENSITY_FULL_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]
    AverageFLuoCells = [PIXEL_INTENSITY_STO_FULL_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum and np.isnan(LED_INTENSITY_[j])==False]

    PSD_PARAM_PER_PROMOTER.append([CellActivityPerProm[j] for j in range(len(CellActivityPerProm)) if X[j]==80])
    ACTIVE_PIXEL_PARAM_PER_PROMOTER.append([AverageFLuoCells[j] for j in range(len(AverageFLuoCells)) if X[j]==80])
    
    Displayed_ = []
    Displayed_SEM_ = []
    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([CellActivityPerProm[j] for j in range(len(CellActivityPerProm)) if np.int(X[j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([CellActivityPerProm[j] for j in range(len(CellActivityPerProm)) if np.int(X[j])==np.int(k)]))
    Displayed_ = np.array(Displayed_)
    Displayed_SEM_ = np.array(Displayed_SEM_)
    ax1.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax1.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    print(PromNum+str(np.nanmean(SignalThresholdPerProm )))
    ax1.text(100, Displayed_[-1], PromNum)
    
    Displayed_ = []
    Displayed_SEM_ = []
    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([PixelPowerPerProm[j] for j in range(len(PixelPowerPerProm)) if np.int(X[j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([PixelPowerPerProm[j] for j in range(len(PixelPowerPerProm)) if np.int(X[j])==np.int(k)]))
    Displayed_ = np.array(Displayed_)
    Displayed_SEM_ = np.array(Displayed_SEM_)
    ax2.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax2.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    print(PromNum+str(np.nanmean(SignalThresholdPerProm )))
    ax2.text(100, Displayed_[-1], PromNum)
    
    Displayed_ = []
    Displayed_SEM_ = []
    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([ActivePixPerProm[j] for j in range(len(ActivePixPerProm)) if np.int(X[j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([ActivePixPerProm[j] for j in range(len(ActivePixPerProm)) if np.int(X[j])==np.int(k)]))
    Displayed_ = np.array(Displayed_)
    Displayed_SEM_ = np.array(Displayed_SEM_)
    ax3.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax3.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    print(PromNum+str(np.nanmean(SignalThresholdPerProm )))
    ax3.text(100, Displayed_[-1], PromNum)
    
    Displayed_ = []
    Displayed_SEM_ = []
    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([SignalThresholdPerProm[j] for j in range(len(SignalThresholdPerProm)) if np.int(X[j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([SignalThresholdPerProm[j] for j in range(len(SignalThresholdPerProm)) if np.int(X[j])==np.int(k)]))
    Displayed_ = np.array(Displayed_)
    Displayed_SEM_ = np.array(Displayed_SEM_)
    ax4.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax4.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    print(PromNum+str(np.nanmean(SignalThresholdPerProm )))
    ax4.text(100, Displayed_[-1], PromNum)

    Displayed_ = []
    Displayed_SEM_ = []
    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([AverageFluo[j] for j in range(len(AverageFluo)) if np.int(X[j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([AverageFluo[j] for j in range(len(AverageFluo)) if np.int(X[j])==np.int(k)]))
    Displayed_ = np.array(Displayed_)
    Displayed_SEM_ = np.array(Displayed_SEM_)
    ax5.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax5.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    print(PromNum+str(np.nanmean(SignalThresholdPerProm )))
    ax5.text(100, Displayed_[-1], PromNum)
    
    Displayed_ = []
    Displayed_SEM_ = []
    for k in [0, 10, 40, 80, 100]:
        Displayed_.append(np.nanmean([AverageFLuoCells[j] for j in range(len(AverageFLuoCells)) if np.int(X[j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([AverageFLuoCells[j] for j in range(len(AverageFLuoCells)) if np.int(X[j])==np.int(k)]))
    Displayed_ = np.array(Displayed_)
    Displayed_SEM_ = np.array(Displayed_SEM_)
    ax6.plot([0, 10, 40, 80, 100],np.array( Displayed_))
    ax6.fill_between([0, 10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    print(PromNum+str(np.nanmean(SignalThresholdPerProm )))
    ax6.text(100, Displayed_[-1], PromNum)
    
     
    
Z = [NOISE_THRESHOLD_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]=='SUSD4(2.4)-GCamp6s']
#plt.plot(X, np.array(Z))
ax1.set_title('CellActivityPerProm')
ax2.set_title('PixelPowerPerProm')
ax3.set_title('ActivePixPerPro ')
ax4.set_title('SignalThresholdPerProm')
ax5.set_title('AverageIntensity')
ax6.set_title('AverageIntensity(activePix)')

plt.figure()
ax = plt.subplot(311)
ax2 = plt.subplot(312)
ax3 = plt.subplot(313)
ax.boxplot(ACTIVE_PIXEL_PARAM_PER_PROMOTER, labels=LIST_, vert=False)
ax2.boxplot(PSD_PARAM_PER_PROMOTER, labels=LIST_, vert=False)


ax3.scatter([LIST_[1] for i in range(len(PSD_PARAM_PER_PROMOTER[1]))], PSD_PARAM_PER_PROMOTER[1], alpha=0.2)
ax3.scatter([LIST_[2] for i in range(len(PSD_PARAM_PER_PROMOTER[2]))], PSD_PARAM_PER_PROMOTER[2], alpha=0.2)
ax3.scatter([LIST_[3] for i in range(len(PSD_PARAM_PER_PROMOTER[3]))], PSD_PARAM_PER_PROMOTER[3], alpha=0.2)
ax3.scatter([LIST_[4] for i in range(len(PSD_PARAM_PER_PROMOTER[4]))], PSD_PARAM_PER_PROMOTER[4], alpha=0.2)

