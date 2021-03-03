# -*- coding: utf-8 -*-
"""
Created on Sat Feb 20 12:13:58 2021

@author: KEVIN-DORGANS
"""

options__= ['calculate intensity epicenter',
 'exclude low intensity pixels']

ACTIVE_PIXELS_ = []
ACTIVE_ROI_ = []
PIXEL_PSD_ = []
PIXEL_PSD_INTENSITY_RATIO_ = []
LED_INTENSITY_ = []
PROMOTER_NAME_ = []
NOISE_THRESHOLD_ = []
SIGNAL_ACORR_ = []

PATHS = load_directory_content_and_sub__()[0]
for i in range(len(PATHS)):
    if ('tif' in PATHS[i])==True and ('30fps' in PATHS[i])==True and ('IRED' in PATHS[i])==False:
        PATH = PATHS[i]
        try:
            
            
            IM_ARRAY_, w, h,  AVERAGE_IMAGE, SAMPLING_FREQUENCY, image_peak_coordinates = preprocess_image(PATHS[i], denoise_method__)
            IM, AVERAGE_FreqDec_Pixels, CLUST_CENTERS, CLUST, p = Image_Decomposition(np.array(IM_ARRAY_), w, h, AVERAGE_IMAGE,5, 12, SAMPLING_FREQUENCY, IMG_NUM, 10, bin_size, method__, denoise_method__, clustering_method_, residuals__, options__    )
            clr = cm.tab20c(np.linspace(0, 1, np.nanmax(CLUST)+1)).tolist()
            
            plt.figure(figsize=(9,6))
            ax1 = plt.subplot(241)
            ax2 = plt.subplot(242, sharex=ax1, sharey=ax1)
            ax3 = plt.subplot(243, sharex=ax1, sharey=ax1)
            ax4 = plt.subplot(244, sharex=ax1, sharey=ax1)
            ax5 = plt.subplot(245)
            ax6 = plt.subplot(246)
            ax7 = plt.subplot(247, sharex=ax1, sharey=ax1)
            
            ax1.imshow(AVERAGE_IMAGE)
            ax2.imshow(IM)
            
            FluoIntensity_FULL = FluoIntensity = np.concatenate(AVERAGE_IMAGE)
            STO_PSD_FULL = STO_PSD = np.concatenate(IM)
            
            
            DIVIDED = np.divide(np.log10(STO_PSD_FULL), np.log10(FluoIntensity_FULL))
            ZSCORE_DIVIDED = sp.stats.zscore(DIVIDED)
            DIVIDED_THRESHOLDED = []
            PSD_THRESHOLDED = []
            THRESHOLDED_ARRAY = []
            ZSCORE_THRESHOLDED = []
            
            for j in range(len(DIVIDED)):
                if ZSCORE_DIVIDED[j]>0:
                    ZSCORE_THRESHOLDED.append(ZSCORE_DIVIDED[j])
                else:
                    ZSCORE_THRESHOLDED.append(np.nan)
            

            for j in range(len(DIVIDED)):
                if ZSCORE_DIVIDED[j]>3:
                    DIVIDED_THRESHOLDED.append(DIVIDED[j])
                    THRESHOLDED_ARRAY.append(IM_ARRAY_[j])
                    PSD_THRESHOLDED.append(STO_PSD_FULL[j])
                else:
                    DIVIDED_THRESHOLDED.append(np.nan)
                    
            Z_THRESHOLD = np.nanmin(DIVIDED_THRESHOLDED)
            
            if len(DIVIDED_THRESHOLDED)==0:
                RandomPixels = [np.random.randint(len(DIVIDED)) for l in range(10)]
                DIVIDED_THRESHOLDED = [DIVIDED[RandomPixels[l]] for l in range(len(RandomPixels))]
                THRESHOLDED_ARRAY = [IM_ARRAY_[RandomPixels[l]] for l in range(len(RandomPixels))]
                PSD_THRESHOLDED= [STO_PSD_FULL[RandomPixels[l]] for l in range(len(RandomPixels))]
                ACTIVE_ROI_.append(0)
            else:
                ACTIVE_ROI_.append(1)
                
            ax3.imshow(np.array(DIVIDED).reshape(-1, w))
            ax4.imshow(np.array(DIVIDED_THRESHOLDED).reshape(-1, w))
            ax7.imshow(np.array(ZSCORE_THRESHOLDED).reshape(-1, w))
            
            MEAN = np.nanmean(THRESHOLDED_ARRAY, axis=0)
            SEM = sp.stats.sem(THRESHOLDED_ARRAY, axis=0)
            ax5.plot(np.linspace(0, 10, len(MEAN)), MEAN, color='black')
            ax5.fill_between(np.linspace(0, 10, len(MEAN)), MEAN+SEM, MEAN-SEM, color='black', alpha=0.1)

            MED_FILT = sp.signal.medfilt(MEAN, np.int(0.5*SAMPLING_FREQUENCY))
            df = pd.DataFrame()
            df['1'] = MEAN-MED_FILT
            df['2'] = MEAN-MED_FILT
            rs = [crosscorr(df['1'], df['1'], lag) for lag in range(-int(SAMPLING_FREQUENCY),int(SAMPLING_FREQUENCY+1))]
            ax6.plot(rs)
            SIGNAL_ACORR_.append(rs)
            LED_INTENSITY_.append(np.int(PATH.split('LED')[-1].split('_')[0]))
            PROMOTER_NAME_.append(PATH.split('EXPORT\\')[1].split('\\')[0])
            PIXEL_PSD_INTENSITY_RATIO_.append(np.nanmean(DIVIDED_THRESHOLDED))
            PIXEL_PSD_.append(PSD_THRESHOLDED)
            ACTIVE_PIXELS_.append(np.count_nonzero(np.nan_to_num(DIVIDED_THRESHOLDED)))
            NOISE_THRESHOLD_.append(Z_THRESHOLD)
            
        except:
            pass
        
        
plt.figure()
ACTIVE_PIXEL_PARAM_PER_PROMOTER = []
PSD_PARAM_PER_PROMOTER = []

for PromNum in LIST_:
    X = [LED_INTENSITY_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum]
    Y = [np.nanmean(PIXEL_PSD_[j]) for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum]
    Z = [PIXEL_PSD_INTENSITY_RATIO_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum]
    ZX = [ACTIVE_PIXELS_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum]
    
    NOISE_THRS = [NOISE_THRESHOLD_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]==PromNum]
    Displayed_ = []
    Displayed_SEM_ = []
    for k in [10, 40, 80, 100]:
        Displayed_.append(np.nanmean([Y[j] for j in range(len(Y)) if np.int(X[j])==np.int(k)]))
        Displayed_SEM_.append(sp.stats.sem([Y[j] for j in range(len(Y)) if np.int(X[j])==np.int(k)]))
    PSD_PARAM_PER_PROMOTER.append([Y[j] for j in range(len(Y)) if np.int(X[j])>10])
    ACTIVE_PIXEL_PARAM_PER_PROMOTER.append([ZX[j] for j in range(len(Y)) if np.int(X[j])>10])
    
    
    Displayed_ = np.array(Displayed_)
    Displayed_SEM_ = np.array(Displayed_SEM_)
    plt.plot([10, 40, 80, 100],np.array( Displayed_))
    plt.fill_between([10, 40, 80, 100], Displayed_SEM_+Displayed_, Displayed_-Displayed_SEM_, alpha=0.1)
    print(PromNum+str(np.nanmean(NOISE_THRS)))
    plt.text(100, Displayed_[-1], PromNum)
Z = [NOISE_THRESHOLD_[j] for j in range(len(LED_INTENSITY_)) if PROMOTER_NAME_[j]=='SUSD4(2.4)-GCamp6s']
#plt.plot(X, np.array(Z))

plt.figure()
ax = plt.subplot(131)
ax2 = plt.subplot(132)
ax3 = plt.subplot(133)
ax.boxplot(PIXEL_PARAM_PER_PROMOTER, labels=LIST_)
ax.set_ylim(0, 1000)

ax2.scatter([LIST_[2] for i in range(len(ACTIVE_PIXEL_PARAM_PER_PROMOTER[2]))], ACTIVE_PIXEL_PARAM_PER_PROMOTER[2], alpha=0.2)
ax2.scatter([LIST_[3] for i in range(len(ACTIVE_PIXEL_PARAM_PER_PROMOTER[3]))], ACTIVE_PIXEL_PARAM_PER_PROMOTER[3], alpha=0.2)
ax2.scatter([LIST_[4] for i in range(len(ACTIVE_PIXEL_PARAM_PER_PROMOTER[4]))], ACTIVE_PIXEL_PARAM_PER_PROMOTER[4], alpha=0.2)

ax2.scatter([LIST_[6] for i in range(len(ACTIVE_PIXEL_PARAM_PER_PROMOTER[6]))], ACTIVE_PIXEL_PARAM_PER_PROMOTER[6], alpha=0.2)


ax3.scatter(ACTIVE_PIXEL_PARAM_PER_PROMOTER[2],PSD_PARAM_PER_PROMOTER[2], alpha=0.2)
ax3.scatter(ACTIVE_PIXEL_PARAM_PER_PROMOTER[3],PSD_PARAM_PER_PROMOTER[3], alpha=0.2)
ax3.scatter(ACTIVE_PIXEL_PARAM_PER_PROMOTER[4],PSD_PARAM_PER_PROMOTER[4], alpha=0.2)
ax3.scatter(ACTIVE_PIXEL_PARAM_PER_PROMOTER[6],PSD_PARAM_PER_PROMOTER[6], alpha=0.2)