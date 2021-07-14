# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 22:08:30 2021

@author: KEVIN-DORGANS
"""


"""
#Ex of importation method from Fiji
SPIKE_LIST_MANUAL_CROP.append(np.array(pd.read_clipboard())[:,1])
SPIKE_LIST_MANUAL_CROP_ID.append('AAV.PHP.eB.5HTr2b-tTA_AAV.PHP.eB.TRE.GCamp6s(10_11vgperml)_IO_2W_200nlx2_i.c._LED40_60fps2021-06-07-203558_00_N256 (IF1-CAM1)')
SPIKE_LIST_MANUAL_CROP_PROMOTER.append('AAV.PHP.eB.5HTr2b-tTA')




#Array_to_csv
PATH = load_directory_content__()[1]
pd.DataFrame(SPIKE_LIST_MANUAL_CROP).to_csv(PATH+'\FijiSomaROIcropped_SPIKES_PER_PROMOTER_rawfluo_db.csv')
pd.DataFrame(SPIKE_LIST_MANUAL_CROP_PROMOTER).to_csv(PATH+'\FijiSomaROIcropped_SPIKES_PER_PROMOTER_rawfluo_PromNum_db.csv')
pd.DataFrame(SPIKE_LIST_MANUAL_CROP_ID).to_csv(PATH+'\FijiSomaROIcropped_SPIKES_PER_PROMOTER_rawfluo_File_db.csv')




#Csv_to_array
PATH = load_directory_content__()[1]
SPIKE_LIST_MANUAL_CROP = []
try:
    temp_ = pd.read_csv(PATH+'\FijiSomaROIcropped_SPIKES_PER_PROMOTER_rawfluo_db.csv', index_col=0, header=0)
    temp_ = np.array(temp_.values)
    for i in range(len(temp_)):
        SPIKE_LIST_MANUAL_CROP.append([temp_[i][j] for j in range(len(temp_[i])) if str(temp_[i][j]) != 'nan'])
    temp_ = pd.read_csv(PATH+'\FijiSomaROIcropped_SPIKES_PER_PROMOTER_rawfluo_PromNum_db.csv')
    SPIKE_LIST_MANUAL_CROP_PROMOTER = np.array(temp_.values[:,1])
    temp_ = pd.read_csv(PATH+'\FijiSomaROIcropped_SPIKES_PER_PROMOTER_rawfluo_File_db.csv')
    SPIKE_LIST_MANUAL_CROP_ID = np.array(temp_.values[:,1])
except:
    print('file not recognized')
try:
    temp_ = pd.read_csv(PATH+'\FijiSomaROIcropped_SPIKES_PER_PROMOTER_rawfluo_Len_db.csv')
    SPIKE_LIST_MANUAL_CROP_LEN = np.array(temp_.values[:,1])
    temp_ = pd.read_csv(PATH+'\FijiSomaROIcropped_SPIKES_PER_PROMOTER_ALL_MicamDelay_db.csv')
    ALL_MicamDelay =  np.array(temp_.values[:,1])
except:
    print('no -LEN- -DELAY- files')
    SPIKE_LIST_MANUAL_CROP_LEN = []
    ALL_MicamDelay = []
"""

def log_func(x, a, b, c):
    return a * np.log2(b * x) + c

def line_func(x, a, b):
    return a * x + b


LED_INTENSITY_ = []
PROMOTER_NAME_ = []
FILE_NAME = []
ALL_SPIKES = []
ALL_PEAK_INDEXES = []
ALL_SPIKE_TIMESTAMPS = []


#FOR COMBINED ePhy/imaging EXPERIMENTS I USED [window=2sec, CaFluoResampling=49]
#FOR COMBINED imaging-onlyu [window=4sec, CaFluoResampling=150]

CalciumEventExtractionWindow = 4 #sec.
CalciumEventResampling = 150 #frame number for all window. Gets slower if bigger.
WindowTimeBeforePeak = 0.5 #Croping time before rise detection (s.)

#Search for all constructions/conditions used and make ref# index
LIST_ = []
for i in range(len(SPIKE_LIST_MANUAL_CROP_PROMOTER)):
    if (SPIKE_LIST_MANUAL_CROP_PROMOTER[i] in LIST_)==False:
        LIST_.append(SPIKE_LIST_MANUAL_CROP_PROMOTER[i])

#Trace-by-trace search for calcium events
for l in range(len(SPIKE_LIST_MANUAL_CROP)):  
    if True:
        PATH = SPIKE_LIST_MANUAL_CROP_ID[l]
        timestamp = int(PATH .split('_00')[0].split('-')[-1])

        try:
            try:
                MEAN = SPIKE_LIST_MANUAL_CROP[l]
                if ('SPIKE_LIST_MANUAL_CROP_LEN' in dir()) == True and len(SPIKE_LIST_MANUAL_CROP_LEN)>0:
                    FluoRecDuration = SPIKE_LIST_MANUAL_CROP_LEN[l]
                else:
                    FluoRecDuration = 10 #Default is 10sec. recording
                    
                SAMPLING_FREQUENCY = len(MEAN)/FluoRecDuration
                
                FILT = [np.nanmedian(MEAN[i+1:i+3])-MEAN[i] for i in range(len(MEAN)-4)]
                FILT = FILT - np.nanmean(FILT)
                FILT_ = FILT
                ZSCORE_FILT = sp.stats.zscore(FILT_)
                PEAKS = sp.signal.find_peaks(ZSCORE_FILT, height=3, distance= SAMPLING_FREQUENCY/1.5)

                FILT = MEAN
                plt.figure()
                plt.title(SPIKE_LIST_MANUAL_CROP_ID[l])
                ax=plt.subplot(131)
                ax2 = plt.subplot(132)
                ax3 = plt.subplot(133)
                
                ax3.plot(sp.stats.zscore(FILT_))
                
                PEAK_LIST = []
                PEAK_LIST_INDEXES = []
                PEAK_TIMESTAMPS = []
                
                ax.plot(MEAN)
                for i in range(len(PEAKS[0])):
                    ax.scatter(PEAKS[0][i], FILT[PEAKS[0][i]], color='red')
                    try:
                        if PEAKS[0][i]+SAMPLING_FREQUENCY*(CalciumEventExtractionWindow-0.5)<len(FILT) and PEAKS[0][i]-SAMPLING_FREQUENCY/2>0:
                            
                            temp_ = FILT[int(PEAKS[0][i]- SAMPLING_FREQUENCY/2): int(PEAKS[0][i]+ SAMPLING_FREQUENCY*(CalciumEventExtractionWindow-0.5))]
                            PEAK_LIST.append(temp_)
                            PEAK_LIST_INDEXES.append(PEAKS[0][i]/SAMPLING_FREQUENCY)
                            PEAK_TIMESTAMPS.append(timestamp)
                        elif int(PEAKS[0][i]- SAMPLING_FREQUENCY/1)<0: #THIS DOESNT WORK
                            missing = abs(int(PEAKS[0][i]- SAMPLING_FREQUENCY/2))
                            temp_ = np.concatenate([np.nan for l in range(missing)], temp_)
                            PEAK_LIST.append(temp_)
                            PEAK_LIST_INDEXES.append(PEAKS[0][i]/SAMPLING_FREQUENCY)
                            PEAK_TIMESTAMPS.append(timestamp)
                    except:
                        print('NoSpikeFound')
                        
                AVG = np.nanmean(PEAK_LIST, axis=0)
                SEM = sp.stats.sem(PEAK_LIST, axis=0)
                ax2.plot(np.linspace(0, FluoRecDuration, len(AVG)), AVG)
                ax2.fill_between(np.linspace(0, FluoRecDuration, len(AVG)), AVG+SEM, AVG-SEM, alpha=0.1)

                try:
                    LED_INTENSITY_.append(np.int(PATH.split('LED')[-1].split('_')[0]))
                except:
                    LED_INTENSITY_.append(np.nan)
                for k in range(len(PEAK_LIST)):
                    
                    PROMOTER_NAME_.append(SPIKE_LIST_MANUAL_CROP_PROMOTER[l])
                    FILE_NAME.append(SPIKE_LIST_MANUAL_CROP_ID[l])
                    ALL_PEAK_INDEXES.append(PEAK_LIST_INDEXES[k] )
                    ALL_SPIKES.append(sp.signal.resample(PEAK_LIST[k], CalciumEventResampling))
                    ALL_SPIKE_TIMESTAMPS.append(PEAK_TIMESTAMPS[k])
            except:
                pass
        except:
            pass
            
"""
plt.figure()
ax = plt.subplot(121)
ax2 = plt.subplot(122)

n_ALL_SPIKES = []
for i in range(len(ALL_SPIKES)):
    n_ALL_SPIKES.append(ALL_SPIKES[i]-np.nanmin(ALL_SPIKES[i][0:int(SAMPLING_FREQUENCY/2)]))

SEM = sp.stats.sem(ALL_SPIKES, axis=0)
MEAN = np.nanmean(ALL_SPIKES, axis=0)
X = np.linspace(0, FluoRecDuration, len(MEAN))
ax.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)
ax.plot(X, np.nanmean(ALL_SPIKES, axis=0))

SEM = sp.stats.sem(n_ALL_SPIKES, axis=0)
MEAN = np.nanmean(n_ALL_SPIKES, axis=0)
X = np.linspace(0, FluoRecDuration, len(MEAN))
ax2.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)
ax2.plot(X, MEAN)
for i in range(len(n_ALL_SPIKES)):
    ax2.plot(X, n_ALL_SPIKES[i], color='black', lw=0.2)
"""

#LIST_ = ['AAV.PHP.eB.5HTr2b-tTA', 'AAV.PHP.S.5HTr2b-tTA', 'AAV9.5HTr2b-tTA','AAV9.5HTr2b(1.8)','AAV9.SUSD4(2.4)']
#LIST_ = ['AAV9.5HTr2b-tTA', 'AAV9.CAG']

ALL_F_ZERO_PROM = []
ALL_DELTA_SPIKE_PROM = []
ALL_RISE_TIME_SPIKE_PROM = []
ALL_WAVE_DECAY_FACTOR_PROM  = []
ALL_SPIKE_WAVES_PROM = []
PROMOTER_NAME_ALL = []
ALL_TIMESTAMPS_PROM = []
ALL_PEAK_TIMES_PROM = []


plt.figure()
for i in range(len(LIST_)):
    ax = plt.subplot(1, len(LIST_), i+1)

    AVG = [ALL_SPIKES[j]*0.0014665418*SAMPLING_FREQUENCY for j in range(len(ALL_SPIKES)) if PROMOTER_NAME_[j] == LIST_[i]]
    PEAK_TIMES = [ALL_PEAK_INDEXES[j] for j in range(len(ALL_SPIKES)) if PROMOTER_NAME_[j] == LIST_[i]]
    SORTED_TIMESTAMPS = [ALL_SPIKE_TIMESTAMPS[j] for j in range(len(ALL_SPIKES)) if PROMOTER_NAME_[j] == LIST_[i]]
    
    if len(AVG)>0:
        
        #DELTA_MAX_SPIKE in fW
        #F_ZERO = baseline = f.intensity just before events in fW
        #DELTA_MAX_SPIKE_INDEX as integer index value 
        #RISE_TIME_SPIKE in sec.
        #DECAY_SLOPE_VALUE in fW
        #DECAY_SLOPE_INDEX as integer index value
        #DECAY_SLOPE_FACTOR in fW/s
        
        PeakDetectionWindow_start = int(0.5*CalciumEventResampling/CalciumEventExtractionWindow +1) #events are aligned at 0.5s.
        PeakDetectionWindow_end = int(1*CalciumEventResampling/CalciumEventExtractionWindow +1) #event rise time cannot be longer than 0.5s.
        BSL_window = int(0.1*CalciumEventResampling/CalciumEventExtractionWindow ) #calcium fluctuates at rates > 2Hz so bsl is taken for 200ms
        BSL_window_start = int(0.5*CalciumEventResampling/CalciumEventExtractionWindow)-BSL_window-2
        DecayWindow_start = int(1*CalciumEventResampling/CalciumEventExtractionWindow)
        DecayWindow_end = int(1.5*CalciumEventResampling/CalciumEventExtractionWindow)
        
        F_ZERO = [np.nanmin(AVG[j][BSL_window_start : PeakDetectionWindow_start]) for j in range(len(AVG))]
        DELTA_MAX_SPIKE = [np.nanmax(AVG[j][PeakDetectionWindow_start : PeakDetectionWindow_end])-F_ZERO[j] for j in range(len(AVG))]
        DELTA_MAX_SPIKE_INDEX = [AVG[j].tolist().index(DELTA_MAX_SPIKE[j]+F_ZERO[j]) for j in range(len(AVG))]
        RISE_TIME_SPIKE = [DELTA_MAX_SPIKE_INDEX[j]/(int(CalciumEventResampling)/CalciumEventExtractionWindow)-0.5 for j in range(len(AVG))]
        DECAY_SLOPE_VALUE =  [AVG[j][int(DELTA_MAX_SPIKE_INDEX[j])+int(CalciumEventResampling/CalciumEventExtractionWindow)]-F_ZERO[j] for j in range(len(AVG))]
        DECAY_SLOPE_INDEX = [int(DELTA_MAX_SPIKE_INDEX[j])+int(CalciumEventResampling)/CalciumEventExtractionWindow for j in range(len(AVG))]
        DECAY_SLOPE_FACTOR = [(DECAY_SLOPE_VALUE[j]-DELTA_MAX_SPIKE[j])/(DECAY_SLOPE_INDEX[j]/(int(CalciumEventResampling)/CalciumEventExtractionWindow)-DELTA_MAX_SPIKE_INDEX[j]/(int(CalciumEventResampling)/CalciumEventExtractionWindow)) for j in range(len(AVG))]

        
        MEAN = np.nanmean(AVG, axis=0)
        SEM = sp.stats.sem(AVG, axis=0)
        X = np.linspace(0, CalciumEventExtractionWindow, len(MEAN))
        
        
        #plt.plot(X, MEAN, label=LIST_[i])
        #plt.fill_between(X, MEAN+SEM, MEAN-SEM, color='black', alpha=0.01)
        print(LIST_[i],'n=',len(AVG), 'max=', np.nanmax(MEAN)-MEAN[10])
        
        ALL_F_ZERO = []
        ALL_DELTA_SPIKE = []
        ALL_SPIKES_PER_PROMOTER = []
        ALL_RISE_TIME_SPIKE = []
        ALL_DECAY_SLOPE_FACTOR = []
        ALL_SPIKE_WAVES = []
        ALL_TIMESTAMPS = []
        ALL_PEAK_TIMES = []
        
        for k in range(len(AVG)):
            #SOME LINE FIT TO SUBTRACT BLEACH DECAY ON SHORT PERIOD
            #MOSTLY AESTHETIC, Ca FLUCTUATES TOO MUCH ON BSL FOR CORRECT FIT OF ANYTHING
            
            y_data = []
            x_data = []
            for l in range(len(AVG[k])):
                if l<int(0.5*CalciumEventResampling/CalciumEventExtractionWindow) or CalciumEventResampling-(CalciumEventResampling*0.15)<l:
                    y_data.append(AVG[k][l])
                    x_data.append(l)
            popt, pcov = curve_fit(line_func, x_data, y_data, maxfev=1000)
            LINE_FIT = line_func(np.linspace(0, len(AVG[k]), len(AVG[k])), *popt)
            AVG[k] = AVG[k] - LINE_FIT            
            AVG[k] = AVG[k] - np.nanmedian(AVG[k][PeakDetectionWindow_start-1-BSL_window:PeakDetectionWindow_start-1])
            MEAN = AVG[k] - F_ZERO[k]
            #!!!FILTER!!! 
            if np.nanmean(AVG[k][int(1.8*CalciumEventResampling/CalciumEventExtractionWindow):int(2*CalciumEventResampling/CalciumEventExtractionWindow)])<(DELTA_MAX_SPIKE[k])/1.4:
            #if True:
                ax.plot(X, AVG[k] , color='black', lw=0.1)
                ax.scatter(RISE_TIME_SPIKE[k]+0.5, DELTA_MAX_SPIKE[k], color='red')
                ALL_F_ZERO.append(F_ZERO[k])
                ALL_DELTA_SPIKE.append(DELTA_MAX_SPIKE[k])
                ALL_SPIKES_PER_PROMOTER.append(MEAN)
                ALL_RISE_TIME_SPIKE.append(RISE_TIME_SPIKE[k])
                ALL_DECAY_SLOPE_FACTOR.append(DECAY_SLOPE_FACTOR[k])
                ALL_SPIKE_WAVES.append(AVG[k])
                ALL_PEAK_TIMES.append(PEAK_TIMES[k])
                ALL_TIMESTAMPS.append(SORTED_TIMESTAMPS[k])
                
    MEAN = np.nanmean(ALL_SPIKE_WAVES, axis=0)
    SEM = sp.stats.sem(ALL_SPIKE_WAVES, axis=0)
    ax.plot(X, np.nanmean(ALL_SPIKE_WAVES, axis=0), label=LIST_[i])
    ax.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)
    
    if len(ALL_SPIKE_WAVES)>0:
        ALL_F_ZERO_PROM.append(ALL_F_ZERO)
        ALL_DELTA_SPIKE_PROM.append(ALL_DELTA_SPIKE)
        ALL_RISE_TIME_SPIKE_PROM.append(ALL_RISE_TIME_SPIKE)
        ALL_WAVE_DECAY_FACTOR_PROM.append(ALL_DECAY_SLOPE_FACTOR)
        ALL_SPIKE_WAVES_PROM.append(ALL_SPIKE_WAVES)
        PROMOTER_NAME_ALL.append([LIST_[i] for k in range(len(RISE_TIME_SPIKE))])
        ALL_TIMESTAMPS_PROM.append(ALL_TIMESTAMPS)
        ALL_PEAK_TIMES_PROM.append(ALL_PEAK_TIMES)
    
plt.xlabel('Time(s.)')
plt.legend()


plt.figure(figsize= (17, 3), num='Recap_LEDPWR')
ax = plt.subplot(171)
ax2 = plt.subplot(172)
ax3 = plt.subplot(173)
ax4 = plt.subplot(174)
ax5 = plt.subplot(175)
ax6 = plt.subplot(176)
ax7 = plt.subplot(177)

X_fit = []
Y_fit = []

for i in range(len(ALL_SPIKE_WAVES_PROM)):
    MEAN = np.nanmean(ALL_SPIKE_WAVES_PROM[i], axis=0)
    SEM = sp.stats.sem(ALL_SPIKE_WAVES_PROM[i], axis=0)
    X = np.linspace(0, CalciumEventExtractionWindow, len(MEAN))
    ax.plot(X, MEAN, lw=1)
    ax.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)
    
for i in range(len(ALL_SPIKE_WAVES_PROM)):
    normed_ = [ALL_SPIKE_WAVES_PROM[i][j]/ALL_F_ZERO_PROM[i][j] for j in range(len(ALL_SPIKE_WAVES_PROM[i]))]
    MEAN = np.nanmean(normed_, axis=0)
    SEM = sp.stats.sem(normed_, axis=0)
    X = np.linspace(0, CalciumEventExtractionWindow, len(MEAN))
    ax2.plot(X, MEAN, lw=1)
    ax2.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)
    
for i in range(len(ALL_DELTA_SPIKE_PROM)):
    ax3.scatter(np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(ALL_DELTA_SPIKE_PROM[i]))
    SEM = sp.stats.sem(ALL_DELTA_SPIKE_PROM[i])
    MEAN = np.nanmean(ALL_DELTA_SPIKE_PROM[i])
    ax3.plot((np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(ALL_F_ZERO_PROM[i])), (MEAN+SEM, MEAN-SEM), color='black', lw=0.1)
    SEM = sp.stats.sem(ALL_F_ZERO_PROM[i])
    MEAN = np.nanmean(ALL_F_ZERO_PROM[i])
    ax3.plot((MEAN+SEM, MEAN-SEM), (np.nanmean(ALL_DELTA_SPIKE_PROM[i]), np.nanmean(ALL_DELTA_SPIKE_PROM[i])), color='black', lw=0.1)

    X_fit.append(np.nanmean(ALL_F_ZERO_PROM[i]))
    Y_fit.append(np.nanmean(ALL_DELTA_SPIKE_PROM[i]))
    
try:
    popt, pcov = curve_fit(log_func , X_fit, Y_fit, maxfev=1000)
    x_ = np.linspace(0, np.nanmax(X_fit), CalciumEventResampling)
    ax3.plot(x_, log_func(x_, *popt), color='black', lw=0.1)
    
    popt, pcov = curve_fit(line_func , [X_fit[i] for i in range(len(X_fit)) if X_fit[i] <500], [Y_fit[i] for i in range(len(Y_fit)) if X_fit[i] <500], maxfev=1000)
    x_ = np.linspace(0, np.nanmax(X_fit), CalciumEventResampling)
    ax3.plot(x_, line_func(x_, *popt), color='black', lw=0.1)
except:
    pass

for i in range(len(ALL_DELTA_SPIKE_PROM)):
    ax4.scatter(np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(np.divide(ALL_DELTA_SPIKE_PROM[i],ALL_F_ZERO_PROM[i])), label=LIST_[i])
    SEM = sp.stats.sem(np.divide(ALL_DELTA_SPIKE_PROM[i],ALL_F_ZERO_PROM[i]))
    MEAN = np.nanmean(np.divide(ALL_DELTA_SPIKE_PROM[i],ALL_F_ZERO_PROM[i]))
    ax4.plot((np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(ALL_F_ZERO_PROM[i])), (MEAN+SEM, MEAN-SEM), color='black', lw=0.1)
    SEM = sp.stats.sem(ALL_F_ZERO_PROM[i])
    MEAN = np.nanmean(ALL_F_ZERO_PROM[i])
    MEAN_Y = np.nanmean(np.divide(ALL_DELTA_SPIKE_PROM[i],ALL_F_ZERO_PROM[i]))
    ax4.plot((MEAN+SEM, MEAN-SEM), (MEAN_Y , MEAN_Y ), color='black', lw=0.1)

for i in range(len(ALL_DELTA_SPIKE_PROM)):
    ax5.scatter(np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(ALL_RISE_TIME_SPIKE_PROM[i]), label=LIST_[i])
    SEM = sp.stats.sem(np.nanmean(ALL_RISE_TIME_SPIKE_PROM[i]))
    MEAN = np.nanmean(np.nanmean(ALL_RISE_TIME_SPIKE_PROM[i]))
    ax5.plot((np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(ALL_F_ZERO_PROM[i])), (MEAN+SEM, MEAN-SEM), color='black', lw=0.1)
    SEM = sp.stats.sem(ALL_F_ZERO_PROM[i])
    MEAN = np.nanmean(ALL_F_ZERO_PROM[i])
    MEAN_Y = np.nanmean(np.nanmean(ALL_RISE_TIME_SPIKE_PROM[i]))
    ax5.plot((MEAN+SEM, MEAN-SEM), (MEAN_Y , MEAN_Y ), color='black', lw=0.1)

for i in range(len(ALL_DELTA_SPIKE_PROM)):
    ax6.scatter(np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(ALL_WAVE_DECAY_FACTOR_PROM[i]), label=LIST_[i])
    SEM = sp.stats.sem(np.nanmean(ALL_WAVE_DECAY_FACTOR_PROM[i]))
    MEAN = np.nanmean(np.nanmean(ALL_WAVE_DECAY_FACTOR_PROM[i]))
    ax6.plot((np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(ALL_F_ZERO_PROM[i])), (MEAN+SEM, MEAN-SEM), color='black', lw=0.1)
    SEM = sp.stats.sem(ALL_F_ZERO_PROM[i])
    MEAN = np.nanmean(ALL_F_ZERO_PROM[i])
    MEAN_Y = np.nanmean(np.nanmean(ALL_WAVE_DECAY_FACTOR_PROM[i]))
    ax6.plot((MEAN+SEM, MEAN-SEM), (MEAN_Y , MEAN_Y ), color='black', lw=0.1)

for i in range(len(ALL_DELTA_SPIKE_PROM)):
    ax7.scatter(np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(np.divide(ALL_WAVE_DECAY_FACTOR_PROM[i], ALL_DELTA_SPIKE_PROM[i])), label=LIST_[i])
    SEM = sp.stats.sem(np.nanmean(np.divide(ALL_WAVE_DECAY_FACTOR_PROM[i], ALL_DELTA_SPIKE_PROM[i])))
    MEAN = np.nanmean(np.nanmean(np.divide(ALL_WAVE_DECAY_FACTOR_PROM[i], ALL_DELTA_SPIKE_PROM[i])))
    ax7.plot((np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(ALL_F_ZERO_PROM[i])), (MEAN+SEM, MEAN-SEM), color='black', lw=0.1)
    SEM = sp.stats.sem(ALL_F_ZERO_PROM[i])
    MEAN = np.nanmean(ALL_F_ZERO_PROM[i])
    MEAN_Y = np.nanmean(np.nanmean(np.divide(ALL_WAVE_DECAY_FACTOR_PROM[i], ALL_DELTA_SPIKE_PROM[i])))
    ax7.plot((MEAN+SEM, MEAN-SEM), (MEAN_Y , MEAN_Y ), color='black', lw=0.1)

ax3.set_xlabel('Somatic F-Zero (fW)')
ax3.set_ylabel('Calcium Spike Amplitude (fW)')
ax4.set_ylabel('Calcium Spike DF/F0')
ax4.set_xlabel('Somatic F-Zero (fW)')
ax5.set_ylabel('Time-to-peak (s.)')
ax6.set_ylabel('Event decay factor (fW/s)')
ax3.set_xlim(0,)
ax3.set_ylim(0,)
plt.legend()
plt.tight_layout()


plt.figure(figsize=(15, 10))
for i in range(len(LIST_)):
    ax = plt.subplot(5, len(LIST_), i+1)
    ax2 = plt.subplot(5, len(LIST_), i+1+len(LIST_))
    ax3 = plt.subplot(5, len(LIST_), i+1+2*len(LIST_))
    ax4 = plt.subplot(5, len(LIST_), i+1+3*len(LIST_))
    ax5 = plt.subplot(5, len(LIST_), i+1+4*len(LIST_))
    
    for j in range(len(LIST_)):
        ax.scatter(np.nanmean(ALL_F_ZERO_PROM[j]), np.nanmean(ALL_RISE_TIME_SPIKE_PROM[j]), color='black', alpha=0.3)
        ax2.scatter(np.nanmean(ALL_F_ZERO_PROM[j]), np.nanmean(ALL_DELTA_SPIKE_PROM[j]), color='black', alpha=0.3)
        ax3.scatter(np.nanmean(ALL_RISE_TIME_SPIKE_PROM[j]), np.nanmean(np.divide(ALL_DELTA_SPIKE_PROM[j], ALL_F_ZERO_PROM[j])), color='black', alpha=0.3)
        ax2.scatter(ALL_F_ZERO_PROM[j], ALL_DELTA_SPIKE_PROM[j], s=5, color='black', alpha=0.1)
        ax3.scatter(ALL_RISE_TIME_SPIKE_PROM[j], np.divide(ALL_DELTA_SPIKE_PROM[j], ALL_F_ZERO_PROM[j]), s=5, color='black', alpha=0.1)
        ax.scatter(ALL_F_ZERO_PROM[j], ALL_RISE_TIME_SPIKE_PROM[j], s=5, color='black', alpha=0.1)
        ax4.scatter(np.nanmean(ALL_F_ZERO_PROM[j]), np.nanmean(np.divide(ALL_DELTA_SPIKE_PROM[j], ALL_F_ZERO_PROM[j])), color='black', alpha=0.3)
        ax4.scatter(ALL_F_ZERO_PROM[j]/np.nanmax(ALL_F_ZERO_PROM[j]), np.divide(ALL_DELTA_SPIKE_PROM[j], ALL_F_ZERO_PROM[j]), s=5, color='black', alpha=0.1)
        ax5.scatter(np.nanmean(ALL_F_ZERO_PROM[j]), np.nanmean(ALL_WAVE_DECAY_FACTOR_PROM[j]), color='black', alpha=0.3)
        ax5.scatter(ALL_F_ZERO_PROM[j], ALL_WAVE_DECAY_FACTOR_PROM[j], s=5, color='black', alpha=0.1)

    
    ax.scatter(ALL_F_ZERO_PROM[i], ALL_RISE_TIME_SPIKE_PROM[i], s=5)
    ax.scatter(np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(ALL_RISE_TIME_SPIKE_PROM[i]))
    ax.set_xlim(0, np.nanmax(np.concatenate(ALL_F_ZERO_PROM)))
    ax.set_ylim(0, np.nanmax(np.concatenate(ALL_RISE_TIME_SPIKE_PROM)))
    ax.set_xlabel('Baseline Fluo. (fW)')
    ax.set_ylabel('Time-to-peak (s.)')
    
    ax2.scatter(ALL_F_ZERO_PROM[i], ALL_DELTA_SPIKE_PROM[i], s=5)
    ax2.scatter(np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(ALL_DELTA_SPIKE_PROM[i]))
    ax2.set_xlim(0, np.nanmax(np.concatenate(ALL_F_ZERO_PROM)))
    ax2.set_ylim(0, np.nanmax(np.concatenate(ALL_DELTA_SPIKE_PROM)))
    ax2.set_xlabel('Baseline Fluo. (fW)')
    ax2.set_ylabel('Peak Intensity (fW)')
    
    ax3.scatter(ALL_RISE_TIME_SPIKE_PROM[i], np.divide(ALL_DELTA_SPIKE_PROM[i], ALL_F_ZERO_PROM[i]), s=5)
    ax3.scatter(np.nanmean(ALL_RISE_TIME_SPIKE_PROM[i]), np.nanmean(np.divide(ALL_DELTA_SPIKE_PROM[i], ALL_F_ZERO_PROM[i])))
    ax3.set_xlim(0, np.nanmax(np.concatenate(ALL_RISE_TIME_SPIKE_PROM)))
    ax3.set_ylabel('F/F0 DFF)')
    ax3.set_xlabel('Time-to-peak (s.)')
    
    ax4.scatter(ALL_F_ZERO_PROM[i], np.divide(ALL_DELTA_SPIKE_PROM[i], ALL_F_ZERO_PROM[i]), s=5)
    ax4.scatter(np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(np.divide(ALL_DELTA_SPIKE_PROM[i], ALL_F_ZERO_PROM[i])))
    ax4.set_xlim(0, np.nanmax(np.concatenate(ALL_F_ZERO_PROM)))
    ax4.set_xlabel('Baseline Fluo. (fW)')
    ax4.set_ylabel('F/F0 DFF)')
    
    ax5.scatter(ALL_F_ZERO_PROM[i], ALL_WAVE_DECAY_FACTOR_PROM[i], s=5)
    ax5.scatter(np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(ALL_WAVE_DECAY_FACTOR_PROM[i]))
    ax5.set_xlim(0, np.nanmax(np.concatenate(ALL_F_ZERO_PROM)))
    ax5.set_ylim(0, np.nanmin(np.concatenate(ALL_WAVE_DECAY_FACTOR_PROM)))
    ax5.set_xlabel('Baseline Fluo. (fW)')
    ax5.set_ylabel('Decay factor (fW/s)')
        
plt.tight_layout()




"""
#THIS GETS BACK ePhy patch-clamped cells at the time of the previously extracted CaSignals

ALL_TIMESTAMPS_ePhy = []
ALL_MEMBRANE_VOLT = []
ALL_ePhy_SIGNALS = []
ALL_sampling_freqs_ephy = []

PATHS_ = load_directory_content__()[0]
for PATH in PATHS_:
    R = easy_IbwLoad(DOWNSAMPLING__=1, _file_ = PATH)
    RAW = R[0]
    SAMPLING_FREQUENCY = 1 / R[1]
    
    timestamp = int(PATH.split('.ibw')[0].split('_')[-1])
    #plt.plot(np.linspace(0, len(RAW)*SAMPLING_FREQUENCY, len(RAW)), RAW)
    
    ALL_TIMESTAMPS_ePhy.append(timestamp)
    ALL_MEMBRANE_VOLT.append(np.nanmedian(RAW))
    ALL_ePhy_SIGNALS.append(RAW)
    ALL_sampling_freqs_ephy.append(SAMPLING_FREQUENCY)


FluoEventID = np.concatenate(ALL_TIMESTAMPS_PROM)
FluoEventTimes = np.concatenate(ALL_PEAK_TIMES_PROM)
AllCalciumWaves = np.concatenate(ALL_SPIKE_WAVES_PROM)
FluoEventRiseTimes = np.concatenate(ALL_RISE_TIME_SPIKE_PROM)
FluoEventDecayConstant = np.concatenate(ALL_WAVE_DECAY_FACTOR_PROM)
FluoEventAmplitude = np.concatenate(ALL_DELTA_SPIKE_PROM)
FluoEventFZERO = np.concatenate(ALL_F_ZERO_PROM)

ePhyEvent_times = []
ePhy_associated_Fluo_trace = []
CalciumEvent_RT = []
CalciumEvent_DC = []
CalciumEvent_F = []
CalciumEvent_FZERO = []

for i in range(len(ALL_TIMESTAMPS_ePhy)):
    try:
        ePhyTimes = [FluoEventTimes[j]+1 for j in range(len(FluoEventID)) if ALL_TIMESTAMPS_ePhy[i]==FluoEventID[j]][0]
        Calcium_signal = [AllCalciumWaves[j] for j in range(len(FluoEventID)) if ALL_TIMESTAMPS_ePhy[i]==FluoEventID[j]][0]
        Gcamp6sRT = [FluoEventRiseTimes[j] for j in range(len(FluoEventID)) if ALL_TIMESTAMPS_ePhy[i]==FluoEventID[j]][0]
        Gcamp6sDC = [FluoEventDecayConstant[j] for j in range(len(FluoEventID)) if ALL_TIMESTAMPS_ePhy[i]==FluoEventID[j]][0]
        Gcamp6sF = [FluoEventAmplitude[j] for j in range(len(FluoEventID)) if ALL_TIMESTAMPS_ePhy[i]==FluoEventID[j]][0]
        Gcamp6sFZERO = [FluoEventFZERO[j] for j in range(len(FluoEventID)) if ALL_TIMESTAMPS_ePhy[i]==FluoEventID[j]][0]
    except:
        ePhyTimes = np.nan
        Calcium_signal = np.nan
        Gcamp6sRT = np.nan
        Gcamp6sDC = np.nan
        Gcamp6sF = np.nan
        Gcamp6sFZERO = np.nan
        
    ePhyEvent_times.append(ePhyTimes)
    ePhy_associated_Fluo_trace.append(Calcium_signal)
    CalciumEvent_RT.append(Gcamp6sRT)
    CalciumEvent_DC.append(Gcamp6sDC)
    CalciumEvent_F.append(Gcamp6sF)
    CalciumEvent_FZERO.append(Gcamp6sFZERO)
    #temp_indexes_ephy = [(FluoEventTimes[j]-1.2)/ALL_sampling_freqs_ephy[j] for j in range(len(ALL_TIMESTAMPS_ePhy)) if ALL_TIMESTAMPS_ePhy[j]==X[i]][0]
    
    #temp_trace_ephy = [ALL_ePhy_SIGNALS[j][int(temp_indexes_ephy)-int(0.5/ALL_sampling_freqs_ephy[j]):int(temp_indexes_ephy)+int(0.5/ALL_sampling_freqs_ephy[j])] for j in range(len(ALL_TIMESTAMPS_ePhy)) if ALL_TIMESTAMPS_ePhy[j]==X[i]][0]
    #plt.plot(temp_trace_ephy)



plt.figure(figsize=(2,6))
ax = plt.subplot(311)
ax2 = plt.subplot(312, sharex=ax)
ax3 = plt.subplot(313)

ePhy_toMean = []
Calcium_toMean = []
ePhy_Spike_Delay = []
ePhy_Spike_Traces = []

FILTERED_CalciumEvent_RT = []
FILTERED_CalciumEvent_DC = []
FILTERED_CalciumEvent_F = []
FILTERED_CalciumEvent_FZERO = []

for i in range(len(ePhyEvent_times)):
    if ePhyEvent_times[i]>0 and np.nanmedian(ALL_ePhy_SIGNALS[i])<-30:
        window_start = int((ePhyEvent_times[i]-0.5+ALL_MicamDelay[i])/ALL_sampling_freqs_ephy[i])
        window_end = int((ePhyEvent_times[i]+CalciumEventExtractionWindow-0.5+ALL_MicamDelay[i])/ALL_sampling_freqs_ephy[i])
        temp_trace_ephy = ALL_ePhy_SIGNALS[i][window_start:window_end]
        
        SpikeExtracted = temp_trace_ephy[int(0.3/ALL_sampling_freqs_ephy[i]):int(0.75/ALL_sampling_freqs_ephy[i])]
        Derivated = [SpikeExtracted[i+5]-SpikeExtracted[i] for i in range(len(SpikeExtracted)-6) if SpikeExtracted[i]!=0]
        SpikePeak = np.nanmax(Derivated)
        SpikePeakIndex = Derivated.index(SpikePeak)
        
        ePhy_Spike_Delay.append(SpikePeakIndex*ALL_sampling_freqs_ephy[i])
        try:
            ePhy_Spike_Traces.append(sp.signal.resample(temp_trace_ephy[SpikePeakIndex-int(0.01/ALL_sampling_freqs_ephy[i]):SpikePeakIndex+int(0.8/ALL_sampling_freqs_ephy[i])], 3000))
        except:
            ePhy_Spike_Traces.append([np.nan for k in range(3000)])
                
        temp_trace_ephy = temp_trace_ephy-np.nanmedian(temp_trace_ephy)
        temp_trace_ephy = temp_trace_ephy/np.nanmax(temp_trace_ephy)
        
        X = np.linspace(0, CalciumEventExtractionWindow, len(temp_trace_ephy))
        ax2.plot(X, temp_trace_ephy+i, color='black', lw=0.1)
        
        temp_trace_calcium = ePhy_associated_Fluo_trace[i]
        temp_trace_calcium = temp_trace_calcium-np.nanmedian(temp_trace_calcium[0:15])
        temp_trace_calcium = temp_trace_calcium/np.nanmax(temp_trace_calcium)
        X = np.linspace(0, CalciumEventExtractionWindow, len(temp_trace_calcium))
        ax.plot(X, temp_trace_calcium, color='black', lw=0.1)
        
        ePhy_toMean.append(sp.signal.resample(temp_trace_ephy, CalciumEventResampling))
        Calcium_toMean.append(temp_trace_calcium)
        FILTERED_CalciumEvent_RT.append(CalciumEvent_RT[i])
        FILTERED_CalciumEvent_DC.append(CalciumEvent_DC[i])
        FILTERED_CalciumEvent_F.append(CalciumEvent_F[i])
        FILTERED_CalciumEvent_FZERO.append(CalciumEvent_FZERO[i])

ax.plot((0.5, 0.5), (-0.2,1), color='black', lw=1.5, alpha=0.5 )
ax2.plot((0.5, 0.5), (0,len(ePhyEvent_times)), color='black', lw=1.5, alpha=0.5 )

MEAN = np.nanmean(np.array(ePhy_Spike_Traces),  axis=0)
SEM = sp.stats.sem(np.array(ePhy_Spike_Traces),  axis=0)
X = np.linspace(-0.1, 0.4, len(MEAN))
for i in range(len(ePhy_Spike_Traces)):
    ax3.plot(X, ePhy_Spike_Traces[i], color='black', lw=0.1)
ax3.plot(X, MEAN)
ax3.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)

MEAN = np.nanmean(np.array(Calcium_toMean),  axis=0)
X = np.linspace(0, CalciumEventExtractionWindow, len(MEAN))
SEM = sp.stats.sem(np.array(Calcium_toMean),  axis=0)
ax.plot(X, MEAN)
ax.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)

ax.set_xlabel('Time(s.)')
ax.set_ylabel('GCamp6s s.events (/Peak)')
ax2.set_xlabel('Time(s.)')
ax2.set_ylabel('cclamp voltage (mV)')
ax3.set_xlabel('Peak-aligned(s.)')
ax3.set_ylabel('Mb. voltage (mV)')
plt.tight_layout()

"""



"""
#CUMULATIVE HISTOGRAMS FOR ALL CALCIUM SPIKE PARAMETERS
plt.figure(figsize=(15, 3))
ax = plt.subplot(151)
ax2 = plt.subplot(152)
ax3 = plt.subplot(153)
ax4 = plt.subplot(154)
ax5 = plt.subplot(155)

for i in range(len(LIST_)):
    ax.hist(ALL_F_ZERO_PROM[i], alpha=0.8, bins=200, cumulative=True, density=True,  histtype='step', label=LIST_[i])
    ax2.hist(ALL_RISE_TIME_SPIKE_PROM[i], alpha=0.8,bins=200, cumulative=True, density=True, histtype='step', label=LIST_[i])
    ax3.hist(ALL_DELTA_SPIKE_PROM[i], alpha=0.8,bins=200,cumulative=True,  density=True, histtype='step', label=LIST_[i])
    ax4.hist(np.add(ALL_DELTA_SPIKE_PROM[i], ALL_F_ZERO_PROM[i]), alpha=0.8, cumulative=True, density=True, histtype='step' , label=LIST_[i])
    ax5.hist(ALL_WAVE_DECAY_FACTOR_PROM[i], alpha=0.8,bins=200, cumulative=True, density=True, histtype='step', label=LIST_[i])

ax.set_xlabel('Baseline intensity (fW)')
ax2.set_xlabel('Time-to-peak (s.)')
ax3.set_xlabel('Peak intensity (fW)')
ax4.set_xlabel('DF/F')
ax5.set_xlabel('Decay constant(fW/s)')
ax.set_xscale('log')
ax3.set_xscale('log')
plt.tight_layout()


#BOXPLOTS FOR ALL CA-SPKE PARAMS
plt.figure(figsize=(15, 3))
ax = plt.subplot(151)
ax2 = plt.subplot(152)
ax3 = plt.subplot(153)
ax4 = plt.subplot(154)
ax5 = plt.subplot(155)

ax.boxplot(ALL_F_ZERO_PROM, labels=LIST_, vert=False)
ax2.boxplot(ALL_RISE_TIME_SPIKE_PROM, labels=LIST_, vert=False)
ax3.boxplot(ALL_DELTA_SPIKE_PROM, labels=LIST_, vert=False)
ax4.boxplot([np.divide(ALL_DELTA_SPIKE_PROM[i], ALL_F_ZERO_PROM[i]) for i in range(len(LIST_))], labels=LIST_, vert=False)
ax5.boxplot(ALL_WAVE_DECAY_FACTOR_PROM, labels=LIST_, vert=False)

ax.set_ylabel('Baseline intensity (fW)')
ax2.set_ylabel('Time-to-peak (s.)')
ax3.set_ylabel('Peak intensity (fW)')
ax4.set_ylabel('DF/F')
ax5.set_ylabel('Decay constant(fW/s)')
ax.set_xscale('log')
ax3.set_xscale('log')
"""

"""
#FIXED GCamp6s Analysis
PATH = load_directory_content__()[1]
ALL_FIXED_TISSUE_CELL_INTENSITIES = []
ALL_FIXED_TISSUE_CELL_INTENSITIES_PROM = []
try:
    temp_ = pd.read_csv(PATH+'\FijiSomaROIcropped_5x_ALL_FIXED_TISSUE_CELL_INTENSITIES_db.csv', index_col=0, header=0)
    temp_ = np.array(temp_.values)
    for i in range(len(temp_)):
        ALL_FIXED_TISSUE_CELL_INTENSITIES.append([temp_[i][j] for j in range(len(temp_[i])) if str(temp_[i][j]) != 'nan'])
    
    temp_ = pd.read_csv(PATH+'\FijiSomaROIcropped_5x_ALL_FIXED_TISSUE_CELL_INTENSITIES_PROM_db.csv')
    ALL_FIXED_TISSUE_CELL_INTENSITIES_PROM = np.array(temp_.values[:,1])
except:
    print('file not recognized')

plt.figure(figsize=(4,2))
ax = plt.subplot(141)
ax2 = plt.subplot(142)
ax3 = plt.subplot(143)
ax4 = plt.subplot(144)

ALL_TEMP = []

for i in range(len(LIST_)):
    temp_ = [ALL_FIXED_TISSUE_CELL_INTENSITIES[j] for j in range(len(ALL_FIXED_TISSUE_CELL_INTENSITIES_PROM)) if ALL_FIXED_TISSUE_CELL_INTENSITIES_PROM[j]==LIST_[i]]
    temp_2 = ALL_F_ZERO_PROM[i]
    ALL_TEMP.append(np.concatenate(temp_))
    
    ax.hist(np.concatenate(temp_), histtype='step', cumulative=True, density=True, bins=200)
    ax2.hist(temp_2, histtype='step', cumulative=True, density=True, bins=200)
    MEAN = np.nanmean(np.concatenate(temp_))
    SEM = sp.stats.sem(np.concatenate(temp_))
    MEAN_2 = np.nanmean(temp_2)
    SEM_2 = sp.stats.sem(temp_2)
    ax3.plot((MEAN, MEAN), (MEAN_2+SEM_2, MEAN_2-SEM_2), color='black')
    ax3.plot((MEAN+SEM, MEAN-SEM), (MEAN_2, MEAN_2), color='black')
    ax3.scatter(np.nanmean(np.concatenate(temp_)), MEAN_2)
ax4.boxplot(ALL_TEMP, labels=LIST_, vert=False)
ax.set_xlabel('Fixed GCamp6s Fluo')
ax2.set_xlabel('InVitro FZero')
ax3.set_xlabel('Fixed GCamp6s Fluo')
ax3.set_ylabel('InVitro FZero')
ax4.set_xlabel('log')
"""



"""
#REVERSE APPROACH : FINDS SPIKES/EVENTS FIRST, TAKES CALCIUM FROM SPIKES

DETECTED_ePHY_SPIKES_WF = []
DETECTED_ePHY_SPIKES_TIMESTAMP = []
DETECTED_ePHY_SPIKE_TIMES = []
ASSOCIATED_CALCIUM_TRACE = []

PCA = sklearn.decomposition.PCA()
KMEANS = sklearn.cluster.KMeans()

resampling_ = 4000
Fluo_resampling_ = 100

for i in range(len(ALL_ePhy_SIGNALS)):
    ePhy_SF = ALL_sampling_freqs_ephy[i]
    ePhy_LEN = len(ALL_ePhy_SIGNALS[i])*ePhy_SF
    
    temp_trace_ephy = sp.signal.resample(ALL_ePhy_SIGNALS[i], int(ePhy_LEN*resampling_))
    
    Derivated = [temp_trace_ephy[i+1]-temp_trace_ephy[i] for i in range(len(temp_trace_ephy)-6) if temp_trace_ephy[i]!=0]
    SpikePeakIndex = sp.signal.find_peaks(Derivated, height=2, distance=int(resampling_/4))[0]
    #SpikePeakIndex = Derivated.index(SpikePeak)
    
    plt.figure()
    plt.plot(temp_trace_ephy)
    for j in range(len(SpikePeakIndex)):
        plt.scatter(SpikePeakIndex[j], temp_trace_ephy[SpikePeakIndex[j]])
        
        temp_ = temp_trace_ephy[SpikePeakIndex[j]-20:SpikePeakIndex[j]+180]
        
        
        FullFluo = [SPIKE_LIST_MANUAL_CROP[k] for k in range(len(SPIKE_LIST_MANUAL_CROP_ID)) if int(SPIKE_LIST_MANUAL_CROP_ID[k].split('_0')[0].split('-')[-1])==ALL_TIMESTAMPS_ePhy[i]][0]
        FullFluoID = [SPIKE_LIST_MANUAL_CROP_ID[k] for k in range(len(SPIKE_LIST_MANUAL_CROP_ID)) if int(SPIKE_LIST_MANUAL_CROP_ID[k].split('_0')[0].split('-')[-1])==ALL_TIMESTAMPS_ePhy[i]][0]
        CameraDelay = [ALL_MicamDelay[k] for k in range(len(SPIKE_LIST_MANUAL_CROP_ID)) if int(SPIKE_LIST_MANUAL_CROP_ID[k].split('_0')[0].split('-')[-1])==ALL_TIMESTAMPS_ePhy[i]][0]
        FullFluoLen = [SPIKE_LIST_MANUAL_CROP_LEN[k] for k in range(len(SPIKE_LIST_MANUAL_CROP_ID)) if int(SPIKE_LIST_MANUAL_CROP_ID[k].split('_0')[0].split('-')[-1])==ALL_TIMESTAMPS_ePhy[i]][0]
        
        TraceSampling = len(FullFluo)/FullFluoLen
        FluoEventPeak = ((SpikePeakIndex[j]/resampling_)-1-CameraDelay)*TraceSampling
        
        FluoEvent_WindowStart = int(FluoEventPeak) - int(TraceSampling)
        FluoEvent_WindowEnd = int(FluoEventPeak) + int(TraceSampling*2)
        if FluoEvent_WindowStart>=0 and FluoEvent_WindowEnd<len(FullFluo):
            temp = sp.signal.resample(FullFluo[FluoEvent_WindowStart:FluoEvent_WindowEnd], Fluo_resampling_)
            temp = temp * 0.0014665418 * Fluo_resampling_
            ASSOCIATED_CALCIUM_TRACE.append(temp)
            DETECTED_ePHY_SPIKES_WF.append(temp_)
            DETECTED_ePHY_SPIKES_TIMESTAMP.append(ALL_TIMESTAMPS_ePhy[i])
            
            DETECTED_ePHY_SPIKE_TIMES.append(SpikePeakIndex[j])



plt.figure()

TO_PCA = []
for i in range(len(DETECTED_ePHY_SPIKE_TIMES)):
    
    temp_ = DETECTED_ePHY_SPIKES_WF[i]-np.nanmedian(DETECTED_ePHY_SPIKES_WF[i][0:int(Fluo_resampling_/5)])
    
    TO_PCA.append(temp_)
    plt.plot(temp_, lw=0.1, color='black')



PCA_fit = PCA.fit_transform(TO_PCA)
KMEAN_CLUSTER = KMEANS.fit_predict(PCA_fit)
X_fit = []
Y_fit = []

plt.figure(figsize=(15, 2))

ax = plt.subplot(171)
ax2 = plt.subplot(172)
ax3 = plt.subplot(173)
ax4 = plt.subplot(174)
ax5 = plt.subplot(175)
ax6 = plt.subplot(176)
ax7 = plt.subplot(177)
for i in range(np.nanmax(KMEAN_CLUSTER)):
    temp_ = [PCA_fit[:,0][j] for j in range(len(KMEAN_CLUSTER)) if KMEAN_CLUSTER[j]==i]
    temp_2 = [PCA_fit[:,1][j] for j in range(len(KMEAN_CLUSTER)) if KMEAN_CLUSTER[j]==i]
    temp_3 = [TO_PCA[j] for j in range(len(KMEAN_CLUSTER)) if KMEAN_CLUSTER[j]==i]
    temp_4 = [ASSOCIATED_CALCIUM_TRACE[j]-np.nanmedian(ASSOCIATED_CALCIUM_TRACE[j][0:25]) for j in range(len(KMEAN_CLUSTER)) if KMEAN_CLUSTER[j]==i]
    temp_5 = [np.nanmax(temp_4[j]) for j in range(len(temp_4))]
    temp_7 = [temp_4[j].tolist().index(np.nanmax(temp_4[j]))*(1/Fluo_resampling_) for j in range(len(temp_4))]
    temp_8 = [np.nanmedian(ASSOCIATED_CALCIUM_TRACE[j][0:int(Fluo_resampling_/4)]) for j in range(len(temp_4))]
    temp_6 = []
    for j in range(len(temp_4)):
        try:
            temp_6.append([k*(0.050/200) for k in range(20, len(temp_3[j])) if temp_3[j][k]>5][-1])
        except:
            temp_6.append(np.nan)
    X_fit.append([temp_5[k] for k in range(len(temp_5)) if temp_6[k]>0 and np.nanmedian(temp_3[k][20:25])>20])
    Y_fit.append([temp_6[k] for k in range(len(temp_5)) if temp_6[k]>0 and np.nanmedian(temp_3[k][20:25])>20])
    
    #CLUSTERED SPIKE TRACES
    ax.scatter(temp_, temp_2)
    MEAN = np.nanmean(temp_3, axis=0)
    SEM = sp.stats.sem(temp_3, axis=0)
    X = np.linspace(0, len(MEAN)/ resampling_, len(MEAN))
    ax2.plot(X, MEAN)
    ax2.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)
    
    #CORRESPONDING CALCIUM TRACES
    SEM = sp.stats.sem(temp_4, axis=0)
    MEAN = np.nanmean(temp_4, axis=0)
    X = np.linspace(0, len(MEAN)/ Fluo_resampling_, len(MEAN))
    ax3.plot(X, MEAN)
    ax3.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)    
    
    ax5.scatter(temp_7, temp_6 ,color='black', alpha=0.1)
    ax5.scatter(np.nanmean(temp_7), np.nanmean(temp_6))
    ax6.scatter(temp_5, temp_6 ,color='black', alpha=0.1)
    ax6.scatter(np.nanmean(temp_5), np.nanmean(temp_6))
    ax7.scatter(temp_8, temp_6 ,color='black', alpha=0.1)
    ax7.scatter(np.nanmean(temp_8), np.nanmean(temp_6))
    print(i, 'CaSpkAmp:',np.nanmean(temp_4), '+/-', sp.stats.sem(temp_4), 'ePhyHW:',np.nanmean(temp_6), '+/-', sp.stats.sem(temp_6))

X_fit = np.concatenate(X_fit)
Y_fit = np.concatenate(Y_fit)

popt, pcov = curve_fit(log_func , X_fit, Y_fit, maxfev=1000)
x_ = np.linspace(-20, np.nanmax(X_fit), len(X_fit))
ax6.plot(x_, log_func(x_, *popt), color='black', lw=0.1)

PCA_model = PCA.fit(TO_PCA)
ax4.plot(PCA_model.components_[0])
ax4.plot(PCA_model.components_[1])
ax4.plot(PCA_model.components_[2])


ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax2.set_xlabel('Time(s)')
ax2.set_ylabel('Mb.potential (mW)')
ax3.set_xlabel('Time(s)')
ax3.set_ylabel('GCamp6s fluo. (fW)')
ax4.set_xlabel('Time(s)')
ax4.set_ylabel('Principal components')
ax5.set_xlabel('Ca.Transient RT(ms)')
ax5.set_ylabel('Spike half-width (ms)')
ax6.set_xlabel('Ca.Transient amp. (fW)')
ax6.set_ylabel('Spike half-width (ms)')
ax7.set_xlabel('Ca.FZero (fW)')
ax7.set_ylabel('Spike half-width (ms)')
plt.tight_layout()
"""