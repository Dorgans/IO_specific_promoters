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

"""

def log_func(x, a, b, c):
    return a * np.log2(b * x) + c

def line_func(x, a, b):
    return a * x + b


LED_INTENSITY_ = []
PROMOTER_NAME_ = []
FILE_NAME = []
ALL_SPIKES = []

for l in range(len(SPIKE_LIST_MANUAL_CROP)):
    
    if True:
        PATH = SPIKE_LIST_MANUAL_CROP_ID[l]
        try:
            
            
            
            try:
                MEAN = SPIKE_LIST_MANUAL_CROP[l]
                SAMPLING_FREQUENCY = len(MEAN)/10
                
                FILT = [MEAN[i+3]/MEAN[i] for i in range(len(MEAN)-3)]
                FILT = FILT - np.nanmean(FILT)
                FILT_ = FILT
                ZSCORE_FILT = sp.stats.zscore(FILT_)
                PEAKS = sp.signal.find_peaks(ZSCORE_FILT, height=3, distance= SAMPLING_FREQUENCY/4)
                FILT = MEAN
                plt.figure()
                plt.title(SPIKE_LIST_MANUAL_CROP_ID[l])
                ax=plt.subplot(131)
                ax2 = plt.subplot(132)
                ax3 = plt.subplot(133)
                
                ax3.plot(sp.stats.zscore(FILT_))
                
                PEAK_LIST = []
                ax.plot(MEAN)
                for i in range(len(PEAKS[0])):
                    ax.scatter(PEAKS[0][i], FILT[PEAKS[0][i]], color='red')
                    try:
                        if PEAKS[0][i]+SAMPLING_FREQUENCY*2.5<len(FILT) and PEAKS[0][i]-SAMPLING_FREQUENCY/2>0:
                            
                            temp_ = FILT[int(PEAKS[0][i]- SAMPLING_FREQUENCY/2): int(PEAKS[0][i]+ SAMPLING_FREQUENCY*2.5)]
                            PEAK_LIST.append(temp_)
                        elif int(PEAKS[0][i]- SAMPLING_FREQUENCY/1)<0:
                            missing = abs(int(PEAKS[0][i]- SAMPLING_FREQUENCY/2))
                            temp_ = np.concatenate([np.nan for l in range(missing)], temp_)
                            PEAK_LIST.append(temp_)
                    except:
                        pass
                AVG = np.nanmean(PEAK_LIST, axis=0)
                SEM = sp.stats.sem(PEAK_LIST, axis=0)
                ax2.plot(np.linspace(0, 10, len(AVG)), AVG)
                ax2.fill_between(np.linspace(0, 10, len(AVG)), AVG+SEM, AVG-SEM, alpha=0.1)

                try:
                    LED_INTENSITY_.append(np.int(PATH.split('LED')[-1].split('_')[0]))
                except:
                    LED_INTENSITY_.append(np.nan)
                for k in range(len(PEAK_LIST)):
                    
                    PROMOTER_NAME_.append(SPIKE_LIST_MANUAL_CROP_PROMOTER[l])
                    FILE_NAME.append(SPIKE_LIST_MANUAL_CROP_ID[l])
                    ALL_SPIKES.append(sp.signal.resample(PEAK_LIST[k], 100)*0.0014665418*SAMPLING_FREQUENCY)
            except:
                pass
        except:
            pass
            

plt.figure()
ax = plt.subplot(121)
ax2 = plt.subplot(122)

n_ALL_SPIKES = []
for i in range(len(ALL_SPIKES)):
    n_ALL_SPIKES.append(ALL_SPIKES[i]-np.nanmin(ALL_SPIKES[i][0:int(SAMPLING_FREQUENCY/2)]))

SEM = sp.stats.sem(ALL_SPIKES, axis=0)
MEAN = np.nanmean(ALL_SPIKES, axis=0)
X = np.linspace(0, 10, len(MEAN))
ax.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)
ax.plot(X, np.nanmean(ALL_SPIKES, axis=0))

SEM = sp.stats.sem(n_ALL_SPIKES, axis=0)
MEAN = np.nanmean(n_ALL_SPIKES, axis=0)
X = np.linspace(0, 10, len(MEAN))
ax2.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)
ax2.plot(X, MEAN)
for i in range(len(n_ALL_SPIKES)):
    ax2.plot(X, n_ALL_SPIKES[i], color='black', lw=0.2)
    


LIST_ = ['AAV.PHP.eB.5HTr2b-tTA', 'AAV.PHP.S.5HTr2b-tTA', 'AAV9.5HTr2b-tTA','AAV9.5HTr2b(1.8)','AAV9.SUSD4(2.4)', 'AAV9.CAG']
X = np.linspace(0, 2.5, len(MEAN))

ALL_F_ZERO_PROM = []
ALL_DELTA_SPIKE_PROM = []

plt.figure()
for i in range(len(LIST_)):
    ALL_F_ZERO = []
    ALL_DELTA_SPIKE = []
    ALL_SPIKES_PER_PROMOTER = []
    
    AVG = [ALL_SPIKES[j] for j in range(len(ALL_SPIKES)) if PROMOTER_NAME_[j] == LIST_[i]]
    F_ZERO = [np.nanmin(AVG[j][5:18]) for j in range(len(AVG))]
    DELTA_MAX_SPIKE = [np.nanmax(AVG[j][10:30])-F_ZERO[j] for j in range(len(AVG))]
    MEAN = np.nanmean(AVG, axis=0)
    SEM = sp.stats.sem(AVG, axis=0)
    X = np.linspace(0, 10, len(MEAN))
    #plt.plot(X, MEAN, label=LIST_[i])
    #plt.fill_between(X, MEAN+SEM, MEAN-SEM, color='black', alpha=0.01)
    print(LIST_[i],'n=',len(AVG), 'max=', np.nanmax(MEAN)-MEAN[10])
    
    for k in range(len(AVG)):
        MEAN = AVG[k]
        """
        y_data = []
        x_data = []
        for l in range(len(AVG[k])):
            if l<40 or l>90:
                y_data.append(AVG[k][l])
                x_data.append(l)
        popt, pcov = curve_fit(line_func, x_data, y_data, maxfev=1000)
        LINE_FIT = line_func(np.linspace(0, len(AVG[k]), len(AVG[k])), *popt)
        MEAN = AVG[k]-LINE_FIT
        """
        MEAN = MEAN-np.nanmin(MEAN[0:20])
        if True:
            #plt.plot(X, MEAN, color='black', lw=0.1)
            ALL_F_ZERO.append(F_ZERO)
            ALL_DELTA_SPIKE.append(DELTA_MAX_SPIKE)
            ALL_SPIKES_PER_PROMOTER.append(MEAN)
    
    MEAN = np.nanmean(ALL_SPIKES_PER_PROMOTER, axis=0)
    SEM = sp.stats.sem(ALL_SPIKES_PER_PROMOTER, axis=0)
    
    plt.plot(X, np.nanmean(ALL_SPIKES_PER_PROMOTER, axis=0), label=LIST_[i])
    plt.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)
    ALL_F_ZERO_PROM.append(ALL_F_ZERO)
    ALL_DELTA_SPIKE_PROM.append(ALL_DELTA_SPIKE)
    
plt.xlabel('Time(s.)')
plt.legend()


plt.figure()

ax = plt.subplot(121)
ax2 = plt.subplot(122)

X_fit = []
Y_fit = []

for i in range(len(ALL_DELTA_SPIKE_PROM)):
    ax.scatter(np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(ALL_DELTA_SPIKE_PROM[i]))
    SEM = sp.stats.sem(ALL_DELTA_SPIKE_PROM[i][0])
    MEAN = np.nanmean(ALL_DELTA_SPIKE_PROM[i][0])
    ax.plot((np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(ALL_F_ZERO_PROM[i])), (MEAN+SEM, MEAN-SEM), color='black', lw=0.1)
    SEM = sp.stats.sem(ALL_F_ZERO_PROM[i][0])
    MEAN = np.nanmean(ALL_F_ZERO_PROM[i][0])
    ax.plot((MEAN+SEM, MEAN-SEM), (np.nanmean(ALL_DELTA_SPIKE_PROM[i]), np.nanmean(ALL_DELTA_SPIKE_PROM[i])), color='black', lw=0.1)


    X_fit.append(np.nanmean(ALL_F_ZERO_PROM[i]))
    Y_fit.append(np.nanmean(ALL_DELTA_SPIKE_PROM[i]))
    
popt, pcov = curve_fit(log_func , X_fit, Y_fit, maxfev=1000)
x_ = np.linspace(0, np.nanmax(X_fit), 100)
ax.plot(x_, log_func(x_, *popt), color='black', lw=0.1)

popt, pcov = curve_fit(line_func , [X_fit[i] for i in range(len(X_fit)) if X_fit[i] <1000], [Y_fit[i] for i in range(len(Y_fit)) if X_fit[i] <1000], maxfev=1000)
x_ = np.linspace(0, np.nanmax(X_fit), 100)
ax.plot(x_, line_func(x_, *popt), color='black', lw=0.1)

for i in range(len(ALL_DELTA_SPIKE_PROM)):
    ax2.scatter(np.nanmean(ALL_F_ZERO_PROM[i]), np.nanmean(np.divide(ALL_DELTA_SPIKE_PROM[i],ALL_F_ZERO_PROM[i])), label=LIST_[i])

ax.set_xlabel('Somatic F-Zero (fW)')
ax.set_ylabel('Calcium Spike Amplitude (fW)')
ax2.set_ylabel('Calcium Spike DF/F0')
ax2.set_xlabel('Somatic F-Zero (fW)')
plt.legend()
ax.set_xlim(0,)
ax.set_ylim(0,)