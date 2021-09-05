# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 17:18:11 2021

@author: KEVIN-DORGANS
"""

"""
#Ex of importation method from Fiji
STO_LIST_MANUAL_CROP.append(np.array(pd.read_clipboard())[:,1])
STO_LIST_MANUAL_CROP_ID.append('5HTr2b(1.8)-GCamp6s_IO_2W_200nlx2_1to2.5_LED40_30fps2021-01-26-200751_00_N256 (IF1-CAM1)')
STO_LIST_MANUAL_CROP_PROMOTER.append('AAV9.5HTr2b(1.8)')

#Array_to_csv
PATH = load_directory_content__()[1]
pd.DataFrame(STO_LIST_MANUAL_CROP).to_csv(PATH+'\FijiSomaROIcropped_STO_PER_PROMOTER_rawfluo_db.csv')
pd.DataFrame(STO_LIST_MANUAL_CROP_PROMOTER).to_csv(PATH+'\FijiSomaROIcropped_STO_PER_PROMOTER_rawfluo_PromNum_db.csv')
pd.DataFrame(STO_LIST_MANUAL_CROP_ID).to_csv(PATH+'\FijiSomaROIcropped_STO_PER_PROMOTER_rawfluo_File_db.csv')

#Csv_to_array
PATH = load_directory_content__()[1]
STO_LIST_MANUAL_CROP = []
try:
    temp_ = pd.read_csv(PATH+'\FijiSomaROIcropped_STO_PER_PROMOTER_rawfluo_db.csv', index_col=0, header=0)
    temp_ = np.array(temp_.values)
    for i in range(len(temp_)):
        STO_LIST_MANUAL_CROP.append([temp_[i][j] for j in range(len(temp_[i])) if str(temp_[i][j]) != 'nan'])
    temp_ = pd.read_csv(PATH+'\FijiSomaROIcropped_STO_PER_PROMOTER_rawfluo_PromNum_db.csv')
    STO_LIST_MANUAL_CROP_PROMOTER = np.array(temp_.values[:,1])
    temp_ = pd.read_csv(PATH+'\FijiSomaROIcropped_STO_PER_PROMOTER_rawfluo_File_db.csv')
    STO_LIST_MANUAL_CROP_ID = np.array(temp_.values[:,1])
except:
    print('file not recognized')

"""
def log_func(x, a, b, c):
    return a * np.log2(-b * x) + c

def line_func(x, a, b):
    return a * x + b

def exp_func(x, a, b, c):
    return a * np.exp(b * x) + c

ALL_STO_TRACES = []
ALL_STO_TRACES_ID = []
ALL_STO_POWER = []
ALL_STO_POWER_f = []
ALL_STO_CCORRCOEFF = []
ALL_STO_POWER_PROMOTER = []
ALL_FZERO = []

LIST_ = []
for i in range(len(STO_LIST_MANUAL_CROP_PROMOTER)):
    if (STO_LIST_MANUAL_CROP_PROMOTER[i] in LIST_)==False:
        LIST_.append(STO_LIST_MANUAL_CROP_PROMOTER[i])

for i in range(len(STO_LIST_MANUAL_CROP)):
    
    MEAN = np.array(STO_LIST_MANUAL_CROP[i], dtype=np.float64)
    SAMPLING_FREQUENCY = len(MEAN)/10
    MEAN = MEAN*0.0014665418*SAMPLING_FREQUENCY

    try:
        MED_FILT = sp.signal.medfilt(MEAN, np.int(0.5*SAMPLING_FREQUENCY))
    except:
        MED_FILT = sp.signal.medfilt(MEAN, np.int(0.5*SAMPLING_FREQUENCY)-1)
    
    df = pd.DataFrame()
    df['1'] = MEAN - MED_FILT
    df['2'] = MEAN - MED_FILT
    rs = [crosscorr(df['1'], df['1'], lag) for lag in range(-int(SAMPLING_FREQUENCY),int(SAMPLING_FREQUENCY+1))]
    rs_peaks = sp.signal.find_peaks(rs, distance=SAMPLING_FREQUENCY/10)[0]
    f, Pxx_spec = signal.welch(MEAN , SAMPLING_FREQUENCY , window='blackman', noverlap = SAMPLING_FREQUENCY/2 , nperseg=SAMPLING_FREQUENCY*7 , nfft=SAMPLING_FREQUENCY*18, average='mean')
    MEAN = MEAN - MED_FILT
    
    Max_Ccorr_peaks = []
    plt.figure(num=STO_LIST_MANUAL_CROP_ID[i].split('_00')[0].split('-')[-1])
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)
    ax1.plot(np.linspace(-1, 1, len(rs)), rs)
    for j in range(len(rs_peaks)):
        ax1.scatter(np.linspace(-1, 1, len(rs))[rs_peaks[j]], rs[rs_peaks[j]], color='red')
        Max_Ccorr_peaks.append(rs[rs_peaks[j]])
    ax2.plot(f, Pxx_spec)
    ax2.set_xlim(2, 20)
    ax2.set_ylim(0, 1000)

    MaxPwr = np.nanmax([Pxx_spec[j] for j in range(len(Pxx_spec)) if f[j]>=2])
    MaxPwr_f = f[Pxx_spec.tolist().index(MaxPwr)]
    ax2.scatter(MaxPwr_f, MaxPwr, color='red')
    
    if 3.5 < MaxPwr_f < 13:
        ALL_STO_POWER.append(MaxPwr)
        ALL_STO_POWER_f.append(MaxPwr_f)
        ALL_STO_CCORRCOEFF.append(np.nanmedian(Max_Ccorr_peaks))
        ALL_STO_POWER_PROMOTER.append(STO_LIST_MANUAL_CROP_PROMOTER[i])
        ALL_FZERO.append(np.nanmedian(STO_LIST_MANUAL_CROP[i][0:int(SAMPLING_FREQUENCY)])*0.0014665418*SAMPLING_FREQUENCY)
        ALL_STO_TRACES.append(STO_LIST_MANUAL_CROP[i])
        ALL_STO_TRACES_ID.append(int(STO_LIST_MANUAL_CROP_ID[i].split('_00')[0].split('-')[-1]))

X_fit = []
Y_fit = []
plt.figure(figsize=(13,3))
ax = plt.subplot(161)
ax2 = plt.subplot(162)
ax3 = plt.subplot(163)
ax4 = plt.subplot(164)
ax5 = plt.subplot(165)
ax6 = plt.subplot(166)

for i in range(len(LIST_)):
    temp_x = []
    temp_y = []
    for j in range(len(ALL_STO_POWER_PROMOTER)):
        if ALL_STO_POWER_PROMOTER[j] == LIST_[i]:
            temp_x.append(ALL_FZERO[j])
            temp_y.append(ALL_STO_POWER[j])
    print(LIST_[i], np.nanmean(temp_y), '+/-', sp.stats.sem(temp_y), 'n=', len(temp_y))
    X_fit.append(np.nanmean(temp_x))
    Y_fit.append(np.nanmean(temp_y))
    MEAN = np.nanmean(temp_x)
    MEAN_2 = np.nanmean(temp_y)
    SEM = sp.stats.sem(temp_x)
    SEM_2 = sp.stats.sem(temp_y)
    SEM_3 = sp.stats.sem(np.divide(temp_y, temp_x))
    
    ax.scatter(MEAN, MEAN_2, label=LIST_[i])
    ax2.scatter(MEAN, MEAN_2/MEAN)
    ax.plot((MEAN-SEM, MEAN+SEM), (MEAN_2, MEAN_2), color='black')
    ax.plot((MEAN, MEAN), (MEAN_2-SEM_2, MEAN_2+SEM_2), color='black')
    ax2.plot((MEAN, MEAN), (MEAN_2/MEAN-SEM_3, MEAN_2/MEAN+SEM_3), color='black')
    ax2.plot((MEAN-SEM, MEAN+SEM), (MEAN_2/MEAN, MEAN_2/MEAN), color='black')
    
    ax5.scatter(MEAN, MEAN_2, label=LIST_[i])
    ax6.scatter(MEAN, MEAN_2/MEAN)
    ax5.plot((MEAN-SEM, MEAN+SEM), (MEAN_2, MEAN_2), color='black')
    ax5.plot((MEAN, MEAN), (MEAN_2-SEM_2, MEAN_2+SEM_2), color='black')
    ax6.plot((MEAN, MEAN), (MEAN_2/MEAN-SEM_3, MEAN_2/MEAN+SEM_3), color='black')
    ax6.plot((MEAN-SEM, MEAN+SEM), (MEAN_2/MEAN, MEAN_2/MEAN), color='black')
     
    
    ax3.hist(temp_y, histtype='step', cumulative=True, density=True, bins=200)
    ax4.hist(np.multiply(np.divide(temp_y, temp_x), 100), histtype='step', cumulative=True, density=True, bins=200)
    ax5.scatter(temp_x, temp_y, s=5)
    ax6.scatter(temp_x, np.divide(temp_y, temp_x), s=5)
popt, pcov = curve_fit(line_func , [X_fit[i] for i in range(len(X_fit)) if X_fit[i] <1000], [Y_fit[i] for i in range(len(Y_fit)) if X_fit[i] <1000], maxfev=1000)
x_ = np.linspace(0, np.nanmax(X_fit), 100)
ax.plot(x_, line_func(x_, *popt), color='black', lw=0.1)
ax.set_xlabel('F_ZERO(fW)')
ax.set_ylabel('STO Power(fW**2/Hz)')
ax2.set_xlabel('F_ZERO(fW)')
ax2.set_ylabel('STO Power(fW**2/Hz)/F_ZERO(fW)')
ax.legend()
ax3.set_xscale('log')
ax4.set_xscale('log')
ax5.set_xscale('log')
ax6.set_xscale('log')
ax5.set_yscale('log')
ax6.set_yscale('log')
ax3.set_xlabel('STO Power(fW**2/Hz)')
ax4.set_xlabel('STO Power(fW**2/Hz)/F_ZERO(fW)')
ax5.set_xlabel('F_ZERO(fW)')
ax6.set_xlabel('F_ZERO(fW)')
ax5.set_ylabel('STO Power(fW**2/Hz)')
ax6.set_ylabel('STO Power(fW**2/Hz)/F_ZERO(fW)')
plt.tight_layout()


"""
#SELECTION FOR PHASE WRAP BASED OF PSD FROM HIGH-PASSED FILTERED
#1 > CALCULATE PSD AGAIN FROM RAW FLUO(bits) AND HIGH PASS FILTER WITH NOISE-PSD (=1bit2)

FULL_WRAPPED_WAVES = []
FULL_WRAPPED_WAVES_PROMOTER = []
ALL_STO_POWER_f_HIGH_PASSED = []
ALL_STO_POWER_f_HIGH_PASSED_PSD = []
ALL_STO_TRACES_HIGH_PASSED = []
ALL_STO_POWER_PROMOTER_HIGH_PASSED = []
ALL_FZERO_HIGH_PASSED = []
ALL_STO_TRACES_HIGH_PASSED_ID = []

freq=1
for j in range(len(ALL_STO_TRACES)):
    ALL = []
    offset =0
    window = int(SAMPLING_FREQUENCY/(freq))
    freq = ALL_STO_POWER_f[j]
    DATA = ALL_STO_TRACES[j]
    SAMPLING_FREQUENCY = len(DATA)/10
    
    X = np.linspace(0, len(DATA), len(DATA))
    try:
        popt, pcov = curve_fit(line_func , X, sp.signal.medfilt(DATA, int(SAMPLING_FREQUENCY)), maxfev=1000)
    except:
        popt, pcov = curve_fit(line_func , X, sp.signal.medfilt(DATA, int(SAMPLING_FREQUENCY)+1), maxfev=1000)
    DATA = DATA - line_func(X, *popt)

    f, Pxx_spec = signal.welch(DATA , SAMPLING_FREQUENCY , window='blackman', noverlap = SAMPLING_FREQUENCY/2 , nperseg=SAMPLING_FREQUENCY*7 , nfft=SAMPLING_FREQUENCY*10, average='median', scaling='spectrum')
    plt.figure(num = (str(ALL_STO_TRACES_ID[j])+' '+ALL_STO_POWER_PROMOTER[j]))
    ax = plt.subplot(121)
    ax2 = plt.subplot(122)
    ax.plot(DATA)
    ax2.plot(f, Pxx_spec)
    Pxx_spec_band_passed = [Pxx_spec[i] for i in range(len(Pxx_spec)) if 3<f[i]<13]
    f_band_passed = [f[i] for i in range(len(Pxx_spec)) if 3<f[i]<13]
    MaxPSD = np.nanmax(Pxx_spec_band_passed)
    MaxFreq = f_band_passed[Pxx_spec_band_passed.index(MaxPSD)]
    
    #MINIMAL VALUE OBSERVED FOR STO_PSD_IMAGING 3-13Hz if ePhy STO>1mV = 0.00024193322287634385
    #MEAN+3*SD STO_PSD_IMAGING 3-13Hz if ePhy STO<1mV =0.0012017975487199199:

    if (MaxPSD*0.0014665418*SAMPLING_FREQUENCY)/ALL_FZERO[j]< 0.00024193322287634385:
        ax2.scatter(MaxFreq, MaxPSD, color='red')
    else:
        ax2.scatter(MaxFreq, MaxPSD, color='blue')
        
        ALL_STO_POWER_f_HIGH_PASSED.append(MaxFreq)
        ALL_STO_POWER_f_HIGH_PASSED_PSD.append(MaxPSD*0.0014665418*SAMPLING_FREQUENCY)
        ALL_STO_TRACES_HIGH_PASSED.append((np.array(ALL_STO_TRACES[j])*0.0014665418*SAMPLING_FREQUENCY).tolist())
        ALL_STO_POWER_PROMOTER_HIGH_PASSED.append(ALL_STO_POWER_PROMOTER[j]) 
        ALL_FZERO_HIGH_PASSED.append(ALL_FZERO[j])
        ALL_STO_TRACES_HIGH_PASSED_ID.append(ALL_STO_TRACES_ID[j])

#2 ISOLATE THE PHASES WITH EVEN SIGNAL CUT FROM ACCURATE PSD PEAK (FROM#1)
#SHIFTS PHASE PEAKS TO MID-PHASE
        
FULL_WRAPPED_WAVES = []
FULL_WRAPPED_WAVES_PROMOTER = []
FULL_WRAPPED_WAVES_FZERO = []

for j in range(len(ALL_STO_TRACES_HIGH_PASSED)):
    ALL = []
    offset =0
    
    freq=1
    DATA = ALL_STO_TRACES_HIGH_PASSED[j]
    SAMPLING_FREQUENCY = len(DATA)/10
    freq = ALL_STO_POWER_f_HIGH_PASSED[j]
    window = int(SAMPLING_FREQUENCY/(freq))
   
    
    
    
    
    for i in range(int(len(DATA)/int(SAMPLING_FREQUENCY/freq))):
        if True:
            temp_ = DATA[offset+ i*window: offset + (i+1)*window]
            if len(temp_)>0 and 0.47 < temp_.index(np.nanmax(temp_))/len(temp_) < 0.53:
                if len(ALL)==0 or len(temp_)==len(ALL[-1]):
                    ALL.append(temp_-np.nanmean(temp_))
            elif len(temp_)>0:
                    new_offset  = int(offset + (temp_.index(np.nanmax(temp_))/len(temp_)-0.5)*len(temp_) )
                    temp_ = DATA[new_offset + i*window: new_offset + (i+1)*window]
                    
            if len(temp_) == window :
                ALL.append(temp_-np.nanmedian(temp_))
                    
        
    MEAN = np.nanmean(ALL, axis=0)
    SEM = sp.stats.sem(ALL, axis=0)
    X = np.linspace(0, 1, len(MEAN))
    
    if len(ALL)>0:
        plt.figure(figsize=(2,2), num=ALL_STO_TRACES_HIGH_PASSED_ID[j])
        for l in range(len(ALL)):
            plt.plot(X, ALL[l], color='black', lw=0.1)
        plt.plot(X, MEAN)
        plt.scatter(X, MEAN)
        plt.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)
        
        if True:
            FULL_WRAPPED_WAVES.append(sp.signal.resample(MEAN, 30))
            FULL_WRAPPED_WAVES_PROMOTER.append(ALL_STO_POWER_PROMOTER_HIGH_PASSED[j])
            FULL_WRAPPED_WAVES_FZERO.append(ALL_FZERO_HIGH_PASSED[j])
            
            
#3 GETS PHASE-SHIFT ERRORS AWAY AND RESAMPLE WAVES
FILTERED_FZEROS_PER_PROMOTER = []
FILTERED_DELTAFWAVE_PER_PROMOTER = []

plt.figure(figsize=(15, 4))
for j in range(len(LIST_)):
    ax = plt.subplot(2, len(LIST_), j+1)
    ax2 = plt.subplot(2, len(LIST_), len(LIST_)+j+1)
    temp_ = []
    temp_fzero = []
    
    for i in range(len(FULL_WRAPPED_WAVES_PROMOTER)):
        WrappedWavePeaks = sp.signal.find_peaks(FULL_WRAPPED_WAVES[i], threshold=-1)[0]
        if FULL_WRAPPED_WAVES_PROMOTER[i]==LIST_[j] and len(WrappedWavePeaks)==1:
            if 0.3 < WrappedWavePeaks[0]/len(FULL_WRAPPED_WAVES[i]) < 0.6:
                temp_.append(FULL_WRAPPED_WAVES[i])
                temp_fzero.append(FULL_WRAPPED_WAVES_FZERO[i])
    MEAN = np.nanmean(temp_, axis=0)
    
    print(sp.signal.find_peaks(MEAN)[0])
    
    X = np.linspace(0, 1, len(MEAN))
    SEM =sp.stats.sem(temp_, axis=0)
    ax.plot(X, MEAN)
    
    ax2. plot(np.linspace(0, 4, 4*len(MEAN)), np.tile(MEAN, 4))
    
    ax.fill_between(X, MEAN+SEM, MEAN-SEM, alpha=0.1)
    #for i in range(len(temp_)):
        #ax.plot(X, temp_[i], lw=0.1, color='black')
    print(LIST_[j], 'n=', len(temp_), 'max=', np.nanmax(MEAN), 'fzero=', np.nanmean(temp_fzero))
    FILTERED_FZEROS_PER_PROMOTER.append(np.nanmean(temp_fzero))
    FILTERED_DELTAFWAVE_PER_PROMOTER.append(np.nanmax(MEAN)-np.nanmin(MEAN))


    ax.set_title(LIST_[j])
    ax.set_xlabel('Normalized Phase')
    ax.set_ylabel('DeltaF (fW)')
    ax2.set_xlabel('4 tiled phases')
    ax2.set_ylabel('DeltaF (fW)')
plt.tight_layout()

"""


"""
#SOME SCRIPT THAT SUPERPOSES STO from IMAGING AND ePHY and Specgram
#IT WORKS WITH "SPIKE-SCRIPT"


imaging_traces = np.concatenate(ALL_SPIKE_WAVES_PROM)
PSD_imaging = []
PSD_ePhy = []
FREQ_imaging = []
FREQ_ePhy = []
ImagingFZERO = []

for i in range(len(SPIKE_LIST_MANUAL_CROP)):
    
    imaging_timestamp = int(SPIKE_LIST_MANUAL_CROP_ID[i].split('_00')[0].split('-')[-1])
    imaging_trace = SPIKE_LIST_MANUAL_CROP[i]
    imaging_length = SPIKE_LIST_MANUAL_CROP_LEN[i]
    imaging_sampling = imaging_length/len(imaging_trace)
    imaging_trace = np.array(imaging_trace)
    
    
    ePhy_trace_id = ALL_TIMESTAMPS_ePhy.index(imaging_timestamp)
    ePhy_sampling = ALL_sampling_freqs_ephy[ePhy_trace_id]
    ePhy_trace = ALL_ePhy_SIGNALS[ePhy_trace_id]
    imaging_delay = ALL_MicamDelay[ePhy_trace_id] + 1
    SAMPLING_FREQUENCY = int(1/ePhy_sampling)
    
    plt.figure(figsize=(5,2), num=imaging_timestamp)
    
    ax = plt.subplot(141)
    ax2 = plt.subplot(142)
    ax3 = plt.subplot(143)
    ax4 = plt.subplot(144)
    ax4.plot(ePhy_trace)
    X = np.linspace(0, imaging_length, len(imaging_trace))
    MEAN = imaging_trace - np.nanmedian(imaging_trace)
    MEAN = MEAN
    #ax.plot(X,  MEAN, alpha=0.1)
    
    try:
        popt, pcov = curve_fit(line_func , X, sp.signal.medfilt(MEAN, int(imaging_sampling)), maxfev=1000)
    except:
        popt, pcov = curve_fit(line_func , X, sp.signal.medfilt(MEAN, int(imaging_sampling)+1), maxfev=1000)
    #ax.plot(X, line_func(X, *popt), color='black', lw=0.1)
    MEAN = MEAN - line_func(X, *popt)
    
    MEAN_IMAGING = MEAN
    #ax.plot(X,  MEAN, color='black')
    
    MEAN = ePhy_trace - np.nanmedian(ePhy_trace)
    MEAN = MEAN
    MEAN = MEAN[int(imaging_delay/ePhy_sampling): int((imaging_length+imaging_delay)/ePhy_sampling)]
    X = np.linspace(0, imaging_length , len(MEAN))
    #ax.plot(X,  MEAN)
    
    PEAKS = sp.signal.find_peaks(MEAN, height=0.01, distance=int((1/ePhy_sampling)/1))[0]
    
    f, Pxx_spec = signal.welch(MEAN , SAMPLING_FREQUENCY , window='bartlett', noverlap = SAMPLING_FREQUENCY/4, nperseg=SAMPLING_FREQUENCY , nfft=SAMPLING_FREQUENCY/0.5, average='mean', scaling='spectrum')
    MaxPwr = np.nanmax([Pxx_spec[j] for j in range(len(Pxx_spec)) if f[j]>=3.5])
    MaxPwr_f = f[Pxx_spec.tolist().index(MaxPwr)]
    MaxePhySTO = Pxx_spec.tolist().index(MaxPwr)
    ax.plot(f, Pxx_spec)
    ax.scatter(MaxPwr_f, MaxPwr, color='red')
    ax3.bar(0, MaxPwr)
    f, Pxx_spec = signal.welch(MEAN_IMAGING, int(1/imaging_sampling) , window='blackman', noverlap = int(1/imaging_sampling)/4 , nperseg=int(1/imaging_sampling)/1 , nfft=int(1/imaging_sampling)/0.5, average='median', scaling='spectrum')
    ax2.plot(f, Pxx_spec)
    
    MaxPwr_IMAGING = np.nanmax([Pxx_spec[j] for j in range(len(Pxx_spec)) if 13>f[j]>3.5])
    MaxPwr_f_IMAGING = f[Pxx_spec.tolist().index(MaxPwr_IMAGING)]
    ax3.bar(1, MaxPwr_IMAGING*0.0014665418*(1/imaging_sampling))
    ax2.scatter(MaxPwr_f_IMAGING , MaxPwr_IMAGING, color='red')
    if 3.6 < MaxPwr_f_IMAGING < 13:
        PSD_ePhy.append(MaxPwr)
        PSD_imaging.append(MaxPwr_IMAGING*0.0014665418*(1/imaging_sampling))
        ImagingFZERO.append(np.nanmedian(SPIKE_LIST_MANUAL_CROP[i])*0.0014665418*(1/imaging_sampling))
        FREQ_imaging.append(MaxPwr_f_IMAGING)
        FREQ_ePhy.append(MaxPwr_f)
        
plt.figure()
ax = plt.subplot(121)
ax2 = plt.subplot(122)
ax.scatter(PSD_ePhy, np.divide(PSD_imaging, ImagingFZERO))
DIV = np.divide(PSD_imaging, ImagingFZERO)

SUP = [DIV[i] for i in range(len(DIV)) if PSD_ePhy[i]>1]
INF = [DIV[i] for i in range(len(DIV)) if PSD_ePhy[i]<1]
ax.set_xscale('log')
ax.set_xlabel('ePhy-STO(mV**2 RMS)')
ax.set_ylabel('GCamp6s-STO (DFF)')

ax2.boxplot([SUP, INF], labels=['ePhy-STO>1mV','ePhy-STO<1mV'])
ax2.set_ylabel('GCamp6s-STO DFF')

"""


"""
#STO DFF vs SPIKE DFF

plt.figure(figsize=(4,4))
for i in range(len(LIST_)):
    MEAN = np.nanmean(np.divide(ALL_DELTA_SPIKE_PROM[i],ALL_F_ZERO_PROM[i]))
    SEM = sp.stats.sem(np.divide(ALL_DELTA_SPIKE_PROM[i],ALL_F_ZERO_PROM[i]))
    MEAN_2 = np.nanmean(np.divide(Y_fit[i], X_fit[i]))
    SEM_2 = sp.stats.sem(np.divide(Y_fit[i], X_fit[i]))
    plt.scatter(MEAN_2, MEAN, label=LIST_[i])
    plt.plot((MEAN_2, MEAN_2), (MEAN + SEM, MEAN-SEM), color='black')
    plt.plot((MEAN_2+SEM_2, MEAN_2-SEM_2), (MEAN, MEAN), color='black')
plt.legend()
plt.xlabel('STO DFF')
plt.ylabel('SPIKE_DFF')
plt.tight_layout()
"""