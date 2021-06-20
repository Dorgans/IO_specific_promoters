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

ALL_STO_POWER = []
ALL_STO_POWER_f = []
ALL_STO_CCORRCOEFF = []
ALL_STO_POWER_PROMOTER = []
ALL_FZERO = []

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
    f, Pxx_spec = signal.welch(MEAN , SAMPLING_FREQUENCY , window='blackman', noverlap = SAMPLING_FREQUENCY/2 , nperseg=SAMPLING_FREQUENCY*4 , nfft=SAMPLING_FREQUENCY*8, average='mean')
    MEAN = MEAN - MED_FILT
    
    Max_Ccorr_peaks = []
    plt.figure()
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)
    ax1.plot(np.linspace(-1, 1, len(rs)), rs)
    for j in range(len(rs_peaks)):
        ax1.scatter(np.linspace(-1, 1, len(rs))[rs_peaks[j]], rs[rs_peaks[j]], color='red')
        Max_Ccorr_peaks.append(rs[rs_peaks[j]])
    ax2.plot(f, Pxx_spec)
    ax2.set_xlim(2, 20)
    ax2.set_ylim(0, 1000)

    MaxPwr = np.nanmax([Pxx_spec[j] for j in range(len(Pxx_spec)) if f[j]>=3])
    MaxPwr_f = f[Pxx_spec.tolist().index(MaxPwr)]
    ax2.scatter(MaxPwr_f, MaxPwr, color='red')
    
    if 3.3 < MaxPwr_f < 13:
        ALL_STO_POWER.append(MaxPwr)
        ALL_STO_POWER_f.append(MaxPwr_f)
        ALL_STO_CCORRCOEFF.append(np.nanmedian(Max_Ccorr_peaks))
        ALL_STO_POWER_PROMOTER.append(STO_LIST_MANUAL_CROP_PROMOTER[i])
        ALL_FZERO.append(np.nanmedian(STO_LIST_MANUAL_CROP[i][0:int(SAMPLING_FREQUENCY)])*0.0014665418*SAMPLING_FREQUENCY)
        
        
X_fit = []
Y_fit = []
plt.figure()
ax = plt.subplot(121)
ax2 = plt.subplot(122)
for i in range(len(LIST_)):
    temp_x = []
    temp_y = []
    for j in range(len(ALL_STO_POWER_PROMOTER)):
        if ALL_STO_POWER_PROMOTER[j] == LIST_[i]:
            temp_x.append(ALL_FZERO[j])
            temp_y.append(ALL_STO_POWER[j])
    X_fit.append(np.nanmean(temp_x))
    Y_fit.append(np.nanmean(temp_y))
    ax.scatter(np.nanmean(temp_x), np.nanmean(temp_y))
    ax2.scatter(np.nanmean(temp_x), np.nanmean(temp_y)/np.nanmean(temp_x))

popt, pcov = curve_fit(line_func , [X_fit[i] for i in range(len(X_fit)) if X_fit[i] <1000], [Y_fit[i] for i in range(len(Y_fit)) if X_fit[i] <1000], maxfev=1000)
x_ = np.linspace(0, np.nanmax(X_fit), 100)
ax.plot(x_, line_func(x_, *popt), color='black', lw=0.1)
