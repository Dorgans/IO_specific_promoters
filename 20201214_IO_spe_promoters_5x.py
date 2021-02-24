# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 21:39:19 2021

@author: dorga
"""
LIST = ['2.5kb',
 '5HTr2b',
 'PDX1',
 'SUSD4',
 '2.5kb',
 '1.3kb',
 '5HTr2b(1.0kb)',
 '5HTr2b(1.8kb)',
 'Igsf9(3.7kb)',
 '5HTr2b(3.0kb)',
 'AV9']
pix_to_um_conversion_factor = 0.3613

plt.figure()
for i in LIST:
    X = []
    Y = []
    
    for j in range(len(FULL_PROMOTER_NAMES)):    
        if i==FULL_PROMOTER_NAMES[j]:
            Y.append(np.divide(FULL_eGFP_SPREAD_IN[j], FULL_eGFP_SPREAD_OUT[j])/np.divide(FULL_tDt_SPREAD_IN[j], FULL_tDt_SPREAD_OUT[j]))
            X.append(np.divide(FULL_eGFP_IO_IN[j], FULL_eGFP_IO_OUT[j])/np.divide(FULL_tDt_IO_IN[j], FULL_tDt_IO_OUT[j]))
    print(str(i)+' : '+str(len(X)))
    plt.scatter(np.nanmean(X), np.nanmean(Y))
    plt.text(np.nanmean(X), np.nanmean(Y), i)
    MEAN = np.nanmean(X)
    SEM = sp.stats.sem(X)
    plt.plot((MEAN+SEM, MEAN-SEM), (np.nanmean(Y),np.nanmean(Y)), color='black')
    MEAN = np.nanmean(Y)
    SEM = sp.stats.sem(Y)
    plt.plot((np.nanmean(X),np.nanmean(X)), (MEAN+SEM, MEAN-SEM),  color='black')
plt.ylabel('IO-SPECIFICITY')
plt.xlabel('EXPRESSION LEVEL IN IO')

plt.figure()
for i in LIST:
    X = []
    Y = []
    X_tDt = []
    Y_tDt = []
    for j in range(len(FULL_PROMOTER_NAMES)):
        if i==FULL_PROMOTER_NAMES[j]:
            X.append(np.divide(FULL_eGFP_IO_IN[j], FULL_eGFP_IO_OUT[j]))
            Y.append(np.divide(FULL_eGFP_SPREAD_IN[j], FULL_eGFP_SPREAD_OUT[j]))
            
            X_tDt.append(np.nanmean(np.divide(FULL_tDt_IO_IN[j], FULL_tDt_IO_OUT[j])))
            Y_tDt.append(np.nanmean(np.divide(FULL_tDt_SPREAD_IN[j], FULL_tDt_SPREAD_OUT[j])))
    plt.scatter(np.nanmean(X), np.nanmean(Y))
    MEAN = np.nanmean(X)
    SEM = sp.stats.sem(X)
    plt.plot((MEAN+SEM, MEAN-SEM), (np.nanmean(Y),np.nanmean(Y)), color='black')
    MEAN = np.nanmean(Y)
    SEM = sp.stats.sem(Y)
    plt.plot((np.nanmean(X),np.nanmean(X)), (MEAN+SEM, MEAN-SEM),  color='black')

    plt.text(np.nanmean(X), np.nanmean(Y), i)
plt.scatter(np.nanmean(X_tDt), np.nanmean(Y_tDt))

MEAN = np.nanmean(X_tDt)
SEM = sp.stats.sem(X_tDt)
plt.plot((MEAN+SEM, MEAN-SEM), (np.nanmean(Y_tDt),np.nanmean(Y_tDt)), color='black')
MEAN = np.nanmean(Y_tDt)
SEM = sp.stats.sem(Y_tDt)
plt.plot((np.nanmean(X_tDt),np.nanmean(X_tDt)), (MEAN+SEM, MEAN-SEM),  color='black')
plt.text(np.nanmean(X_tDt), np.nanmean(Y_tDt), i)
plt.xlabel('IO-strength')
plt.ylabel('IO-SPECIFICITY')

plt.figure()
ax = plt.subplot(411)
ax2 = plt.subplot(412)
ax3 = plt.subplot(413)
ax4 = plt.subplot(414)

Y = []
for i in LIST:
    X = []
    for j in range(len(FULL_PROMOTER_NAMES)):
        if i==FULL_PROMOTER_NAMES[j]:
            X.append(FULL_eGFP_SPREAD_IN[j]*pix_to_um_conversion_factor)
    Y.append(X)
ax.boxplot(Y, labels=np.array(LIST))
ax.set_xlabel('Promoter ID')
ax.set_ylabel('IO Surface (um2)')

Y = []
for i in LIST:
    X = []
    for j in range(len(FULL_PROMOTER_NAMES)):
        if i==FULL_PROMOTER_NAMES[j]:
            X.append(FULL_eGFP_IO_IN[j])
    Y.append(X)
ax2.boxplot(Y, labels=np.array(LIST))
ax2.set_xlabel('Promoter ID')
ax2.set_ylabel('Average staining intensity outside IO')

Y = []
for i in LIST:
    X = []
    for j in range(len(FULL_PROMOTER_NAMES)):
        if i==FULL_PROMOTER_NAMES[j]:
            X.append(FULL_eGFP_SPREAD_OUT[j]*pix_to_um_conversion_factor)
    Y.append(X)
ax3.boxplot(Y, labels=np.array(LIST))
ax3.set_xlabel('Promoter ID')
ax3.set_ylabel('IO Surface (um2)')

Y = []
for i in LIST:
    X = []
    for j in range(len(FULL_PROMOTER_NAMES)):
        if i==FULL_PROMOTER_NAMES[j]:
            X.append(FULL_eGFP_IO_OUT[j])
    Y.append(X)
ax4.boxplot(Y, labels=np.array(LIST))
ax4.set_xlabel('Promoter ID')
ax4.set_ylabel('Average staining intensity outside IO')