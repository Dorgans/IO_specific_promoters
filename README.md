# Designing AAV vectors for in vivo monitoring of the subtle calcium fluctuations of inferior olive network

ARTICLE:
<a href='https://www.frontiersin.org/articles/10.3389/fncel.2022.825056/abstract'>https://www.frontiersin.org/articles/10.3389/fncel.2022.825056/abstract</a>

RAW DATA:
<a href='https://www.ebi.ac.uk/biostudies/studies/S-BIAD408'>https://www.ebi.ac.uk/biostudies/studies/S-BIAD408</a>

 <it>Code used for comparing expression levels and signal properties of transgenes expressed via various viral vectors from Kiyoto Kurima.</it>

- #### Mostly confocal image analysis
<img src="https://user-images.githubusercontent.com/46438160/126939760-7795f9b8-161c-474d-a118-ddb932a70ef6.png" alt="IOn IHC (MAP2, Alpha-tub, NeuN)" width="300" height="300">

- #### Or in-vitro / in-vivo calcium imaging analysis
<img src="https://user-images.githubusercontent.com/46438160/126941800-4bfd64ff-befc-4868-ab29-530474372caf.png" alt="AAV.PHP.eB, dual-AAV, GCamp7" height="200">


##### Datasets are available for analysis

<a href='https://github.com/Dorgans/IO_specific_promoters/blob/master/20200304_ANATOMY_HighMagnificationAnalysis.py'>60x Anatomy (Figure2)<a>
 <br/>
<a href='https://github.com/Dorgans/IO_specific_promoters/blob/master/20210328_IHC_NeuroAstroIdentification.py'> IHC Neuron/Astro identification (Figure1)</a>
 <br/>
<a href='https://github.com/Dorgans/IO_specific_promoters/blob/master/20210619_SOMA_SPIKE_ANALYSIS.py'>Spike-related analysis (Figure6,7,8)</a>
 <br/>
 
    rem
    code sequences can be launched individually and require the corresponding datasets uploaded in the repository
    each chunk of code (1) imports (2) calculates (3) plots data
    some code is just tentative code but I kept it. I will clear everything un-used when project is closed.
    DATA folder does not contain raw data but extracted data from raw (such as average soma R.O.I intensity, timeseries average for soma pixels ..)
    will link here the RAW data repositor(ies) if there are

<br/>
<br/>

##### Code / data correspondance

    20210726
    ###### 20200828_IO_spe_promoters_main.py               >> datasets "IO_Spe_Prom_ANATOMY_5x" 
    ###### 20201214_IO_spe_promoters_60x.py                >> datasets "IO_Spe_Prom_ANATOMY_60x_IOPr", "IO_Spe_Prom_ANATOMY_60x_IOM" 
    ###### 20210328_IHC_NeuroAstroIdentification.py        >> datasets "IHC_CellSurface"
    ###### 20210619_SOMA_SPIKE_ANALYSIS.py                 >>  datasets "FIG_InVitro_IOSPIKEActivity_enGCamp6s\CSV_COMPIL_FOR pyScript_[...]" 
    ###### 20210619_SOMA_STO_ANALYSIS.py                   >> datasets "FIG_InVitro_IOSTOActivity_enGCamp6s\CSV_COMPIL_FOR pyScript_[...]" 
    ###### 20210220_STO_PSD_AREA.py                        >> UNUSED BUT WAS INFORMATIVE
    ###### 20200304_ANATOMY_HighMagnificationAnalysis.py   >> UNUSED BUT WAS INFORMATIVE
