The repository consists of MATLAB scripts for two main purposes: 
1) Basic processing of auditory recordings adapted mainly from Philipp Eugster's original code, including
ProcessRawdata/FindWavOnset.m: find onset time from .wav files for task A, B, CA based on threshold and store onset time for each trace
The following three scripts are in different folders and best ran using the Main.m so that parameters are configured correctly:
	Initialize/Set_Variables_YS.m: set global variables
	ReadRawdata/ReadRawData.m: read Praat files, extrace pitch traces in Hz and realign pitch traces to vocalization onset time. Since some traces have negative onset time, the time axes now starts from -1s to 2s assuming t=0 is the vocalization onset.
	ProcessRawdata/Prepare_Data_Cent.m: convert pitch traces to Cent, mark outlier traces, and compute mean and variances for each trace
Plot_Rawdata_Cent_PitchTraces.m: plot pitch traces of all tasks for individual subjects

2) Model fitting and statistic analysis of pitch traces by Yaqing Su, including:
FitLME_all.m: Linear mixed effect model with subject ID as random intercept for all tasks, store models
FitLME_noCB.m: LME without task CB
CompareLMEtoData.m: evaluate LME model and plot reconstructed traces along data (needs updating)
FitIndividualShiftSigmoidInterval.m: fit a sigmoid model for individual subject's traces in task CB
ComparePerformanceBvsC.m: use model parameters to evaluate each subject's task C performace against task B performance
