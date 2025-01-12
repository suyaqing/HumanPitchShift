
%%% Main file from which the program is run. Only file that needs input/manipulation by user.
%%% The parameters in the section "PROGRAM FLOW" need user input. The parameters in the "PARAMETER" section are pre-set but can be changed if wanted.



% Author: Philipp Eugster <eugsteph@student.ethz.ch>
% Created: 18.12.2022
% Modified by Yaqing Su


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% START CHANGE HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PROGRAM FLOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set ans save all variables used for reading in files and extract
%parameters


%%% Set "true" if rawdata is read in (may because WINDOW_OUTLIER is changed). If set to "true" set the "datapath" variable to the FOLDER, where the rawdata (PRAAT txt files) are stored.
%%% If set to "false", the data is read in from the PARAMETERS_CENT.m file directly. In that case set the "datapath" variable to the FOLDER, where the PARAMETER_CENT.m file is stored.
READ_RAWDATA = true;
datapath = 'D:\Pitch adaptation\dataFromMengli\'; %Make sure to include the last backslash symbol, so the path points into the folder!

%%% Set to "true" if for each participant the rawdata for the different tests should be ploted. The format and resolution can be set below.
PRODUCE_PLOTS = true;

%%% If single plots per participant and Test should be created, set to "true". If all tests per participant should be plotted onto one page, set to "false"; The format and resolution can be set below.
%%% Only works if PRODUCE_PLOTS is set to "true".
SINGLE_PLOTS = false;

%%% Store the path, where all X.m files and potential Plots should be saved. If plots are created they will be stored in a Subfolder of "savepath" called "\Plots\All" which will be created in the script "Plot_Rawdata_All.m"
savepath = 'D:\Pitch adaptation\test\'; %Make sure to include the last backslash symbol, so the path points into the folder!



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%% Set Resolution and Format for Plots
RESOLUTION  = 200; %Resolution in dots per inch (DPI), specified as a whole number that is greater than or equal to 1.
FORMAT = ".pdf"; %Possible formats are: '.jpg', '.png', '.tif', '.gif', '.pdf', '.emf' or '.eps'
CONTENT_TYPE = 'auto'; % Type of content to store when saving as an EMF, EPS, or PDF file. Specify the value as one of these options:
                         % 'auto' — MATLAB controls whether the content is a vector graphic or an image.
                         % 'vector' — Stores the content as a vector graphic that can scale to any size. If you are saving a PDF file, embeddable fonts are included in the file.
                         % 'image' — Rasterizes the content into one or more images within the file.



%Set all parameters
%%% Set Window for outlier detection. Window 1 to 5 are used for tests A, B & C and WINDOW_OUTLIER_D for test D.
WINDOW_OUTLIER_1 = (40:80)+100;
WINDOW_OUTLIER_2 = (60:100)+100;
WINDOW_OUTLIER_3 = (80:120)+100;
WINDOW_OUTLIER_4 = (100:140)+100;
WINDOW_OUTLIER_5 = (120:150)+100;
WINDOW_OUTLIER_D = (300:370)+100;

%,'f02' SUBJECT f02 removed, see readme in TestA_Phil.m
%'f20NC' SUBJECT f20NC removed, see readme in TestA_Phil.m
%If both will be included, copy sbjnames from "all_analysis" file from
%Mengli.
sbjnames = {'m02NC','f03NC','f04','m03','f05','m04',...
    'm05', 'f06', 'f07', 'f08','f09', 'm06', 'f10','f11NC','m07','f12NC',...
    'm08','m09NC','f13NC','m10NC','m11NC','m12NC','m13NC','f14NC','m14NC',...
    'm15', 'm16NC','m17NC','m18NC','m19','m20NC','f15','m21NC', 'f16','f17',...
    'f18NC','m22','m23','f19NC', 'f20','m24NC', 'm25NC','m26','f21NC',...
    'f22NC','f23','f24','f25NC','f26NC','m27','f27NC'};

% sbjnames = {'m02NC','f03NC','f04','m03','f05','m04',...
%     'm05', 'f06', 'f07', 'f08'};
% sbjnames = {'f08'};
% 
% sbjnames = {'m02NC','f03NC'};

nsubj = length(sbjnames); %Stores the number of subjects.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Only change this parameter if the sample rate for the Praat extraction of the pitch values changed.

TIMESTEPS_PER_SECOND = 100; %Datapoints in extractet pitch file (.txt) per second (used in all tests).
TPS = TIMESTEPS_PER_SECOND;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STOPT CHANGE HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Creates a script that saves the provided datapath and savepath
fid = fopen( 'getPath.m', 'wt' );
fprintf(fid, 'savepath = ''%s'';\n', savepath);
fprintf(fid, 'datapath = ''%s'';\n', datapath);
fclose(fid);

%%% Sets variables that are fixed and should not be changed.


%%%Add all the subfolders to the matlab-path
% addpath("Initialize\" ,genpath("HypothesisTests\"), genpath("PlotRawdata\"), "ProcessRawdata\", "ReadRawdata\") 
%%
%%%& Calls scripts accordingly.
if READ_RAWDATA == true
    Set_Variables_YS; 
    ReadRawData;
    Prepare_Data_Cent;
else
    load(savepath + "PARAMETERS_CENT.mat");
    Set_Variables_YS
end

if PRODUCE_PLOTS == true
    if SINGLE_PLOTS == true
        Plot_Rawdata_Cent_A;
        Plot_Rawdata_Cent_B;
        Plot_Rawdata_Cent_C;
        Plot_Rawdata_Cent_D;
    end
    if SINGLE_PLOTS == false
        Plot_Rawdata_Cent_PitchTraces;
    end
end


%%% Load already read in rawdata:
% filename = [savepath 'variables' '.mat'];
% load(filename);
% 
% filename = [savepath 'RAWDATA' '.mat'];
% load(filename);





