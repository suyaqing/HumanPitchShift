%%% This File reads raw .wav audio recordings and plot for visual
%%% inspection

% Author: Yaqing Su
% Final changes: 10.01.2024


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% For each Participant the data for the 4 different Tests is read in.
%%% candidates EGG files are excluded.
sbjnames = {'m02NC','f03NC','f04','m03','f05','m04',...
    'm05', 'f06', 'f07', 'f08','f09', 'm06', 'f10','f11NC','m07','f12NC',...
    'm08','m09NC','f13NC','m10NC','m11NC','m12NC','m13NC','f14NC','m14NC',...
    'm15', 'm16NC','m17NC','m18NC','m19','m20NC','f15','m21NC', 'f16','f17',...
    'f18NC','m22','m23','f19NC', 'f20','m24NC', 'm25NC','m26','f21NC',...
    'f22NC','f23','f24','f25NC','f26NC','m27','f27NC'};
% sbjnames = {'f04'};
nsubj = length(sbjnames); 
datapath = 'D:\Pitch adaptation\dataFromMengli\';
FS = 44100;
TARGET = [0 400]; % The two targets used in the experiment per participant in cent.
SHIFT = [-100 -50 -25 0 25 50 100]; % The different shifts in cent that are used in Test CB and also in D.

SESSIONS_C = 5; % Number of sessions in task C in Experiment perfomred by Mengli

RESOLUTION  = 300; %Resolution in dots per inch (DPI), specified as a whole number that is greater than or equal to 1.
FORMAT = ".pdf"; %Possible formats are: '.jpg', '.png', '.tif', '.gif', '.pdf', '.emf' or '.eps'
CONTENT_TYPE = 'auto';
savepath = 'D:\Pitch adaptation\test\'; 
folder = savepath + "\Plots\ALL\";

disp("Reading in Data...");
for b=1:nsubj
    fprintf('process: %s ...\n',sbjnames{b});
       
    wavdir = [datapath sbjnames{b} filesep];%readpathwav;
    f0dir = [datapath 'ana_' sbjnames{b} filesep]; % stores path to already extraced f0 pitch traces.
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ IN TEST A %%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get names and directories of files to be read in (and EGG files to be excluded)
    fname = sprintf('%s_A1*.txt',sbjnames{b});
    fnEGG = sprintf('%s_A1*_EGG*.txt',sbjnames{b});
    wavname = sprintf('%s_A1*.wav',sbjnames{b});
    EGGwavname = sprintf('%s_A1*_EGG.wav',sbjnames{b});
    d = dir([f0dir fname]);
    dEGG = dir([f0dir fnEGG]);
    dwav = dir([wavdir wavname]);
    dEGGwav = dir([wavdir EGGwavname]);
   
    %Check if data is found
    if isempty(d)
        fprintf('No data found for: %s ...\n',sbjnames{b});
    %    
    else
        fn = {dwav.name}';
        fEGG = {dEGGwav.name}';
        fn = setdiff(fn,fEGG); % delete EGG filenames from all filenames
        ntrials = size(fn,1);
        testA(b).Sbjname = sbjnames(b);
        testA(b).raw = nan(FS,ntrials);
%         testA(b).rawenv = nan(FS,ntrials);
        for f=1:ntrials
            [y, ~] = audioread([wavdir fn{f}]);
            testA(b).raw(:,f) = y(1:FS);
%             testA(b).rawenv(:,f) = envelope(y(1:FS), 4410, 'rms');
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ IN TEST B %%%%%%%%%%%%%%%%%%%%%%%%%%

    testB(b).Sbjname = sbjnames(b);
%     testB(b).Ref_freq = ref_freq(b);
    
    % get conditions
    d = dir([wavdir sbjnames{b} '_B1_list.txt']);
    fid = fopen([wavdir d.name],'r');
    if fid ~= -1
        C = textscan(fid,'%d\t%d','headerlines',4);
        fclose(fid);
        ntrials = double(max(C{1}));
        target = double(C{2});
        % get data
        fname = sprintf('%s_B1*.txt',sbjnames{b});
        fnEGG = sprintf('%s_B1*_EGG*.txt',sbjnames{b});
        wavname = sprintf('%s_B1*.wav',sbjnames{b});
        EGGwavname = sprintf('%s_B1*_EGG*.wav',sbjnames{b});
        d = dir([f0dir fname]);
        dEGG = dir([f0dir fnEGG]);
        dwav = dir([wavdir wavname]);
        dEGGwav = dir([wavdir EGGwavname]);
        fn = {dwav.name}';
        fEGG = {dEGGwav.name}';
        fn = setdiff(fn,fEGG);
        ntrialsEGG = size(fEGG,1);
        testB(b).raw = nan(FS,ntrials);
%         testB(b).rawenv = nan(FS,ntrials);
        for f=1:ntrials
            [y, ~] = audioread([wavdir fn{f}]);
            testB(b).raw(:,f) = y(1:FS);
%             testB(b).rawenv(:,f) = envelope(y(1:FS), 4410, 'rms');
        end
        testB(b).target = target;  
    else
        testB(b).raw = [];
        testB(b).target = [];
    end
    
%     testB(31).raw(:,18) = nan;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%% READ IN TEST C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    testC(b).Sbjname = sbjnames(b);
%     testC(b).Ref_freq = ref_freq(b);


    targetA = []; targetB = []; shift = [];
    fnA = []; fnB = [];
    fnAEGG = []; fnBEGG = [];
    
    for s=1:SESSIONS_C
        % get conditions
        d = dir([wavdir sbjnames{b} '_C' num2str(s) '_list.txt']);
        fid = fopen([wavdir d.name],'r');
        if fid == -1
            break
        end
        if fid ~= -1
            C = textscan(fid,'%d\t%d\t%d\t%d','headerlines',4);
            fclose(fid);
            ntrialsA = double(max(C{1}));
            shift = [shift; double(C{4})];
            targetA = [targetA; double(C{2})];
            targetB = [targetB; double(C{3})];

            
            fnameRA = sprintf('%s_C%d*A_real.txt',sbjnames{b},s);
            fnameRB = sprintf('%s_C%d*B_real.txt',sbjnames{b},s);
            fnameRAEGG = sprintf('%s_C%d*A_EGG.txt',sbjnames{b},s);
            fnameRBEGG = sprintf('%s_C%d*B_EGG.txt',sbjnames{b},s);
            wavnameRA = sprintf('%s_C%d*A_real.wav',sbjnames{b},s);
            wavnameRB = sprintf('%s_C%d*B_real.wav',sbjnames{b},s);
            EGGwavnameRA = sprintf('%s_C%d*A_EGG.wav',sbjnames{b},s);
            EGGwavnameRB = sprintf('%s_C%d*B_EGG.wav',sbjnames{b},s);
            dA = dir([f0dir fnameRA]);
            dAEGG = dir([f0dir fnameRAEGG]);
            dAwav = dir([wavdir wavnameRA]);
            dAEGGwav = dir([wavdir EGGwavnameRA]);
            fnA = [fnA; {dAwav.name}'];
            fnAEGG = [fnAEGG; {dAEGGwav.name}'];
            fnA = setdiff(fnA,fnAEGG);
            ntrialsA = size(dA,1);
            ntrialsAEGG = size(dAEGG,1);

            dB = dir([f0dir fnameRB]);
            dBEGG = dir([f0dir fnameRBEGG]);
            dBwav = dir([wavdir wavnameRB]);
            dBEGGwav = dir([wavdir EGGwavnameRB]);
            fnB = [fnB; {dBwav.name}'];
            fnBEGG = [fnBEGG; {dBEGGwav.name}'];
            fnB = setdiff(fnB,fnBEGG);
            ntrialsB = size(dB,1);
            ntrialsBEGG = size(dBEGG,1);
        else
            ntrialsA = 0;
            ntrialsAEGG = 0;
            ntrialsB = 0;
            ntrialsBEGG = 0;
        end
    end
    if fid == -1
        testC(b).rawA = [];
        testC(b).rawAEGG = [];
        testC(b).rawB = [];
        testC(b).rawBEGG = [];
        testC(b).mB = []; 
        testC(b).mBEGG = [];
        testC(b).mA = [];   
        testC(b).mAEGG = [];   
        testC(b).targetA = []; 
        testC(b).targetB = [];  
        testC(b).shift = [];        
    else
        testC(b).rawA = nan(FS,ntrialsA*s);
        testC(b).rawB = nan(FS,ntrialsB*s); 
%         testC(b).rawenvA = nan(FS,ntrialsA*s);
%         testC(b).rawenvB = nan(FS,ntrialsB*s);
     
        
        for f=1:ntrialsA*s
            [y, ~] = audioread([wavdir fnA{f}]);
            testC(b).rawA(:,f) = y(1:FS);
%             testC(b).rawenvA(:,f) = envelope(y(1:FS), 4410, 'rms');
            
        end
        
        for f = 1:ntrialsB*s
            [y, ~] = audioread([wavdir fnB{f}]);
            testC(b).rawB(:,f) = y(1:FS);
%             testC(b).rawenvB(:,f) = envelope(y(1:FS), 4410, 'rms');
            
        end
        testC(b).targetA = targetA;
        testC(b).targetB = targetB;
        testC(b).shift = shift;
    end  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fig=figure;
%     sgtitle(sbjnames{b})
%     subplot(221)
%     plot(testA(b).rawenv);
%     title("Test A")
%     xlabel("Sample")
%     
%     subplot(222)
%     plot(testB(b).rawenv);
%     title("Test B")
%     xlabel("Sample")
%     
%     subplot(223)
%     plot(testC(b).rawenvA);
%     title("Test CA")
%     xlabel("Sample")
%     
%     subplot(224)
%     plot(testC(b).rawenvB);
%     title("Test CB")
%     xlabel("Sample")
%     
%     fig.WindowState = 'maximized';
%     exportgraphics(fig, folder + sbjnames{b} + "_RawAudioTraces" + FORMAT, "Resolution",RESOLUTION,'ContentType',CONTENT_TYPE)
%     close(fig);
%  
end
%%
% 
% getPath;
% save(savepath + "RAWDATA.mat",'testA','testB','testC','testD');
% clear all
% getPath;
% 
% load(savepath + "RAWDATA.mat");
% load(savepath + "variables.mat");
% 
% 
% disp("All data read in & saved at: ");
% disp(savepath + "RAWDATA.mat");

