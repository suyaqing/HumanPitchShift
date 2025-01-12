%%% This File reads in the Praat txt-files (rawdata) and is taken partially from Mengli Feng.
%%% The code was altered, so that it potentially works with Praat rawdata using a different sample rate.

% Author: Mengli Feng
% Adapted by: Philipp Eugster <eugsteph@student.ethz.ch>
% Final changes: 18.12.2022

% adapted by Yaqing Su, extending time axes to -1s and align to threshold-based onsets
% last changes: 26.02.2024


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% For each Participant the data for the 4 different Tests is read in.
%%% candidates EGG files are excluded.

ons = load([savepath 'OnsetTimes.mat']);
r = FS/TIMESTEPS_PER_SECOND; % convert sample index from raw to pitch
disp("Reading in Data...");
for b=1:nsubj
    fprintf('process: %s ...\n',sbjnames{b});
       
    wavdir = [datapath sbjnames{b} filesep];%readpathwav;
    f0dir = [datapath 'ana_' sbjnames{b} filesep]; % stores path to already extraced f0 pitch traces.
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ IN TEST A %%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get names and directories of files to be read in (and EGG files to be excluded)
    fname = sprintf('%s_A1*.txt',sbjnames{b});
    fnEGG = sprintf('%s_A1*_EGG*.txt',sbjnames{b});
    d = dir([f0dir fname]);
    dEGG = dir([f0dir fnEGG]);
   
    %Check if data is found
    if isempty(d)
        fprintf('No data found for: %s ...\n',sbjnames{b});
    %    
    else
        fn = {d.name}';
        fEGG = {dEGG.name}';
        fn = setdiff(fn,fEGG); % delete EGG filenames from all filenames
        ntrials = size(fn,1);
        testA(b).Sbjname = sbjnames(b);
        testA(b).Ref_freq = ref_freq(b);
        testA(b).raw_HZ = nan(LENGTH,ntrials);
        testA(b).Onsets = ceil(ons.testA(b).onsets/r);
        for f=1:ntrials
            onset = testA(b).Onsets(f);
            start = TPS - onset + 2; % t=0 is at the 101th sample
            dat = readpraatdata_YS([f0dir fn{f}],ref_freq(b),REC_DURATION,TIMESTEPS_PER_SECOND); %Reads in the praat data, removes small voice fragments and corrects for pitch-shift
            % pitch traces are aligned to the start of audio recording,
            % realign now to the start of vocalization
            if start>0
                testA(b).raw_HZ(start:start+length(dat)-1,f) = dat;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ IN TEST B %%%%%%%%%%%%%%%%%%%%%%%%%%

    testB(b).Sbjname = sbjnames(b);
    testB(b).Ref_freq = ref_freq(b);
    
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
        d = dir([f0dir fname]);
        dEGG = dir([f0dir fnEGG]);
        fn = {d.name}';
        fEGG = {dEGG.name}';
        fn = setdiff(fn,fEGG);
        ntrialsEGG = size(fEGG,1);
        testB(b).raw = nan(LENGTH,ntrials);
        testB(b).Onsets = ceil(ons.testB(b).onsets/r);
        for f=1:ntrials
            onset = testB(b).Onsets(f);
            start = TPS - onset + 2; % t=0 is at the 101th sample
            dat = readpraatdata_YS([f0dir fn{f}],ref_freq(b),REC_DURATION,TIMESTEPS_PER_SECOND);
            if start > 0
                testB(b).raw(start:start+length(dat)-1,f) = dat;
            end

        end
        testB(b).target = target;
        
    else
        testB(b).raw = [];
        testB(b).target = [];
        testB(b).Onsets = [];
    end
    
%     testB(31).raw(:,18) = nan;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%% READ IN TEST C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    testC(b).Sbjname = sbjnames(b);
    testC(b).Ref_freq = ref_freq(b);


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
            dA = dir([f0dir fnameRA]);
            dAEGG = dir([f0dir fnameRAEGG]);
            fnA = [fnA; {dA.name}'];
            fnAEGG = [fnAEGG; {dAEGG.name}'];
            fnA = setdiff(fnA,fnAEGG);
            ntrialsA = size(dA,1);
            ntrialsAEGG = size(dAEGG,1);
            dB = dir([f0dir fnameRB]);
            dBEGG = dir([f0dir fnameRBEGG]);
            fnB = [fnB; {dB.name}'];
            fnBEGG = [fnBEGG; {dBEGG.name}'];
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
        testC(b).OnsetsA = [];
        testC(b).OnsetsB = [];
    else
        testC(b).rawA = nan(TPS+TIMESTEPS_PER_SECOND*REC_DURATION,ntrialsA*s);
        testC(b).rawB = nan(TPS+TIMESTEPS_PER_SECOND*REC_DURATION,ntrialsB*s);        
        testC(b).mA = nan(ntrialsA*s);        
        testC(b).mB = nan(ntrialsB*s);
        testC(b).OnsetsA = ceil(ons.testC(b).onsetsA/r);
        testC(b).OnsetsB = ceil(ons.testC(b).onsetsB/r);
        if ntrialsAEGG == ntrialsA
            testC(b).rawAEGG = nan(TIMESTEPS_PER_SECOND*REC_DURATION,ntrialsA*s);
            testC(b).mAEGG = nan(ntrialsA*s);     
        else
            testC(b).rawAEGG  = [];
            testC(b).mAEGG  =[]; 
        end
        if ntrialsBEGG == ntrialsB
            testC(b).rawBEGG = nan(TIMESTEPS_PER_SECOND*REC_DURATION,ntrialsB*s);
            testC(b).mBEGG = nan(ntrialsB*s);
        else
            testC(b).rawBEGG = [];
            testC(b).mBEGG=[];
        end
        for f=1:ntrialsA*s
            onset = testC(b).OnsetsA(f);
            start = TPS - onset + 2; % t=0 is at the 101th sample
            dat = readpraatdata_YS([f0dir fnA{f}],ref_freq(b),REC_DURATION,TIMESTEPS_PER_SECOND);
            if start > 0
                testC(b).rawA(start:start+length(dat)-1,f) = dat;
            end
%             if ntrialsAEGG == ntrialsA
%                 testC(b).rawAEGG(:,f) = readpraatdata_YS([f0dir fnAEGG{f}],ref_freq(b),REC_DURATION,TIMESTEPS_PER_SECOND);
%                 if isempty(find(testC(b).rawAEGG(:,f)>1, 1))
%                     temp = testC(b).rawAEGG(:,f);
%                     temp(~isnan(temp))=nan;
%                     testC(b).rawAEGG(:,f) = temp;
%                 end
%             endr
            % wipe out empty data (empty data(zeros) are changed to ones in readpraatdata)
            if isempty(find(testC(b).rawA(:,f)>1, 1))
                temp = testC(b).rawA(:,f);
                temp(~isnan(temp))=nan;
                testC(b).rawA(:,f) = temp;
            end
        end
        
        for f = 1:ntrialsB*s
            onset = testC(b).OnsetsB(f);
            start = TPS - onset + 2; % t=0 is at the 101th sample
            dat = readpraatdata_YS([f0dir fnB{f}],ref_freq(b),REC_DURATION,TIMESTEPS_PER_SECOND);
            if start > 0
                testC(b).rawB(start:start+length(dat)-1,f) = dat;
            end
%             testC(b).rawB(:,f) = readpraatdata_YS([f0dir fnB{f}],ref_freq(b),REC_DURATION,TIMESTEPS_PER_SECOND);
            if isempty(find(testC(b).rawB(:,f)>1, 1))
                temp = testC(b).rawB(:,f);
                temp(~isnan(temp))=nan;
                testC(b).rawB(:,f) = temp;
            end
%             if ntrialsBEGG == ntrialsB
%                 testC(b).rawBEGG(:,f) = readpraatdata_YS([f0dir fnBEGG{f}],ref_freq(b),REC_DURATION,TIMESTEPS_PER_SECOND);
%                 if isempty(find(testC(b).rawAEGG(:,f)>1, 1))
%                     temp = testC(b).rawAEGG(:,f);
%                     temp(~isnan(temp))=nan;
%                     testC(b).rawAEGG(:,f) = temp;
%                 end
%             end       
        end
        testC(b).targetA = targetA;
        testC(b).targetB = targetB;
        testC(b).shift = shift;
    end  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ IN TEST D %%%%%%%%%%%%%%%%%%%%%%%%%%


%     target = []; shift = []; fn = [];
% 
%     for s=1:SESSIONS_D
%         % get conditions
%         d = dir([wavdir sbjnames{b} '_D' num2str(s) '_list.txt']);
%         fid = fopen([wavdir d.name],'r');
%         if fid == -1
%             break
%         else
%             C = textscan(fid,'%d\t%d\t%d','headerlines',4);
%             fclose(fid);
%             ntrials = double(max(C{1}));
%             target = [target; double(C{2})];
%             shift = [shift; double(C{3})];
%             fname = sprintf('%s_D%d*.mat',sbjnames{b},s);
%             d = dir([wavdir fname]);
%             fn = [fn; {d.name}'];
%         end
%     end
%     if fid == -1
%         testD(b).tuning = [];
%         testD(b).raw = [];
%         testD(b).target = [];
%         testD(b).shift = [];
%     else
%         testD(b).tuning = nan(TIMESTEPS_PER_SECOND*REC_DURATION_D,ntrials*(s));
%         testD(b).raw = nan(TIMESTEPS_PER_SECOND*REC_DURATION_D,ntrials*(s));
%         for f=1:ntrials*(s)
%             [~,~,testD(b).tuning(:,f)] = readknobdata([wavdir fn{f}]);
%             testD(b).raw(:,f) = testD(b).tuning(:,f) + shift(f) + target(f);
%             testD(b).raw(:,f) = 2.^(testD(b).raw(:,f)/1200)*ref_freq(b);
%         end
%         testD(b).target = target;
%         testD(b).shift = shift;
%     end
end
%%
getPath;
save(savepath + "RAWDATA_newOn.mat",'testA','testB','testC');
clear all
getPath;

load(savepath + "RAWDATA_newOn.mat");
load(savepath + "variables.mat");


disp("All data read in & saved at: ");
disp(savepath + "RAWDATA.mat");

