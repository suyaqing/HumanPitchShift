%%% This File reads raw .wav audio recordings find vocal onset

% Author: Yaqing Su
% Log of changes: 17.01.2024 using fixed threshold


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
savepath = 'D:\Pitch adaptation\test\';
FS = 44100;
TARGET = [0 400]; % The two targets used in the experiment per participant in cent.
SHIFT = [-100 -50 -25 0 25 50 100]; % The different shifts in cent that are used in Test CB and also in D.
% thresCB = 1e-3;
SESSIONS_C = 5; % Number of sessions in task C in Experiment perfomred by Mengli

RESOLUTION  = 100; %Resolution in dots per inch (DPI), specified as a whole number that is greater than or equal to 1.
FORMAT = ".pdf"; %Possible formats are: '.jpg', '.png', '.tif', '.gif', '.pdf', '.emf' or '.eps'
CONTENT_TYPE = 'auto';
% savepath = 'D:\Pitch adaptation\test\'; 
folder = savepath + "\Plots\ALL\";
load([savepath 'VocalThresCB.mat'])
fs = 44100;
%%
disp("Reading in Data...");
for b=1:nsubj
    thresCB = 0.001;
%     thresCB = thres(b).thresCB;
    fprintf('process: %s ...\n',sbjnames{b});
       
    wavdir = [datapath sbjnames{b} filesep];%readpathwav;
%     wavdir2 = [datapath 'Yaqing\' sbjnames{b} filesep];
%     f0dir = [datapath 'ana_' sbjnames{b} filesep]; % stores path to already extraced f0 pitch traces.
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ IN TEST A %%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get names and directories of files to be read in (and EGG files to be excluded)
%     fname = sprintf('%s_A1*.txt',sbjnames{b});
%     fnEGG = sprintf('%s_A1*_EGG*.txt',sbjnames{b});
    wavname = sprintf('%s_A1*.wav',sbjnames{b});
    EGGwavname = sprintf('%s_A1*_EGG.wav',sbjnames{b});
%     d = dir([f0dir fname]);
%     dEGG = dir([f0dir fnEGG]);
    dwav = dir([wavdir wavname]);
    dEGGwav = dir([wavdir EGGwavname]);
   
    %Check if data is found
    if isempty(dwav)
        fprintf('No data found for: %s ...\n',sbjnames{b});
    %    
    else
        fn = {dwav.name}';
        fEGG = {dEGGwav.name}';
        fn = setdiff(fn,fEGG); % delete EGG filenames from all filenames
        ntrials = size(fn,1);
        onsets = -ones(1, ntrials);
        testA(b).Sbjname = sbjnames(b);

%         testA(b).raw = nan(FS,ntrials);
%         testA(b).rawenv = nan(FS,ntrials);
        for f=1:ntrials
            [y, ~] = audioread([wavdir fn{f}]);
            [sig, onsets(f)] = apply_threshold(thresCB, y, 1);
%             audiowrite([wavdir2 fn{f}], sig, FS);
%             testA(b).rawenv(:,f) = envelope(y(1:FS), 4410, 'rms');
        end

        testA(b).onsets = onsets;
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
%         fname = sprintf('%s_B1*.txt',sbjnames{b});
%         fnEGG = sprintf('%s_B1*_EGG*.txt',sbjnames{b});
        wavname = sprintf('%s_B1*.wav',sbjnames{b});
        EGGwavname = sprintf('%s_B1*_EGG*.wav',sbjnames{b});
%         d = dir([f0dir fname]);
%         dEGG = dir([f0dir fnEGG]);
        dwav = dir([wavdir wavname]);
        dEGGwav = dir([wavdir EGGwavname]);
        fn = {dwav.name}';
        fEGG = {dEGGwav.name}';
        fn = setdiff(fn,fEGG);
        onsets = -ones(1, ntrials);
%         ntrialsEGG = size(fEGG,1);
%         testB(b).raw = nan(FS,ntrials);
%         testB(b).rawenv = nan(FS,ntrials);
        for f=1:ntrials
            [y, ~] = audioread([wavdir fn{f}]);
            [sig, onsets(f)] = apply_threshold(thresCB, y, 1);
%             audiowrite([wavdir2 fn{f}], sig, FS);
%             testB(b).raw(:,f) = y(1:FS);
%             testB(b).rawenv(:,f) = envelope(y(1:FS), 4410, 'rms');
        end
        testB(b).target = target; 
        testB(b).onsets = onsets;
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

            
%             fnameRA = sprintf('%s_C%d*A_real.txt',sbjnames{b},s);
%             fnameRB = sprintf('%s_C%d*B_real.txt',sbjnames{b},s);
%             fnameRAEGG = sprintf('%s_C%d*A_EGG.txt',sbjnames{b},s);
%             fnameRBEGG = sprintf('%s_C%d*B_EGG.txt',sbjnames{b},s);
            wavnameRA = sprintf('%s_C%d*A_real.wav',sbjnames{b},s);
            wavnameRB = sprintf('%s_C%d*B_real.wav',sbjnames{b},s);
            EGGwavnameRA = sprintf('%s_C%d*A_EGG.wav',sbjnames{b},s);
            EGGwavnameRB = sprintf('%s_C%d*B_EGG.wav',sbjnames{b},s);
%             dA = dir([f0dir fnameRA]);
%             dAEGG = dir([f0dir fnameRAEGG]);
            dAwav = dir([wavdir wavnameRA]);
            dAEGGwav = dir([wavdir EGGwavnameRA]);
            fnA = [fnA; {dAwav.name}'];
            fnAEGG = [fnAEGG; {dAEGGwav.name}'];
            fnA = setdiff(fnA,fnAEGG);
            ntrialsA = size(dAwav,1);
            ntrialsAEGG = size(dAEGGwav,1);

%             dB = dir([f0dir fnameRB]);
%             dBEGG = dir([f0dir fnameRBEGG]);
            dBwav = dir([wavdir wavnameRB]);
            dBEGGwav = dir([wavdir EGGwavnameRB]);
            fnB = [fnB; {dBwav.name}'];
            fnBEGG = [fnBEGG; {dBEGGwav.name}'];
            fnB = setdiff(fnB,fnBEGG);
            ntrialsB = size(dBwav,1);
            ntrialsBEGG = size(dBEGGwav,1);

            onsetsA = -ones(1, ntrialsA);
            onsetsB = -ones(1, ntrialsB);
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
%         testC(b).rawA = nan(FS,ntrialsA*s);
%         testC(b).rawB = nan(FS,ntrialsB*s); 
%         testC(b).rawenvA = nan(FS,ntrialsA*s);
%         testC(b).rawenvB = nan(FS,ntrialsB*s);
     
        
        for f=1:ntrialsA*s
            [y, ~] = audioread([wavdir fnA{f}]);
            [sig, onsetsA(f)] = apply_threshold(thresCB, y, 1);
%             audiowrite([wavdir2 fnA{f}], sig, FS);
%             testC(b).rawA(:,f) = y(1:FS);
%             testC(b).rawenvA(:,f) = envelope(y(1:FS), 4410, 'rms');
            
        end
        
        for f = 1:ntrialsB*s
            [y, ~] = audioread([wavdir fnB{f}]);
%             [sig, onsetsB(f)] = apply_threshold(thresCB, y, 0);
%             audiowrite([wavdir2 fnB{f}], y, FS);
            onsetsB(f) = 1;
%             testC(b).rawB(:,f) = y(1:FS);
%             testC(b).rawenvB(:,f) = envelope(y(1:FS), 4410, 'rms');
            
        end
        testC(b).targetA = targetA;
        testC(b).targetB = targetB;
        testC(b).onsetsA = onsetsA;
        testC(b).onsetsB = onsetsB;
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

save([savepath 'OnsetTimes.mat'], 'testA', 'testB', 'testC');


function [sig, onset] = apply_threshold(thres, rawsig, shift)
% rawsig has the dimension of sample*trial
rawsig = rawsig(:);
len = length(rawsig);
sig = zeros(size(rawsig));
% frameAve = sig;
winAve = sig;

fl= 64;
wl = 100;

nf = floor(len/fl);
for fr = 1:nf-(wl-1)
    start = (fr-1)*fl;
    if start < len-fl*wl
        idx = (start+1):(start+fl*wl);
        winAve(idx) = mean(abs(rawsig(idx)));
    end
end
winAve(idx(end)+1:len)= mean(abs(rawsig(idx)));


for ntrial = 1:size(rawsig, 2)
    indOn = find(winAve>=thres, 1);
    if ~isempty(indOn)
        onset = indOn + shift*fl*(wl-1); % onset is the index here
        nvalid = len-onset+1;
        sig(1:nvalid) = rawsig(onset:len, ntrial);
    else % the whole signal is zero
        onset = nan;
    end
    
end

end

function winAve = compute_winAve(rawsig, fl, wl)
len = length(rawsig);
sig = zeros(size(rawsig));
% frameAve = sig;
winAve = sig;

nf = floor(len/fl);
% for fr = 1:nf
%     idx = ((fr-1)*fl+1):fr*fl;
%     frameAve(idx, :) = mean(abs(rawsig(idx, :)));
% end
% lastFr = (nf*fl+1):len;
% frameAve(lastFr, :) = mean(abs(rawsig(lastFr, :)));
for fr = 1:nf-(wl-1)
    start = (fr-1)*fl;
    if start < len-fl*wl
        idx = (start+1):(start+fl*wl);
        winAve(idx, :) = ones(fl*wl, 1)*mean(abs(rawsig(idx, :)));
    end
end
winAve(idx(end)+1:len, :)= ones(length(idx(end)+1:len), 1)*mean(abs(rawsig(idx, :)));


end