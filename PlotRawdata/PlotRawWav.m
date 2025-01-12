% plot raw audio traces separated by target
% author: Yaqing Su
% last update: 14/02/2024
sbjnames = {'m02NC','f03NC','f04','m03','f05','m04',...
    'm05', 'f06', 'f07', 'f08','f09', 'm06', 'f10','f11NC','m07','f12NC',...
    'm08','m09NC','f13NC','m10NC','m11NC','m12NC','m13NC','f14NC','m14NC',...
    'm15', 'm16NC','m17NC','m18NC','m19','m20NC','f15','m21NC', 'f16','f17',...
    'f18NC','m22','m23','f19NC', 'f20','m24NC', 'm25NC','m26','f21NC',...
    'f22NC','f23','f24','f25NC','f26NC','m27','f27NC'};
% sbjnames = {'f04'};
nsubj = length(sbjnames); 
datapath = 'D:\Pitch adaptation\dataFromMengli\Yaqing\';
savepath = 'D:\Pitch adaptation\testYS\';
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

% fs = 44100;
col_B = [[0 0.4470 0.7410];[0.4940 0.1840 0.5560];]; 
ws = 441;

% testA = [];
% testB = [];
% testC = [];
%%
fs_win = FS/ws;
for b=1:nsubj

    wavdir = [datapath sbjnames{b} filesep];
%     f0dir = [datapath 'ana_' sbjnames{b} filesep]; % stores path to already extraced f0 pitch traces.
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ IN TEST A %%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get names and directories of files to be read in (and EGG files to be excluded)

    wavname = sprintf('%s_A1*.wav',sbjnames{b});
    EGGwavname = sprintf('%s_A1*_EGG.wav',sbjnames{b});

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
        testA(b).ma = [];

%         testA(b).raw = nan(FS,ntrials);
%         testA(b).rawenv = nan(FS,ntrials);
        for f=1:ntrials
            [y, ~] = audioread([wavdir fn{f}]);
            testA(b).ma = [testA(b).ma moving_average(y, ws)];
        end
        tt = (0:length(testA(b).ma)-1)/fs_win;

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% READ IN TEST B %%%%%%%%%%%%%%%%%%%%%%%%%%

    testB(b).Sbjname = sbjnames(b);
    testB(b).ma0 = [];
    testB(b).ma400 = [];
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
        
        for f=1:ntrials
            [y, ~] = audioread([wavdir fn{f}]);
            if target(f)==0
                testB(b).ma0 = [testB(b).ma0 moving_average(y, ws)];
            else
                testB(b).ma400 = [testB(b).ma400 moving_average(y, ws)];
            end
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
    testC(b).maA0 = []; testC(b).maA400 = [];
    testC(b).maB0 = []; testC(b).maB400 = [];

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


            wavnameRA = sprintf('%s_C%d*A_real.wav',sbjnames{b},s);
            wavnameRB = sprintf('%s_C%d*B_real.wav',sbjnames{b},s);
            EGGwavnameRA = sprintf('%s_C%d*A_EGG.wav',sbjnames{b},s);
            EGGwavnameRB = sprintf('%s_C%d*B_EGG.wav',sbjnames{b},s);

            dAwav = dir([wavdir wavnameRA]);
            dAEGGwav = dir([wavdir EGGwavnameRA]);
            fnA = [fnA; {dAwav.name}'];
            fnAEGG = [fnAEGG; {dAEGGwav.name}'];
            fnA = setdiff(fnA,fnAEGG);
            ntrialsA = size(dAwav,1);
            ntrialsAEGG = size(dAEGGwav,1);


            dBwav = dir([wavdir wavnameRB]);
            dBEGGwav = dir([wavdir EGGwavnameRB]);
            fnB = [fnB; {dBwav.name}'];
            fnBEGG = [fnBEGG; {dBEGGwav.name}'];
            fnB = setdiff(fnB,fnBEGG);
            ntrialsB = size(dBwav,1);
            ntrialsBEGG = size(dBEGGwav,1);

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
        
        for f=1:ntrialsA*s
            [y, ~] = audioread([wavdir fnA{f}]);
            if targetA(f)==0
                testC(b).maA0 = [testC(b).maA0 moving_average(y, ws)];
            else
                testC(b).maA400 = [testC(b).maA400 moving_average(y, ws)];
            end
            
        end
        
        for f = 1:ntrialsB*s
            [y, ~] = audioread([wavdir fnB{f}]);
            if targetA(f)==0
                testC(b).maB0 = [testC(b).maB0 moving_average(y, ws)];
            else
                testC(b).maB400 = [testC(b).maB400 moving_average(y, ws)];
            end
%                         
        end
        testC(b).targetA = targetA;
        testC(b).targetB = targetB;
        
        testC(b).shift = shift;
    end  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    range = 0.5*fs_win:1.5*fs_win;
    fig=figure;
    sgtitle(sbjnames{b})
    subplot(241)
    plot(tt, testA(b).ma, Color=col_B(1, :));
    title("Test A")
    xlabel("Time (s)")

    
    subplot(242)
    [~, pB] = ttest2(testB(b).ma0', testB(b).ma400');
    hold on
    pl1 = plot(tt, testB(b).ma0, Color=col_B(1, :));
    pl2 = plot(tt, testB(b).ma400, Color=col_B(2, :));
    hold off
    legend([pl1(1), pl2(1)], '0 Cent', '400 Cent')
    title(["Test B"])
    xlabel("Time (s)")
    subplot(246)
    plot(tt, -log10(pB));
    title("t-test")
    xlabel("Time (s)")
    ylabel("-log10(p)")

    
    subplot(243)
    [~, pCA] = ttest2(testC(b).maA0', testC(b).maA400');
    hold on
    pl1 = plot(tt, testC(b).maA0, Color=col_B(1, :));
    pl2 = plot(tt, testC(b).maA400, Color=col_B(2, :));
    title(["Test CA"])
    xlabel("Time (s)")
    legend([pl1(1), pl2(1)], '0 Cent', '400 Cent')
    subplot(247)
    plot(tt, -log10(pCA));
    title("t-test")
    xlabel("Time (s)")
    ylabel("-log10(p)")
    
    subplot(244)
    [~, pCB] = ttest2(testC(b).maB0', testC(b).maB400');
    hold on
    pl1 = plot(tt, testC(b).maB0, Color=col_B(1, :));
    pl2 = plot(tt, testC(b).maB400, Color=col_B(2, :));
    title(["Test CB"])
    xlabel("Time (s)")
    legend([pl1(1), pl2(1)], "0 Cent", "400 Cent")
    subplot(248)
    plot(tt, -log10(pCB));
    title("t-test")
    xlabel("Time (s)")
    ylabel("-log10(p)")
%     
    fig.WindowState = 'maximized';
    exportgraphics(fig, folder + sbjnames{b} + "_RawAudio_Target" + FORMAT, "Resolution",RESOLUTION,'ContentType',CONTENT_TYPE)
    close(fig);
%  
end

function sig_out = moving_average(sig_in, ws)
len = length(sig_in);
% sig_out = zeros(size(sig_in));
% frameAve = sig;

nw = floor(len/ws);
sig_out = zeros(nw, size(sig_in, 2));

for w = 1:nw
    start = (w-1)*ws;
%     if start < len-nw*ws
        idx = (start+1):(start+ws);
        sig_out(w, :) = mean(abs(sig_in(idx, :)));
%     end
end

end