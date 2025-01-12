%%% This File processes the rawdata and creates a "masterfile" called PARAMETERS_HZ.mat.
%%% Next, Outliers are detected and saved (are excluded from all further calculations).
%%% Parameters are extraced very similar to the PARAMETER_CENT.mat file, explained here: <Calculations_PARAMETERS_CENT.pdf>
%%% The handling of the PARAMETERS_HZ.mat file is similar to the PARAMETERS_CENT.mat file explained here <ReadMe.txt>

% Author: Philipp Eugster <eugsteph@student.ethz.ch>
% Final changes: 21.12.2022





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEST_A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%APPEND PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Append Subject names, the Reference Frequency for each participant and the rawdata for testA in to the PARAMETERS_HZ file.

for s = 1:nsubj
    PARAMETERS_HZ(s).Sbjnames = testA(s).Sbjname;
    PARAMETERS_HZ(s).ref_freq = ref_freq(s);
    PARAMETERS_HZ(s).A_HZ = testA(s).raw_HZ;
    PARAMETERS_HZ(s).TARGET(1) = ref_freq(s); 
    PARAMETERS_HZ(s).TARGET(2) = 2^(1/3) * ref_freq(s); %Equal to the frequency in HZ 4 semitones higher than the ref_freq
    for i = 1:length(SHIFT) % SHIFT =  [-100   -50   -25     0    25    50   100]
        PARAMETERS_HZ(s).TARGET(i) = 2^(SHIFT(i)/1200) * ref_freq(s);
    end

end


%%% Make sure, that at least 10 trials are present at each timepoint. Otherwise the mean and var have no power.

for s = 1:nsubj
    for i = 1:LENGTH
        if sum(~isnan(PARAMETERS_HZ(s).A_HZ(i,:))) < 10
            PARAMETERS_HZ(s).A_HZ(i:end,:) = nan;
            continue
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% REMOVE UNNATURAL JUMPS %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Unnatural jumps in pitch are removed (Can happen at the end of a trial).

%ManualRemovalOutliers_TestA; 



%%%%%%%%%%%%%%%%%%%%%%%%%%DETECT OUTLIERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% DETECT OUTLIERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Detect Outliers by using the isoutlier() matlab function for five different time Windows.
%%% The time windows are specified in the Main.m File.
%%% A trial is flagged as an outlier if it is detected in two or more of the time windows.

for s = 1:nsubj
    outlr1 = isoutlier(mean(PARAMETERS_HZ(s).A_HZ(WINDOW_OUTLIER_1,:),1,'omitnan'),'median');
    outlr2 = isoutlier(mean(PARAMETERS_HZ(s).A_HZ(WINDOW_OUTLIER_2,:),1,'omitnan'),'median');
    outlr3 = isoutlier(mean(PARAMETERS_HZ(s).A_HZ(WINDOW_OUTLIER_3,:),1,'omitnan'),'median');
    outlr4 = isoutlier(mean(PARAMETERS_HZ(s).A_HZ(WINDOW_OUTLIER_4,:),1,'omitnan'),'median');
    outlr5 = isoutlier(mean(PARAMETERS_HZ(s).A_HZ(WINDOW_OUTLIER_5,:),1,'omitnan'),'median');
    outlr_sum = outlr1 + outlr2 + outlr3 + outlr4 + outlr5;
    outlr = outlr_sum >= 2;
    PARAMETERS_HZ(s).A_outlr = outlr;
    PARAMETERS_HZ(s).A_nOutlr = sum(PARAMETERS_HZ(s).A_outlr);
end


%%%%%%%%%%%%%%%%%% CALCULATE MEAN AND VARIANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calcualte MEAN and VARIANCE for testA and append to parameters File
%%% The detected Outliers are exluded for the calculation!.

for s = 1:nsubj
    PARAMETERS_HZ(s).A_mean = mean(PARAMETERS_HZ(s).A_HZ(:,~PARAMETERS_HZ(s).A_outlr),2,'omitnan');
    PARAMETERS_HZ(s).A_variance = var(PARAMETERS_HZ(s).A_HZ(:,~PARAMETERS_HZ(s).A_outlr),0,2,'omitnan');
end

%%%%%%%% CALCULATE VOCAL TREND, VOCAL ERROR AND LATENT VOCAL TARGET %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%IN HERTZ
%LATENT_VOCAL_TARGET = Is the reference frequency (ref_freq)
%VOCAL_TREND =  Mean - ref_freq
%VOCAL_NOISE = Mean value subtracted from the raw values (16 values per
%timebin) ==> mean of Vocal noise = 0
%VOCAL_ERROR = Variance

for s = 1:nsubj
    PARAMETERS_HZ(s).A_latentVocalTarget = ref_freq(s);
    PARAMETERS_HZ(s).A_vocalTrend = PARAMETERS_HZ(s).A_mean - ref_freq(s);
    PARAMETERS_HZ(s).A_vocalNoise = PARAMETERS_HZ(s).A_HZ - PARAMETERS_HZ(s).A_mean;
    PARAMETERS_HZ(s).A_vocalNoise(:,PARAMETERS_HZ(s).A_outlr) = nan; %Set the corresponding noise terms of the outliers to NaN (otherwise one gets a size-conflict for the cell arrays in TestB, C and D since higher dimensional cellarrays are used there.
    PARAMETERS_HZ(s).A_vocalError = PARAMETERS_HZ(s).A_variance;
end


disp("Extracted Parameters for Test A...");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Sort data, such that the third dimension corresponds to the two different targets (idx = 0 for Target 0 Cent & idx = 1 for Target 400 Cent).


for s = 1:nsubj
    for t = 1:length(TARGET)
        idx = testB(s).target == TARGET(t);
        PARAMETERS_HZ(s).B_HZ(:,:,t) = testB(s).raw(:,idx);
    end
end

%%% Make sure, that at least 10 trials are present at each timepoint per target. Otherwise the mean and var have no power.
for s = 1:nsubj
    for t = 1:length(TARGET)
        for i = 1:LENGTH
            if sum(~isnan( PARAMETERS_HZ(s).B_HZ(i,:,t))) < 10
                PARAMETERS_HZ(s).B_HZ(i:end,:,t) = nan;
                continue
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% REMOVE UNNATURAL JUMPS %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Unnatural jumps in pitch are removed (Can happen at the end of a trial).

%ManualRemovalOutliers_TestB;


%%%%%%%%%%%%%%%%%%%%%%% DELETE OUTLIERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Detect Outliers by using the isoutlier() matlab function for five different time Windows.
%%% The time windows are specified in the Main.m File.
%%% A trial is flagged as an outlier if it is detected in two or more of the time windows.

for s = 1:nsubj
    for t = 1:length(TARGET)
        outlr1 = isoutlier(mean(PARAMETERS_HZ(s).B_HZ(WINDOW_OUTLIER_1,:,t),1,'omitnan'),'median');
        outlr2 = isoutlier(mean(PARAMETERS_HZ(s).B_HZ(WINDOW_OUTLIER_2,:,t),1,'omitnan'),'median');
        outlr3 = isoutlier(mean(PARAMETERS_HZ(s).B_HZ(WINDOW_OUTLIER_3,:,t),1,'omitnan'),'median');
        outlr4 = isoutlier(mean(PARAMETERS_HZ(s).B_HZ(WINDOW_OUTLIER_4,:,t),1,'omitnan'),'median');
        outlr5 = isoutlier(mean(PARAMETERS_HZ(s).B_HZ(WINDOW_OUTLIER_5,:,t),1,'omitnan'),'median');
        outlr_sum = outlr1 + outlr2 + outlr3 + outlr4 + outlr5;
        outlr = outlr_sum >= 2;
        PARAMETERS_HZ(s).B_outlr(:,:,t) = outlr;
        PARAMETERS_HZ(s).B_nOutlr(:,:,t) = sum(PARAMETERS_HZ(s).B_outlr(:,:,t));
    end
end



%%%%%%%%%%%%%%%%%% CALCULATE MEAN AND VARIANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calcualte MEAN and VARIANCE for testB and append to parameters File. Outliers are excluded in calculation.

for s = 1:nsubj
    for t = 1:length(TARGET)
        PARAMETERS_HZ(s).B_mean(:,t) = mean(PARAMETERS_HZ(s).B_HZ(:,~PARAMETERS_HZ(s).B_outlr(:,:,t),t),2,'omitnan');
        PARAMETERS_HZ(s).B_variance(:,t) = var(PARAMETERS_HZ(s).B_HZ(:,~PARAMETERS_HZ(s).B_outlr(:,:,t),t),0,2,'omitnan');
    end
end    


%%%%%%%%%%%%%%%%%%%CALCULATE PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% B_sensoryTargetEffect = Mean - target - vocal trend. The vocal trend is only measured without a target in Test A. So for both targets in test B, the same vocal trend is used.
%%% B_vocalTargetNoise = Mean value subtracted from the raw values (16 values per timebin and target).
%%% B_vocalTargetError = Variance of the vocalTargetNoise.


for s = 1:nsubj
    for t = 1:length(TARGET)
        PARAMETERS_HZ(s).B_sensoryTargetEffect(:,t) = PARAMETERS_HZ(s).B_mean(:,t) - PARAMETERS_HZ(s).TARGET(t) - PARAMETERS_HZ(s).A_vocalTrend(:,1); % assumption that the vocal trend is same for 0 and 400 CENT
        PARAMETERS_HZ(s).B_vocalTargetNoise(:,:,t) = PARAMETERS_HZ(s).B_HZ(:,:,t) - PARAMETERS_HZ(s).B_mean(:,t); % Cannot only append vocalTargetNoise data for trials that are not outliers because this will result in a size confilct.
        % Have to set outliers to nan instead to always have an 200x16 sized cell array.
        PARAMETERS_HZ(s).B_vocalTargetNoise(:,PARAMETERS_HZ(s).B_outlr(:,:,t),t) = nan; 
        PARAMETERS_HZ(s).B_vocalTargetError(:,t) = var(PARAMETERS_HZ(s).B_vocalTargetNoise(:,:,t),0,2,'omitnan');
    end
end        
        

disp("Extracted Parameters for Test B in _HZs...");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test CA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Sort data, such that the third dimension corresponds to the two different targets (idx = 1 for Target 0 Cent & idx = 2 for Target 400 Cent).

for s=1:nsubj
    for t = 1:length(TARGET)
        PARAMETERS_HZ(s).CA_HZ(:,:,t) = testC(s).rawA(:,find(testC(s).targetA == TARGET(t)));
    end
end

%%% Make sure, that at least 20 trials are present at each timepoint per target. Otherwise the mean and var have no power.

for s = 1:nsubj
    for t = 1:length(TARGET)
        for i = 1:LENGTH
            if sum(~isnan( PARAMETERS_HZ(s).CA_HZ(i,:,t))) < 20
                PARAMETERS_HZ(s).CA_HZ(i:end,:,t) = nan;
                continue
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% REMOVE UNNATURAL JUMPS %%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Unnatural jumps in pitch are removed (Can happen at the end of a trial).

%ManualRemovalOutliers_TestCA;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELETE OUTLIERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Detect Outliers by using the isoutlier() matlab function for five different time Windows.
%%% The time windows are specified in the Main.m File.
%%% A trial is flagged as an outlier if it is detected in two or more of the time windows.

for s = 1:nsubj
    for t = 1:length(TARGET)
        outlr1 = isoutlier(mean(PARAMETERS_HZ(s).CA_HZ(WINDOW_OUTLIER_1,:,t),1,'omitnan'),'median');
        outlr2 = isoutlier(mean(PARAMETERS_HZ(s).CA_HZ(WINDOW_OUTLIER_2,:,t),1,'omitnan'),'median');
        outlr3 = isoutlier(mean(PARAMETERS_HZ(s).CA_HZ(WINDOW_OUTLIER_3,:,t),1,'omitnan'),'median');
        outlr4 = isoutlier(mean(PARAMETERS_HZ(s).CA_HZ(WINDOW_OUTLIER_4,:,t),1,'omitnan'),'median');
        outlr5 = isoutlier(mean(PARAMETERS_HZ(s).CA_HZ(WINDOW_OUTLIER_5,:,t),1,'omitnan'),'median');
        outlr_sum = outlr1 + outlr2 + outlr3 + outlr4 + outlr5;
        outlr = outlr_sum >= 2;
        PARAMETERS_HZ(s).CA_outlr(:,:,t) = outlr;
        PARAMETERS_HZ(s).CA_nOutlr(:,:,t) = sum(PARAMETERS_HZ(s).CA_outlr(:,:,t));
    end
end


%%%%%%%%%%%%%%% CALCULATE MEAN, VARIANCE & PARAMETERS%%%%%%%%%%%%%%%%%%%%%%

%%% Calcualte MEAN and VARIANCE for testCA and append to parameters File.  Outliers are excluded in calculation.
%%% CA_vocalFeebackEffect = Mean - target - B_sensoryTargetEffect - A_vocalTrend. The vocal trend is only measured without a target in Test A. So for both targets in test CA, the same vocal trend is used (1 value per timebin and target).
%%% CA_feedbackNoise = Mean value subtracted from the raw values (70 values per timebin and target).
%%% CA_feedbackError = Variance of the CA_feedbackNoise (1 value per timebin and target).

for s=1:nsubj
    for t = 1:length(TARGET)
        PARAMETERS_HZ(s).CA_mean(:,t) = mean(PARAMETERS_HZ(s).CA_HZ(:,~PARAMETERS_HZ(s).CA_outlr(:,:,t),t),2,'omitnan');
        PARAMETERS_HZ(s).CA_variance(:,t) =  var(PARAMETERS_HZ(s).CA_HZ(:,~PARAMETERS_HZ(s).CA_outlr(:,:,t),t),0,2,'omitnan');
        PARAMETERS_HZ(s).CA_vocalFeebackEffect(:,t) =  PARAMETERS_HZ(s).CA_mean(:,t) - PARAMETERS_HZ(s).TARGET(t) - PARAMETERS_HZ(s).A_vocalTrend(:,1) - PARAMETERS_HZ(s).B_sensoryTargetEffect(:,t);
        PARAMETERS_HZ(s).CA_feedbackNoise(:,:,t) = PARAMETERS_HZ(s).CA_HZ(:,:,t) - PARAMETERS_HZ(s).CA_mean(:,t);
        PARAMETERS_HZ(s).CA_feedbackNoise(:,PARAMETERS_HZ(s).CA_outlr(:,:,t),t) = nan; %Have to set outliers to nan instead to always have an 200x16 sized cell array.
        PARAMETERS_HZ(s).CA_feedbackError(:,t) = var(PARAMETERS_HZ(s).CA_feedbackNoise(:,:,t),0,2,'omitnan');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test CB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Sort data, such that the third dimension corresponds to the two different targets (idx = 1 for Target 0 Cent & idx = 2 for Target 400 Cent).
%%% And the fourth dimension to the different SHIFTS (from idx = 1 to idx = 7 in SHIFT array [-100   -50   -25     0    25    50   100]).
%%% Convert dada to cent using the participants respective ref_freq.

for s=1:nsubj
    for t = 1:length(TARGET)
        for f = 1:length(SHIFT)
            PARAMETERS_HZ(s).CB_HZ(:,:,t,f) = testC(s).rawB(:,find(testC(s).targetB == TARGET(t) & testC(s).shift == SHIFT(f)));  
        end
    end

end

%%% Make sure, that at least 7 trials are present at each timepoint per target and shift. Otherwise the mean and var have no power.
for s = 1:nsubj
    for t = 1:length(TARGET)
        for f = 1:length(SHIFT)
            for i = 1:LENGTH
                if sum(~isnan( PARAMETERS_HZ(s).CB_HZ(i,:,t,f))) < 7
                    PARAMETERS_HZ(s).CB_HZ(i:end,:,t,f) = nan;
                    continue
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% REMOVE UNNATURAL JUMPS %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Unnatural jumps in pitch are removed (Can happen at the end of a trial).

%ManualRemovalOutliers_TestCB;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELETE OUTLIERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Detect Outliers by using the isoutlier() matlab function for five different time Windows.
%%% The time windows are specified in the Main.m File.
%%% A trial is flagged as an outlier if it is detected in two or more of the time windows.

for s = 1:nsubj
    for t = 1:length(TARGET)
        for f = 1:length(SHIFT)   
            outlr1 = isoutlier(mean(PARAMETERS_HZ(s).CB_HZ(WINDOW_OUTLIER_1,:,t,f),1,'omitnan'),'median');
            outlr2 = isoutlier(mean(PARAMETERS_HZ(s).CB_HZ(WINDOW_OUTLIER_2,:,t,f),1,'omitnan'),'median');
            outlr3 = isoutlier(mean(PARAMETERS_HZ(s).CB_HZ(WINDOW_OUTLIER_3,:,t,f),1,'omitnan'),'median');
            outlr4 = isoutlier(mean(PARAMETERS_HZ(s).CB_HZ(WINDOW_OUTLIER_4,:,t,f),1,'omitnan'),'median');
            outlr5 = isoutlier(mean(PARAMETERS_HZ(s).CB_HZ(WINDOW_OUTLIER_5,:,t,f),1,'omitnan'),'median');
            outlr_sum = outlr1 + outlr2 + outlr3 + outlr4 + outlr5;
            outlr = outlr_sum >= 2;
            PARAMETERS_HZ(s).CB_outlr(:,:,t,f) = outlr;
            PARAMETERS_HZ(s).CB_nOutlr(:,:,t,f) = sum(PARAMETERS_HZ(s).CB_outlr(:,:,t,f));
        end
    end
end


%%%%%%%%%%%%%%%%% CALCULATE MEAN, VARIANCE & PARAMETERS %%%%%%%%%%%%%%%%%%%

%%% Calcualte MEAN and VARIANCE for testC and append to parameters File. Outliers are excluded in calculation.
%%% CB_shiftedFeedbackEffect = Mean - target - CA_vocalFeebackEffect - B_sensoryTargetEffect - A_vocalTrend. The vocal trend is only measured without a target in Test A. So for both targets in test CB, the same vocal trend is used (1 value per timebin and target and shift).
%%% For all shift conditions the same values of CA_vocalFeebackEffect, B_sensoryTargetEffect, and A_vocalTrend are used.
%%% CB_shiftedFeedbackNoise = Mean value subtracted from the raw values (10 values per timebin and target and shift).
%%% CB_shiftedFeedbackError = Variance of the CB_shiftedFeedbackNoise (1 value per timebin and target and shift).

for s=1:nsubj
    for t = 1:length(TARGET)
         for f = 1:length(SHIFT)
            PARAMETERS_HZ(s).CB_mean(:,:,t,f) = mean(PARAMETERS_HZ(s).CB_HZ(:,~PARAMETERS_HZ(s).CB_outlr(:,:,t,f),t,f),2,'omitnan');
            PARAMETERS_HZ(s).CB_variance(:,:,t,f) = var(PARAMETERS_HZ(s).CB_HZ(:,~PARAMETERS_HZ(s).CB_outlr(:,:,t,f),t,f),0,2,'omitnan');
            PARAMETERS_HZ(s).CB_shiftedFeedbackEffect(:,:,t,f) =  PARAMETERS_HZ(s).CB_mean(:,:,t,f) - PARAMETERS_HZ(s).TARGET(t) - PARAMETERS_HZ(s).A_vocalTrend(:,1) - PARAMETERS_HZ(s).B_sensoryTargetEffect(:,t) - PARAMETERS_HZ(s).CA_vocalFeebackEffect(:,t);
            PARAMETERS_HZ(s).CB_shiftedFeedbackNoise(:,:,t,f) =  PARAMETERS_HZ(s).CB_HZ(:,:,t,f) - PARAMETERS_HZ(s).CB_mean(:,:,t,f);
            PARAMETERS_HZ(s).CB_shiftedFeedbackNoise(:,PARAMETERS_HZ(s).CB_outlr(:,:,t,f),t,f) = nan; %Have to set outliers to nan instead to always have an 200x16 sized cell array.
            PARAMETERS_HZ(s).CB_shiftedFeedbackError(:,:,t,f) = var(PARAMETERS_HZ(s).CB_shiftedFeedbackNoise(:,:,t,f),0,2,'omitnan');
         end
    end
end



disp("Extracted Parameters for Test C in _HZ...");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA CLEANUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%%% Correct a wrong assignment of a shift condition for participant f04

f04_idx = find(sbjnames == "f04");
if ~isempty(f04_idx)
    testD(f04_idx).shift(116) = 100;
end

%%% Sort data, such that the third dimension corresponds to the two different targets (idx = 1 for Target 0 Cent & idx = 2 for Target 400 Cent).
%%% And the fourth dimension to the different SHIFTS (from idx = 1 to idx = 7 in SHIFT array [-100   -50   -25     0    25    50   100]).
%%% The data in test D is in the units CENT.

for s=1:nsubj
    for t=1:length(TARGET)
        for f=1:length(SHIFT)
            idx = find(testD(s).target == TARGET(t) & testD(s).shift == SHIFT(f));    
            PARAMETERS_HZ(s).D_cent(:,:,t,f) = testD(s).tuning(:,idx) + TARGET(t) + SHIFT(f); %Stores the frequency (sncl target and shift) in CENT.
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERT CENT TO HZ %%%%%%%%%%%%%%%%%%%%%%%%%

for s=1:nsubj
    for t=1:length(TARGET)
        for f=1:length(SHIFT)
            PARAMETERS_HZ(s).D_HZ(:,:,t,f) =2.^(PARAMETERS_HZ(s).D_cent(:,:,t,f) ./ 1200) * ref_freq(s);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELETE OUTLIERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Detect Outliers by using the isoutlier() matlab function and one time window.
%%% The time window is specified in the Main.m script (WINDOW_OUTLIER_D).
%%% A trial is flagged as an outlier if it is detected in the time window.


for s = 1:nsubj
    for t = 1:length(TARGET)
        for f = 1:length(SHIFT)   
            outlr = isoutlier(mean(PARAMETERS_HZ(s).D_HZ(WINDOW_OUTLIER_D,:,t,f),1,'omitnan'),'median');
            PARAMETERS_HZ(s).D_outlr(:,:,t,f) = outlr;
            PARAMETERS_HZ(s).D_nOutlr(:,:,t,f) = sum(PARAMETERS_HZ(s).D_outlr(:,:,t,f));
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%% CALCULATE MEAN, VARIANCE & PARAMETERS %%%%%%%%%%%%

%%% Calcualte MEAN and VARIANCE for testD and append to parameters File. Outliers are excluded in calculation.
%%% D_manualFeedback: Mean(tEnd) - target.
%%% D_meanEffectOfInitialCondition = Mean(t) - manualFeedback / mean(t=0) - manualFeedback ==> gives a number between 1 and 0
%%% D_manualFeedbackNoise = Rawdata - Mean
%%% D_manualFeedbackError = variance of manualFeedbackNoise.

for s = 1:nsubj
    for t=1:length(TARGET)
        for f=1:length(SHIFT)
            PARAMETERS_HZ(s).D_mean(:,:,t,f) = mean(PARAMETERS_HZ(s).D_HZ(:,~PARAMETERS_HZ(s).D_outlr(:,:,t,f),t,f),2,'omitnan');
            PARAMETERS_HZ(s).D_variance(:,:,t,f) = var(PARAMETERS_HZ(s).D_HZ(:,~PARAMETERS_HZ(s).D_outlr(:,:,t,f),t,f),0,2,'omitnan');
            PARAMETERS_HZ(s).D_manualFeedback(t,f) = PARAMETERS_HZ(s).D_mean(end,:,t,f) - PARAMETERS_HZ(s).TARGET(t);
            diff = PARAMETERS_HZ(s).D_mean(1,:,t,f) - PARAMETERS_HZ(s).D_manualFeedback(t,f);
            PARAMETERS_HZ(s).D_meanEffectOfInitialCondition(:,:,t,f) = ( PARAMETERS_HZ(s).D_mean(:,:,t,f) - PARAMETERS_HZ(s).D_manualFeedback(t,f) ) / diff;
            PARAMETERS_HZ(s).D_manualFeedbackNoise(:,:,t,f) = PARAMETERS_HZ(s).D_HZ(:,:,t,f) - PARAMETERS_HZ(s).D_mean(:,:,t,f);
            PARAMETERS_HZ(s).D_manualFeedbackNoise(:,PARAMETERS_HZ(s).D_outlr(:,:,t,f),t,f) = nan;
            PARAMETERS_HZ(s).D_manualFeedbackError(:,:,t,f) = var(PARAMETERS_HZ(s).D_manualFeedbackNoise(:,:,t,f),0,2,'omitnan');
        end
    end
end

disp("Extracted Parameters for Test D in Hertz...");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Save PARMAETERS_HZ.mat file at the choosen location (Main.m script)
%%% Clear all variables and load the needed Parameters again.

filename = savepath + "PARAMETERS_HZ.mat";
save(filename,'PARAMETERS_HZ');
disp("Parameter in HZ saved at:");
disp(filename);
clear all;
getPath;
load(savepath + "variables.mat");
load(savepath + "PARAMETERS_HZ.mat");
if exist(savepath + "PARAMETERS_CENT.mat", 'file')
        load(savepath + "PARAMETERS_CENT.mat");
end

