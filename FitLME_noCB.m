% Fit a LME to all data at each time point except for task CB. 
% Subject numbers as a random effect of intercept

savepath = 'D:\Pitch adaptation\test\';
file_CENT = savepath + "PARAMETERS_CENT.mat";
file_var = savepath + "variables.mat";
load(file_CENT)
load(file_var)
% EXCLUDE = {'f12NC', 'f18NC','f19NC','f25NC','m11NC','m12NC','m18NC','m23','m27'};
EXCLUDE = {'None'};
% construct input table for the GLM, one subject after another
varNames = {'pitch', 'T0', 'T400', 'F0', 'F400', 'Subj'};
nvar = length(varNames);
tA = 16;
tB = 16;
tCA = 70;
tCB = 10;
offset = 50; % number of data points before t=0
% try a single time point here to start
% t = 100;
%%
tic
data_alltime = [];
glms = struct([]);
% since the majority of pitch traces don't start at 0 need to find a way to
% eliminate time points without sufficient data, otherwise the design
% matrix will not have full column rank
for t = (1:170)+offset
    data_allsub = [];
    for s = 1:nsubj
        data = PARAMETERS_CENT(s);
        if ~any(ismember(EXCLUDE, data.Sbjnames))
            data_A = zeros(tA-sum(data.A_outlr), nvar);
            data_A(:, 1) = data.A_cent(t, ~data.A_outlr)';
    %         data_A(:, end) = s;
    
            data_B0 = zeros(tB-sum(data.B_outlr(:, :, 1)), nvar);
            data_B0(:, 1) = data.B_cent(t, ~data.B_outlr(:, :, 1), 1)';
            data_B0(:, 2) = 1; 
    %         data_B0(:, end) = s;
    
            data_B400 = zeros(tB-sum(data.B_outlr(:, :, 2)), nvar);
            data_B400(:, 1) = data.B_cent(t, ~data.B_outlr(:, :, 2), 2)'-400;
            data_B400(:, 3) = 1; 
    %         data_B400(:, end) = s;
    
            data_CA0 = zeros(tCA-sum(data.CA_outlr(:, :, 1)), nvar);
            data_CA0(:, 1) = data.CA_cent(t, ~data.CA_outlr(:, :, 1), 1)';
            data_CA0(:, [2,4]) = 1; 
    %         data_CA0(:, end) = s;
    
            data_CA400 = zeros(tCA-sum(data.CA_outlr(:, :, 2)), nvar);
            data_CA400(:, 1) = data.CA_cent(t, ~data.CA_outlr(:, :, 2), 2)'-400;
            data_CA400(:, [3,5]) = 1; 
    %         data_CA400(:, end) = s;
            
            data_subj = [data_A; data_B0; data_B400; data_CA0; data_CA400];
            data_subj(:, end) = s;
            data_allsub = [data_allsub; data_subj];
        end
    end
    data_alltime(:, :, t-offset) = data_allsub;
    Tbl = table(data_allsub(:, 1), data_allsub(:, 2), data_allsub(:, 3), data_allsub(:, 4), ...
    data_allsub(:, 5), data_allsub(:, 6), 'VariableNames', varNames);
    try
        glme = fitglme(Tbl, ...
        'pitch ~ 1 + T0 + T400 + F0 + F400 + (1|Subj)');
        glms(t-offset).glme = glme;

    catch
        glms(t-offset).glme = [];

    end
end
toc
%%
% Tbl = table(data_allsub(:, 1), data_allsub(:, 2), data_allsub(:, 3), data_allsub(:, 4), ...
%     data_allsub(:, 5), data_allsub(:, 6), data_allsub(:, 7), data_allsub(:, 8), ...
%     data_allsub(:, 9), data_allsub(:, 10), data_allsub(:, 11), data_allsub(:, 12), ...
%     'VariableNames', varNames);

fn = savepath + "MaxiGLMMNoCB_RandIntSubj_ALLSUB.mat";
save(fn, "glms", "offset")

% glme = fitglme(Tbl, ...
%     'pitch ~ 1 + T0 + T400 + F0 + F400 + S100n + S50n + S25n + S25 + S50 + S100 + (1|Subj)');