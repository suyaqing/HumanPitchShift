% Fit a LME to all data at each time point. Different shifts and
% targets are modeled as separate effects. Pitch trend is the  constant in
% the model. Subject numbers as a random effect

savepath = 'D:\Pitch adaptation\test\';
file_CENT = savepath + "PARAMETERS_CENT_newOn.mat";
file_var = savepath + "variables.mat";
load(file_CENT)
load(file_var)
%EXCLUDE = {'f12NC', 'f18NC','f19NC','f25NC','m11NC','m12NC','m18NC','m23','m27'};
EXCLUDE = {'NONE'};
% construct input table for the GLM, one subject after another
varNames = {'pitch', 'T0', 'T400', 'F0', 'F400', 'S100n', 'S50n', 'S25n', ...
    'S25', 'S50', 'S100', 'Subj'};
nvar = length(varNames);
tA = 16;
tB = 16;
tCA = 70;
tCB = 10;
offset = 100; % number of data points before t=0
% try a single time point here to start
% t = 100;
%%
tic
excl_otlr = 0; % whether to exclude outliers
data_alltime = [];
glms = struct([]);
% since the majority of pitch traces don't start at 0 need to find a way to
% eliminate time points without sufficient data, otherwise the design
% matrix will not have full column rank
for t = (1:120)+offset
    data_allsub = [];
    glme = [];
    for s = 1:nsubj
        data = PARAMETERS_CENT(s);
        if ~any(ismember(EXCLUDE, data.Sbjnames))
            if excl_otlr
                data_A = zeros(tA-sum(data.A_outlr), nvar);
                data_A(:, 1) = data.A_cent(t, ~data.A_outlr)';
       
                data_B0 = zeros(tB-sum(data.B_outlr(:, :, 1)), nvar);
                data_B0(:, 1) = data.B_cent(t, ~data.B_outlr(:, :, 1), 1)';
        
                data_B400 = zeros(tB-sum(data.B_outlr(:, :, 2)), nvar);
                data_B400(:, 1) = data.B_cent(t, ~data.B_outlr(:, :, 2), 2)'-400;
       
                data_CA0 = zeros(tCA-sum(data.CA_outlr(:, :, 1)), nvar);
                data_CA0(:, 1) = data.CA_cent(t, ~data.CA_outlr(:, :, 1), 1)';
        
                data_CA400 = zeros(tCA-sum(data.CA_outlr(:, :, 2)), nvar);
                data_CA400(:, 1) = data.CA_cent(t, ~data.CA_outlr(:, :, 2), 2)'-400;
            else
                data_A = zeros(tA, nvar);
                data_A(:, 1) = data.A_cent(t, :)';
       
                data_B0 = zeros(tB, nvar);
                data_B0(:, 1) = data.B_cent(t, :, 1)';
        
                data_B400 = zeros(tB, nvar);
                data_B400(:, 1) = data.B_cent(t, :, 2)'-400;
                       
                data_CA0 = zeros(tCA, nvar);
                data_CA0(:, 1) = data.CA_cent(t, :, 1)';
                         
                data_CA400 = zeros(tCA, nvar);
                data_CA400(:, 1) = data.CA_cent(t, :, 2)'-400;
                 
            end
            data_B0(:, 2) = 1; 
            data_B400(:, 3) = 1; 
            data_CA0(:, [2,4]) = 1;
            data_CA400(:, [3,5]) = 1;

            data_CB = [];
            for target = 1:2
                for shift = 1:7
                    if excl_otlr
                        data_temp = zeros(tCB - sum(data.CB_outlr(:, :, target, shift)), nvar);
                        data_temp(:, 1) = data.CB_cent(t, ~data.CB_outlr(:, :, target, shift), target, shift)';
                    else
                        data_temp = zeros(tCB, nvar);
                        data_temp(:, 1) = data.CB_cent(t, :, target, shift)';
                    end
                    if target == 2
                        data_temp(:, 1) = data_temp(:, 1) - 400;
                    end
                    % target
                    data_temp(:, target + 1) = 1;
                    % feedback
                    data_temp(:, target + 3) = 1;
                    % shift in feedback, assume independent from effect of feedback
                    if shift < 4
                        data_temp(:, shift + 5) = 1;
                    else
                        if shift > 4
                            data_temp(:, shift + 4) = 1;
                        end
                    end
                    data_CB = [data_CB; data_temp];
                end
            end
            data_subj = [data_A; data_B0; data_B400; data_CA0; data_CA400; data_CB];
            data_subj(:, end) = s;
            data_allsub = [data_allsub; data_subj];
        end
    end
    data_alltime(:, :, t-offset) = data_allsub;
    Tbl = table(data_allsub(:, 1), data_allsub(:, 2), data_allsub(:, 3), data_allsub(:, 4), ...
    data_allsub(:, 5), data_allsub(:, 6), data_allsub(:, 7), data_allsub(:, 8), ...
    data_allsub(:, 9), data_allsub(:, 10), data_allsub(:, 11), data_allsub(:, 12), ...
    'VariableNames', varNames);
    try
        glme = fitglme(Tbl, ...
        'pitch ~ 1 + T0 + T400 + F0 + F400 + S100n + S50n + S25n + S25 + S50 + S100 + (1|Subj)');
        glms(t-offset).glme = glme;
%         [psi, ~,stats] = covarianceParameters(glme);
%         glms(t-offset).coeff = glme.Coefficients;
%         glms(t-offset).cov = glme.CoefficientCovariance;
%         glms(t-offset).rsq = glme.Rsquared;
%         glms(t-offset).crit = glme.ModelCriterion;
%         glms(t-offset).random = stats{1};
%         glms(t-offset).psi = psi;
    catch
        glms(t-offset).glme = glme;
%         glms(t-offset).coeff = [];
%         glms(t-offset).cov = [];
%         glms(t-offset).rsq = [];
%         glms(t-offset).crit = [];
%         glms(t-offset).random = [];
    end
end
toc
%%
% Tbl = table(data_allsub(:, 1), data_allsub(:, 2), data_allsub(:, 3), data_allsub(:, 4), ...
%     data_allsub(:, 5), data_allsub(:, 6), data_allsub(:, 7), data_allsub(:, 8), ...
%     data_allsub(:, 9), data_allsub(:, 10), data_allsub(:, 11), data_allsub(:, 12), ...
%     'VariableNames', varNames);
if excl_otlr
    fn = savepath + "MaxiGLMM_RandIntSubj.mat";
else
    fn = savepath + "MaxiGLMM_RandIntSubj_NoOL.mat";
end
save(fn, "glms")

% glme = fitglme(Tbl, ...
%     'pitch ~ 1 + T0 + T400 + F0 + F400 + S100n + S50n + S25n + S25 + S50 + S100 + (1|Subj)');