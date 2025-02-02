% compare the performance of task C and task B using stats from fitted
% curves
clear
savepath = 'D:\Pitch adaptation\test\';
GLMMdata = savepath + "MaxiGLMM_RandIntSubj.mat";
Sigmoidfit = 'taskC_sigmoid_fit_newWin.mat';
file_CENT = savepath + "PARAMETERS_CENT_newOn.mat";
file_var = savepath + "variables.mat";
load(file_CENT)
load(file_var)
load(GLMMdata)
load(Sigmoidfit)
% EXCLUDE = {'f12NC', 'f18NC','f19NC','f25NC','m11NC','m12NC','m18NC','m23','m27'};
EXCLUDE = {'None'};

folder = savepath + "\Plots\GLMM\";
RESOLUTION  = 300; %Resolution in dots per inch (DPI), specified as a whole number that is greater than or equal to 1.
FORMAT = ".pdf"; %Possible formats are: '.jpg', '.png', '.tif', '.gif', '.pdf', '.emf' or '.eps'
CONTENT_TYPE = 'auto'; 

col_A = [0 0.4470 0.7410];
col_B = [[0 0.4470 0.7410];[0.4940 0.1840 0.5560];]; 
col_CB = [[0.4940 0.1840 0.5560];[0 0.4470 0.7410];[0.3010 0.7450 0.9330];[0.4660 0.6740 0.1880];[0.9290 0.6940 0.1250];[0.8500 0.3250 0.0980];[0.6350 0.0780 0.1840];]; %RGB color code, nice for 7 shifts.
fmt = '%.2e';
offset = 0;
%%
glms_stable = glms(offset+51:offset+130);
T = length(glms_stable);
% tt = ([1:T]-(100-offset))/100 - 1/100; tt = tt(:);
Nparam = 11;
%nsub = length(stats_all_sub);
% taskB_ini = zeros(2, nsub); taskB_stab = zeros(2, nsub);
% taskCA_ini = zeros(2, nsub); taskCA_stab = zeros(2, nsub);
coefs_all = nan(Nparam, T);
rand_interc = nan(nsubj, T);
rsqC = rsq; % squared bias in task C sigmoid fit
slopeC = abs(log(-slope)); % how close the slope of task C sigmoid is to 1
slopeC_norm = slope./stds;
rmseB = nan(size(rmse));
bias2B = nan(size(rsqC));
rmseB_targ = nan(size(rsqC));
rmseCA_targ = nan(size(rsqC));
meanCA = nan(size(rsqC));

for t = 1:T
    % fisrt get fixed and random effect coefficients
    glm = glms_stable(t).glme;
    if ~isempty(glm)
        coefs_all(:, t) = glm.Coefficients.Estimate;
        [B, BNames, ~] = randomEffects(glm);
        subj_idx = str2double(table2array(BNames(:, 2)));
        rand_interc(subj_idx, t) = B;
    end
end

subj_idx = str2double(table2array(BNames(:, 2)));
rand_interc = rand_interc';
betas = coefs_all';
indT = (100+50)+(1:T); % note there's a 100ms offset in data
pA = nan(nsubj, 1); pB = nan(nsubj, 2);
% pCA = nan(nsubj, 2); pCB = nan(nsubj, 2, 6);
for s = 1:nsubj
    data = PARAMETERS_CENT(subj_idx(s)); 
    sname = data.Sbjnames;
    if ~isempty(sname)
%         sv = sv+1; %valid subj
        % first get all betas into a matrix       
        itc = rand_interc(:, subj_idx(s));
        mean_B = data.B_mean(indT, :); 
        data_B{1} = data.B_cent(indT, ~data.B_outlr(:, :, 1), 1); % time*trial
        data_B{2} = data.B_cent(indT, ~data.B_outlr(:, :, 1), 2);
        data_CA{1} = data.CA_cent(indT, ~data.CA_outlr(:, :, 1), 1); % time*trial
        data_CA{2} = data.CA_cent(indT, ~data.CA_outlr(:, :, 1), 2);
        meanCA(s, 1) = mean(mean(data_CA{1}));
        meanCA(s, 2) = mean(mean(data_CA{2})) - 400;
        betas_B = zeros(T, 2); %pB = zeros(1,2);
        for targ = 1:2
            betas_B(:, targ) = betas(:, 1) + betas(:, targ+1) + itc;
            
            if targ == 2
                betas_B(:, targ) = betas_B(:, targ) + 400;
            end
%             [~, pB(s, targ)] = kstest2(betas_B(:, targ), mean_B(:, targ));
            [rmseB(s, targ), bias2B(s, targ), rmseB_targ(s, targ)] = individual_rmse_task(betas_B(:, targ), data_B{targ}, targ);
            [~, ~, rmseCA_targ(s, targ)] = individual_rmse_task(betas_B(:, targ), data_CA{targ}, targ);
        end
    end
end
%%
figure
subplot(221)
hold on
plot(rsqC(:, 1), rmseB_targ(:, 1), "color", col_B(1, :), "Marker", 'o', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
plot(rsqC(:, 2), rmseB_targ(:, 2), "color", col_B(2, :), "Marker", '*', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
[r0,p0] = corr(rsqC(:, 1), rmseB_targ(:, 1),'Type','Spearman', 'Rows','complete');
[r400,p400] = corr(rsqC(:, 2), rmseB_targ(:, 2),'Type','Spearman','Rows','complete');
legend(['0 Cent, r=' num2str(r0, fmt) ', p=' num2str(p0, fmt)], ...
    ['400 Cent, r=' num2str(r400, fmt) ', p=' num2str(p400, fmt)], 'FontSize',12)
% xlim([0, lim])
% ylim([0, lim])
ax = gca; 
ax.FontSize = 12; 
xlabel("Task C R squared", 'FontSize',14)
ylabel("Task B RMS error to target", 'FontSize',14)
% title("RMS error of fit", 'FontSize',14)

subplot(222)
hold on
plot(dprime(:, 1), rmseB_targ(:, 1), "color", col_B(1, :), "Marker", 'o', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
plot(dprime(:, 2), rmseB_targ(:, 2), "color", col_B(2, :), "Marker", '*', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
[r0,p0] = corr(dprime(:, 1), rmseB_targ(:, 1),'Type','Spearman', 'Rows','complete');
[r400,p400] = corr(dprime(:, 2), rmseB_targ(:, 2),'Type','Spearman','Rows','complete');
legend(['0 Cent, r=' num2str(r0, fmt) ', p=' num2str(p0, fmt)], ...
    ['400 Cent, r=' num2str(r400, fmt) ', p=' num2str(p400, fmt)], 'FontSize',12)
% xlim([0, lim])
% ylim([0, lim])
xlabel("Task C d-prime", 'FontSize',14)
ylabel("Task B RMS error to target", 'FontSize',14)
ax = gca; 
ax.FontSize = 12; 
% title("Squared bias of fit", 'FontSize',14)
        
subplot(223)
hold on
plot(slopeC_norm(:, 1), rmseB_targ(:, 1), "color", col_B(1, :), "Marker", 'o', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
plot(slopeC_norm(:, 2), rmseB_targ(:, 2), "color", col_B(2, :), "Marker", '*', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
[r0,p0] = corr(slopeC_norm(:, 1), rmseB_targ(:, 1),'Type','Spearman', 'Rows','complete');
[r400,p400] = corr(slopeC_norm(:, 2), rmseB_targ(:, 2),'Type','Spearman','Rows','complete');
legend(['0 Cent, r=' num2str(r0, fmt) ', p=' num2str(p0, fmt)], ...
    ['400 Cent, r=' num2str(r400, fmt) ', p=' num2str(p400, fmt)], 'FontSize',12)
% legend('0 Cent', '400 Cent',  'FontSize',12)
% xlim([0, lim])
% ylim([0, 8])
ax = gca; 
ax.FontSize = 12; 
xlabel("Task C normalized slope", 'FontSize',14)
ylabel("Task B RMS error to target", 'FontSize',14)
% title("Task C performance", 'FontSize',14)


subplot(224)
hold on
plot(slopeC(:, 1), rmseB_targ(:, 1), "color", col_B(1, :), "Marker", 'o', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
plot(slopeC(:, 2), rmseB_targ(:, 2), "color", col_B(2, :), "Marker", '*', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
[r0,p0] = corr(slopeC(:, 1), rmseB_targ(:, 1),'Type','Spearman', 'Rows','complete');
[r400,p400] = corr(slopeC(:, 2), rmseB_targ(:, 2),'Type','Spearman','Rows','complete');
legend(['0 Cent, r=' num2str(r0, fmt) ', p=' num2str(p0, fmt)], ...
    ['400 Cent, r=' num2str(r400, fmt) ', p=' num2str(p400, fmt)], 'FontSize',12)
% legend('0 Cent', '400 Cent',  'FontSize',12)
% xlim([0, lim])
% ylim([0, 8])
ax = gca; 
ax.FontSize = 12; 
xlabel("Task C log slope", 'FontSize',14)
ylabel("Task B RMS error to target", 'FontSize',14)
% title("Task C performance", 'FontSize',14)

% suptitle('Task B performance vs Task C controller')
%%

figure;
subplot(131)
hold on
plot(interc(:, 1), meanCA(:, 1), "color", col_B(1, :), "Marker", 'o', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
plot(interc(:, 2), meanCA(:, 2), "color", col_B(2, :), "Marker", '*', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
[r0,p0] = corr(interc(:, 1), meanCA(:, 1),'Type','Pearson', 'Rows','complete');
[r400,p400] = corr(interc(:, 2), meanCA(:, 2),'Type','Pearson','Rows','complete');
plot([-600, 100], [-600, 100], 'k--')
legend(['0 Cent, r=' num2str(r0, fmt) ', p=' num2str(p0, fmt)], ...
    ['400 Cent, r=' num2str(r400, fmt) ', p=' num2str(p400, fmt)], 'FontSize',12)
% legend('0 Cent', '400 Cent',  'FontSize',12)
xlim([-600, 100])
ylim([-600, 100])
ax = gca; 
ax.FontSize = 12; 
xlabel("Task C intercept", 'FontSize',14)
ylabel("Task CA mean", 'FontSize',14)

subplot(132)
hold on
plot(abs(interc(:, 1)), rmseCA_targ(:, 1), "color", col_B(1, :), "Marker", 'o', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
plot(abs(interc(:, 2)), rmseCA_targ(:, 2), "color", col_B(2, :), "Marker", '*', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
[r0,p0] = corr(abs(interc(:, 1)), rmseCA_targ(:, 1),'Type','Pearson', 'Rows','complete');
[r400,p400] = corr(abs(interc(:, 2)), rmseCA_targ(:, 2),'Type','Pearson','Rows','complete');
plot([0, 600], [0, 600], 'k--')
legend(['0 Cent, r=' num2str(r0, fmt) ', p=' num2str(p0, fmt)], ...
    ['400 Cent, r=' num2str(r400, fmt) ', p=' num2str(p400, fmt)], 'FontSize',12)
% legend('0 Cent', '400 Cent',  'FontSize',12)
% xlim([0, lim])
% ylim([0, 8])
ax = gca; 
ax.FontSize = 12; 
xlabel("Task C absolute intercept", 'FontSize',14)
ylabel("Task CA RMS error to target", 'FontSize',14)

subplot(133)
hold on
plot(abs(interc(:, 1)), rmseB_targ(:, 1), "color", col_B(1, :), "Marker", 'o', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
plot(abs(interc(:, 2)), rmseB_targ(:, 2), "color", col_B(2, :), "Marker", '*', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
[r0,p0] = corr(abs(interc(:, 1)), rmseB_targ(:, 1),'Type','Pearson', 'Rows','complete');
[r400,p400] = corr(abs(interc(:, 2)), rmseB_targ(:, 2),'Type','Pearson','Rows','complete');
legend(['0 Cent, r=' num2str(r0, fmt) ', p=' num2str(p0, fmt)], ...
    ['400 Cent, r=' num2str(r400, fmt) ', p=' num2str(p400, fmt)], 'FontSize',12)
% legend('0 Cent', '400 Cent',  'FontSize',12)
xlim([0, 600])
% ylim([0, 8])
ax = gca; 
ax.FontSize = 12; 
xlabel("Task C absolute intercept", 'FontSize',14)
ylabel("Task B RMS error to target", 'FontSize',14)

function [rmse, bias2, rmse_targ] = individual_rmse_task(fit, data, targ)
ntrial = size(data, 2);
N = size(data, 1)*ntrial;
fit_expand = fit*ones(1, ntrial);
tsse = sum(sum((data-fit_expand).^2)); % total sum of square error
rmse = sqrt(tsse/N);
if targ==1
    bias2 = mean(fit.^2);
    rmse_targ = sqrt(sum(sum(data.^2))/N);
else
    bias2 = mean((fit-400).^2);
    targ_expand = 400*ones(1, ntrial);
    rmse_targ = sqrt(sum(sum((data-targ_expand).^2))/N);
end
end
