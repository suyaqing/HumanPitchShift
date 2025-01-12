% compute and plot residual from GLME fitting
% for each time point, compute: 1) model fit, 2) mean abs residual, 3) variance
% of residual. Plot separately for each task

% last updated by Yaqing Su 11-March-2024
clear
savepath = 'D:\Pitch adaptation\test\';
GLMMdata = savepath + "MaxiGLMMnoCB_RandIntSubj.mat";
file_CENT = savepath + "PARAMETERS_CENT.mat";
file_var = savepath + "variables.mat";
load(file_CENT)
load(file_var)
load(GLMMdata)
EXCLUDE = {'f12NC', 'f18NC','f19NC','f25NC','m11NC','m12NC','m18NC','m23','m27'};

folder = savepath + "\Plots\GLMM\";
RESOLUTION  = 300; %Resolution in dots per inch (DPI), specified as a whole number that is greater than or equal to 1.
FORMAT = ".pdf"; %Possible formats are: '.jpg', '.png', '.tif', '.gif', '.pdf', '.emf' or '.eps'
CONTENT_TYPE = 'auto'; 

col_A = [0 0.4470 0.7410];
col_B = [[0 0.4470 0.7410];[0.4940 0.1840 0.5560];]; 
col_CB = [[0.4940 0.1840 0.5560];[0 0.4470 0.7410];[0.3010 0.7450 0.9330];[0.4660 0.6740 0.1880];[0.9290 0.6940 0.1250];[0.8500 0.3250 0.0980];[0.6350 0.0780 0.1840];]; %RGB color code, nice for 7 shifts.
fmt = '%.3e';
% for each participant, plot performance in task CA against performance in
% task B
%%
% definition #1 of performance: averaged upper-lower?
T = length(glms);
tt = ([1:T]-(100-offset))/100 - 1/100; tt = tt(:);
% knee = 0.4;
% Nknee = knee*100;
nsub = 42;
Nparam = 5;
%nsub = length(stats_all_sub);
% taskB_ini = zeros(2, nsub); taskB_stab = zeros(2, nsub);
% taskCA_ini = zeros(2, nsub); taskCA_stab = zeros(2, nsub);
coefs_all = nan(Nparam, T);
rand_interc = nan(nsub, T);
glm = glms(1);


taskA_fit = nan(T, 1); taskA_res = nan(T, 1);
taskB_fit = nan(T, 2); taskB_res = nan(T, 2);
taskCA_fit = nan(T, 2); taskCA_res = nan(T, 2);
% taskCB_fit = nan(T, 2, 6); taskCB_res = nan(T, 2, 6);

for t = 1:T
    glm = glms(t).glme;
    if ~isempty(glm)
        mufit = fitted(glm);
        var = glm.Variables;
        T0 = var.T0; T400 = var.T400;
        F0 = var.F0; F400 = var.F400;
%         S100 = var.S100; S50 = var.S50; S25 = var.S25;
%         S100n = var.S100n; S50n = var.S50n; S25n = var.S25n;
%         S = [S100n S50n S25n S25 S50 S100];
%         sft = S100|S50|S25|S100n|S50n|S25n;
        subj = var.Subj;
        res_all = abs(var.pitch - mufit);
        ind_A = find(~T0&~T400);
        ind_B0 = find(T0&~F0); ind_B400 = find(T400&~F400);
        ind_CA0 = find(F0); ind_CA400 = find(F400); % for fit without CB
%         ind_CA0 = find(F0&~sft); ind_CA400 = find(F400&~sft);
    
        taskA_fit(t) = mean(mufit(ind_A), 'omitnan');
        taskA_res(t) = mean(res_all(ind_A), 'omitnan');
        taskB_fit(t, 1) = mean(mufit(ind_B0), 'omitnan');
        taskB_res(t, 1) = mean(res_all(ind_B0), 'omitnan');
        taskB_fit(t, 2) = mean(mufit(ind_B400), 'omitnan')+400;
        taskB_res(t, 2) = mean(res_all(ind_B400), 'omitnan');
        taskCA_fit(t, 1) = mean(mufit(ind_CA0), 'omitnan');
        taskCA_res(t, 1) = mean(res_all(ind_CA0), 'omitnan');
        taskCA_fit(t, 2) = mean(mufit(ind_CA400), 'omitnan')+400;
        taskCA_res(t, 2) = mean(res_all(ind_CA400), 'omitnan');
%         for shift = 1:6
%             ind0 = find(T0&S(:, shift));
%             taskCB_fit(t, 1, shift) = mean(mufit(ind0), 'omitnan');
%             taskCB_res(t, 1, shift) = mean(res_all(ind0), 'omitnan');
%             ind400 = find(T400&S(:, shift));
%             taskCB_fit(t, 2, shift) = mean(mufit(ind400), 'omitnan')+400;
%             taskCB_res(t, 2, shift) = mean(res_all(ind400), 'omitnan');
%         end
    end
end
%%
figure;

sp = subplot(231);
hold on
plot(tt, taskA_fit, Color=col_A, LineWidth=2)
% tvalid = find(~isnan(mean_A));
% patch([tt(tvalid); flipud(tt(tvalid))], ...
%     [mean_A(tvalid)-sd_A(tvalid);  flipud(mean_A(tvalid)+sd_A(tvalid))], col_A, ...
%     'FaceAlpha',0.2, 'EdgeColor','none')
% plot(tt, betas_A, Color=col_A, LineWidth=2, LineStyle="--")
yline(0, LineStyle="--");
xline(0, LineStyle="--");
yl = ylim;
yl = [min(-100, yl(1)), max(100, yl(2))];
ylim(sp, yl);
xlim(sp, [-0.4, 1.2])
title("Task A Model Mean Fit", 'FontSize',15)
xlabel('Time (s)', 'FontSize',15)
ylabel('Magnitude (cent)', 'FontSize',15)
clear ylim

sp = subplot(234);
hold on
plot(tt, taskA_res, Color=col_A, LineWidth=2)
yline(0, LineStyle="--");
xline(0, LineStyle="--");
% yl = ylim;
% yl = [min(-100, yl(1)), max(100, yl(2))];
% ylim(sp, yl);
xlim(sp, [-0.4, 1.2])
title("Task A Model Mean Residual", 'FontSize',15)
xlabel('Time (s)', 'FontSize',15)
ylabel('Magnitude (cent)', 'FontSize',15)
clear ylim

sp = subplot(232);
hold on
plot(tt, taskB_fit(:, 1), Color=col_B(1, :), LineWidth=2)
plot(tt, taskB_fit(:, 2), Color=col_B(2, :), LineWidth=2)
yline(0, LineStyle="--");
yline(400, LineStyle="--");
xline(0, LineStyle="--");
yl = ylim;
yl = [min(-100, yl(1)), max(500, yl(2))];
ylim(sp, yl);
xlim(sp, [-0.4, 1.2])
title("Task B Model Mean Fit", 'FontSize',15)
legend("0 Cent","400 Cent", 'FontSize',12)
xlabel('Time (s)', 'FontSize',15)
ylabel('Magnitude (cent)', 'FontSize',15)
clear ylim

sp = subplot(235);
hold on
plot(tt, taskB_res(:, 1), Color=col_B(1, :), LineWidth=2)
plot(tt, taskB_res(:, 2), Color=col_B(2, :), LineWidth=2)
yl = ylim;
yl = [min(0, yl(1)), max(100, yl(2))];
ylim(sp, yl);
xlim(sp, [-0.4, 1.2])
xline(0, LineStyle="--");
title("Task B Model Mean Residual", 'FontSize',15)
legend("0 Cent","400 Cent", 'FontSize',12)
xlabel('Time (s)', 'FontSize',15)
ylabel('Magnitude (cent)', 'FontSize',15)
clear ylim

sp = subplot(233);
hold on
plot(tt, taskCA_fit(:, 1), Color=col_B(1, :), LineWidth=2)
plot(tt, taskCA_fit(:, 2), Color=col_B(2, :), LineWidth=2)
yline(0, LineStyle="--");
yline(400, LineStyle="--");
xline(0, LineStyle="--");
yl = ylim;
yl = [min(-100, yl(1)), max(500, yl(2))];
ylim(sp, yl);
xlim(sp, [-0.4, 1.2])
title("Task CA Model Mean Fit", 'FontSize',15)
legend("0 Cent","400 Cent", 'FontSize',12)
xlabel('Time (s)', 'FontSize',15)
ylabel('Magnitude (cent)', 'FontSize',15)
clear ylim

sp = subplot(236);
hold on
plot(tt, taskCA_res(:, 1), Color=col_B(1, :), LineWidth=2)
plot(tt, taskCA_res(:, 2), Color=col_B(2, :), LineWidth=2)
yl = ylim;
yl = [min(0, yl(1)), max(100, yl(2))];
ylim(sp, yl);
xlim(sp, [-0.4, 1.2])
xline(0, LineStyle="--");
title("Task CA Model Mean Residual", 'FontSize',15)
legend("0 Cent","400 Cent", 'FontSize',12)
xlabel('Time (s)', 'FontSize',15)
ylabel('Magnitude (cent)', 'FontSize',15)
clear ylim

% sp = subplot(244);
% hold on
% for targ = 1:2
%     for sh = 1:6
%         plot(tt, taskCB_fit(:, targ, sh), Color=col_CB(sh, :), LineWidth=2)
%     end
% end
% yline(0, LineStyle="--");
% yline(400, LineStyle="--");
% yl = ylim;
% yl = [min(-200, yl(1)), max(600, yl(2))];
% ylim(sp, yl);
% xlim(sp, [0, 1.2])
% title("Task CB Model Mean Fit", 'FontSize',15)
% xlabel('Time (s)', 'FontSize',15)
% ylabel('Magnitude (cent)', 'FontSize',15)
% legend("-100", "-50", "-25","25","50", "100",'FontSize',12,'Location', 'west');
% 
% sp = subplot(248);
% hold on
% for targ = 1:2
%     for sh = 1:6
%         if targ==1
%             plot(tt, taskCB_res(:, targ, sh), Color=col_CB(sh, :), LineWidth=2)
%         else
%             plot(tt, taskCB_res(:, targ, sh), Color=col_CB(sh, :), LineWidth=2, LineStyle="--")
%         end
%     end
% end
% yl = ylim;
% yl = [min(0, yl(1)), max(100, yl(2))];
% ylim(sp, yl);
% xlim(sp, [0, 1.2])
% title("Task CB Model Mean Residual", 'FontSize',15)
% xlabel('Time (s)', 'FontSize',15)
% ylabel('Magnitude (cent)', 'FontSize',15)
% legend("-100", "-50", "-25","25","50", "100",'FontSize',12,'Location', 'west');