% plot hypothesis tests for individual subjects using the GLMM fitted
% across all subjects + random intercept for each
clear
savepath = 'D:\Pitch adaptation\test\';
GLMMdata = savepath + "MaxiGLMMNoC_RandIntSubj_ALLSUB.mat";
Sigmoidfit = 'taskC_sigmoid_fit.mat';
file_CENT = savepath + "PARAMETERS_CENT.mat";
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
fmt = '%.3e';
% for each participant, plot performance in task CA against performance in
% task B
%%
% definition #1 of performance: averaged upper-lower?
T = length(glms);
tt = ([1:T]-(100-offset))/100 - 1/100; tt = tt(:);
% knee = 0.4;
% Nknee = knee*100;
% nsubj = 42;
Nparam = 3;
%nsub = length(stats_all_sub);
% taskB_ini = zeros(2, nsub); taskB_stab = zeros(2, nsub);
% taskCA_ini = zeros(2, nsub); taskCA_stab = zeros(2, nsub);
coefs_all = nan(Nparam, T);
rand_interc = nan(nsubj, T);

for t = 1:T
    % fisrt get fixed and random effect coefficients
    glm = glms(t).glme;
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
indT = (100-offset)+(1:T);
pA = nan(nsubj, 1); pB = nan(nsubj, 2);
% pCA = nan(nsubj, 2); pCB = nan(nsubj, 2, 6);
for s = 1:5
    data = PARAMETERS_CENT(subj_idx(s)); 
    sname = data.Sbjnames;
    if ~isempty(sname)
%         sv = sv+1; %valid subj
        % first get all betas into a matrix       
        itc = rand_interc(:, subj_idx(s));
                % get all mean pitch traces
        mean_A = data.A_mean(indT); mean_B = data.B_mean(indT, :); 
%         mean_CA = data.CA_mean(indT, :); 
%         mean_CB = squeeze(data.CB_mean(indT, :, :, :)); % including 0 shift
%         mean_CB = mean_CB(:, :, [1:3, 5:7]); % remove 0 shift
        sd_A = sqrt(data.A_variance(indT));
        sd_B = sqrt(data.B_variance(indT, :));
%         sd_CA = sqrt(data.CA_variance(indT, :));
%         var_CB = squeeze(data.CB_variance(indT, :, :, :)); % including 0 shift
%         sd_CB = sqrt(var_CB(:, :, [1:3, 5:7])); 
        % reconstruct pitch traces from betas and KS test
        betas_B = zeros(T, 2); %pB = zeros(1,2);
%         betas_CA = zeros(T, 2); %pCA = zeros(1,2);
%         betas_CB = zeros(T, 2, 6); %pCB = zeros(2,6);

        betas_A = betas(:, 1) + itc; 
        [~, pA(s)] = kstest2(betas_A, mean_A);
        for targ = 1:2
            betas_B(:, targ) = betas(:, 1) + betas(:, targ+1) + itc;
            if targ == 2
                betas_B(:, targ) = betas_B(:, targ) + 400;
            end
            [~, pB(s, targ)] = kstest2(betas_B(:, targ), mean_B(:, targ));
%             betas_CA(:, targ) = betas(:, 1) + betas(:, targ+1) + betas(:, targ+3) + itc;
%             if targ == 2
%                 betas_CA(:, targ) = betas_CA(:, targ) + 400;
%             end
%             [~, pCA(s, targ)] = kstest2(betas_CA(:, targ), mean_CA(:, targ));
%             for shft = 1:6
%                 betas_CB(:, targ, shft) = betas(:, 1) + betas(:, targ+1) + ...
%                     betas(:, targ+3) + betas(:, shft+5) + itc;
%                 if targ == 2
%                     betas_CB(:, targ, shft) = betas_CB(:, targ, shft) + 400;
%                 end
%                 [~, pCB(s, targ, shft)] = kstest2(betas_CB(:, targ, shft), mean_CB(:, targ, shft));
%             end
        end
        % plot individual figure
        fig = figure();
        sp = subplot(131);
        hold on
        plot(tt, mean_A, Color=col_A, LineWidth=2)
        tvalid = find(~isnan(mean_A));
        patch([tt(tvalid); flipud(tt(tvalid))], ...
            [mean_A(tvalid)-sd_A(tvalid);  flipud(mean_A(tvalid)+sd_A(tvalid))], col_A, ...
            'FaceAlpha',0.2, 'EdgeColor','none')
        plot(tt, betas_A, Color=col_A, LineWidth=2, LineStyle="--")
        yline(0, LineStyle="--");
        xline(0, LineStyle="--");
        yl = ylim;
        yl = [min(-100, yl(1)), max(100, yl(2))];
        ylim(sp, yl);
        xlim(sp, [-.4, 1.2])
        legend("Mean","SD", "Fit", 'FontSize',10 );
        title("Task A", 'FontSize',15)
        xlabel('Time (s)', 'FontSize',15)
        ylabel('Magnitude (cent)', 'FontSize',15)
        clear ylim

%         sp=subplot(245);
%         hold on
%         plot(tt, betas_A, Color=col_A, LineWidth=2)
%         yline(0, LineStyle="--");
%         yl = ylim;
%         yl = [min(-100, yl(1)), max(100, yl(2))];
%         ylim(sp, yl);
%         title("Reconstructed Task A", 'FontSize',15)
%         xlabel('Time (s)', 'FontSize',15)
%         ylabel('Magnitude (cent)', 'FontSize',15)
%         legend(['p=' num2str(pA, fmt)], 'FontSize',10)
%         clear ylim

        % subplot 2: target
        sp=subplot(132);
        hold on
        for targ = 1:2
            plot(tt, mean_B(:, targ), Color=col_B(targ, :), LineWidth=2)
            tvalid = find(~isnan(mean_B(:, targ)));
            patch([tt(tvalid); flipud(tt(tvalid))], [mean_B(tvalid, targ)-sd_B(tvalid, targ); ...
            flipud(mean_B(tvalid, targ)+sd_B(tvalid, targ))], ...
            col_B(targ, :),'FaceAlpha',0.2, 'EdgeColor','none')
            plot(tt, betas_B(:, targ), Color=col_B(targ, :), LineWidth=2, LineStyle="--")
        end
        yline(0, LineStyle="--");
        yline(400, LineStyle="--");
        yl = ylim;
        yl = [min(-200, yl(1)), max(600, yl(2))];
        ylim(sp, yl);
        xlim(sp, [-.4, 1.2])
        xline(0, LineStyle="--");
        title("Task B", 'FontSize',15)
        xlabel('Time (s)', 'FontSize',15)
        ylabel('Magnitude (cent)', 'FontSize',15)
        legend("Mean 0 Cent","SD 0 Cent","Fit 0 Cent", "Mean 400 Cent","SD 400 Cent","Fit 400 Cent",...
            'FontSize',10,'NumColumns',2)
        hold off
        clear ylim

%         sp=subplot(246);
%         hold on
%         for targ = 1:2
%             plot(tt, betas_B(:, targ), Color=col_B(targ, :), LineWidth=2)
%         end
%         yline(0, LineStyle="--");
%         yline(400, LineStyle="--");
%         yl = ylim;
%         yl = [min(-200, yl(1)), max(600, yl(2))];
%         ylim(sp, yl);
%         title("Reconstructed Task B", 'FontSize',15)
%         xlabel('Time (s)', 'FontSize',15)
%         ylabel('Magnitude (cent)', 'FontSize',15)
%         legend(['p=' num2str(pB(1), fmt)], ['p=' num2str(pB(2), fmt)],'FontSize',10,'NumColumns',2)
%         clear ylim

%         % subplot 3: feedback
%         sp=subplot(133);
%         hold on
%         for targ = 1:2
%             plot(tt, mean_CA(:, targ), Color=col_B(targ, :), LineWidth=2)
%             tvalid = find(~isnan(mean_CA(:, targ)));
%             patch([tt(tvalid); flipud(tt(tvalid))], [mean_CA(tvalid, targ)-sd_CA(tvalid, targ); ...
%             flipud(mean_CA(tvalid, targ)+sd_CA(tvalid, targ))], ...
%             col_B(targ, :),'FaceAlpha',0.2, 'EdgeColor','none')
%             plot(tt, betas_CA(:, targ), Color=col_B(targ, :), LineWidth=2, LineStyle="--")
%         end
%         yline(0, LineStyle="--");
%         yline(400, LineStyle="--");
%         xline(0, LineStyle="--");
%         yl = ylim;
%         yl = [min(-200, yl(1)), max(600, yl(2))];
%         ylim(sp, yl);
%         xlim(sp, [-.4, 1.2])
%         title("Task CA", 'FontSize',15)
%         xlabel('Time (s)', 'FontSize',15)
%         ylabel('Magnitude (cent)', 'FontSize',15)
%         legend("Mean 0 Cent","SD 0 Cent","Fit 0 Cent", "Mean 400 Cent","SD 400 Cent","Fit 400 Cent",...
%             'FontSize',10,'NumColumns',2)
%         clear ylim

%         sp=subplot(247);
%         hold on
%         for targ = 1:2
%             plot(tt, betas_CA(:, targ), Color=col_B(targ, :), LineWidth=2)
%         end
%         yline(0, LineStyle="--");
%         yline(400, LineStyle="--");
%         yl = ylim;
%         yl = [min(-200, yl(1)), max(600, yl(2))];
%         ylim(sp, yl);
%         title("Reconstructed Task CA", 'FontSize',15)
%         xlabel('Time (s)', 'FontSize',15)
%         ylabel('Magnitude (cent)', 'FontSize',15)
%         legend(['p=' num2str(pCA(1), fmt)], ['p=' num2str(pCA(2), fmt)],'FontSize',10,'NumColumns',2)
%         clear ylim

%         sp=subplot(144);
%         hold on
%         for targ = 1:2
%             for sh = 1:6
%                 plot(tt, mean_CB(:, targ, sh), Color=col_CB(sh, :), LineWidth=2)
%                 tvalid = find(~isnan(mean_CB(:, targ, sh)));
%                 patch([tt(tvalid); flipud(tt(tvalid))], [mean_CB(tvalid, targ, sh)-sd_CB(tvalid, targ, sh); ...
%                 flipud(mean_CB(tvalid, targ, sh)+sd_CB(tvalid, targ, sh))], ...
%                 col_CB(sh, :),'FaceAlpha',0.2, 'EdgeColor','none')
%                 plot(tt, betas_CB(:, targ, sh), Color=col_CB(sh, :), LineWidth=2, LineStyle="--")
%             end
%         end
%         yline(0, LineStyle="--");
%         yline(400, LineStyle="--");
%         yl = ylim;
%         yl = [min(-200, yl(1)), max(600, yl(2))];
%         ylim(sp, yl);
%         xlim(sp, [0, 1.2])
%         title("Task CB", 'FontSize',15)
%         xlabel('Time (s)', 'FontSize',15)
%         ylabel('Magnitude (cent)', 'FontSize',15)
%         legend("-100","", "Fit", "-50", "", "Fit","-25","", "Fit","25","", "Fit","50","", ...
%             "Fit","100","", "Fit",'FontSize',10,'Location', 'west');
%         clear ylim

%         sp=subplot(248);
%         hold on
%         for targ = 1:2
%             for sh = 1:6
%                 plot(tt, betas_CB(:, targ, sh), Color=col_CB(sh, :), LineWidth=2)
%             end
%         end
%         yline(0, LineStyle="--");
%         yline(400, LineStyle="--");
%         yl = ylim;
%         yl = [min(-200, yl(1)), max(600, yl(2))];
%         ylim = yl;
%         title("Reconstructed Task CB", 'FontSize',15)
%         xlabel('Time (s)', 'FontSize',15)
%         ylabel('Magnitude (cent)', 'FontSize',15)
%         legend(['p=' num2str(pCB(1,1), fmt), ', ' num2str(pCB(2,1), fmt)], ...
%             ['p=' num2str(pCB(1,2), fmt), ', ' num2str(pCB(2,2), fmt)], ...
%             ['p=' num2str(pCB(1,3), fmt), ', ' num2str(pCB(2,3), fmt)], ...
%             ['p=' num2str(pCB(1,4), fmt), ', ' num2str(pCB(2,4), fmt)], ...
%             ['p=' num2str(pCB(1,5), fmt), ', ' num2str(pCB(2,5), fmt)], ...
%             ['p=' num2str(pCB(1,6), fmt), ', ' num2str(pCB(2,6), fmt)], ...
%             'Location', 'west','FontSize',10)
%         clear ylim
%         legend("-100", "-50", "-25", "25","50","100",'FontSize',10);

        sgtitle(data.Sbjnames)
        fig.WindowState = 'maximized';
        exportgraphics(fig, folder + data.Sbjnames + "_GLMMvsIndDataNoCB" + FORMAT, "Resolution",RESOLUTION,'ContentType',CONTENT_TYPE)
%         close(fig);
    end
end

%% plot parameters over time
% effect of target vs. feedback
h = figure();
tiledlayout(1,2);
nexttile
% subplot(121)
scatter(betas(31:end, 2), betas(31:end, 4), [], tt(31:end), 'fill')
xlabel('Effect of target', 'FontSize', 12)
ylabel('Effect of feedback', 'FontSize', 12)
title('Target 0 Cent', 'FontSize', 12)

% subplot(122)
nexttile
scatter(betas(31:end, 3), betas(31:end, 5), [], tt(31:end), 'fill')
xlabel('Effect of target', 'FontSize', 12)
ylabel('Effect of feedback', 'FontSize', 12)
title('Target 400 Cent', 'FontSize', 12)
c = colorbar;
c.Ticks = -0.2:0.2:1.2;
c.Label.String = 'Time (s)';
c.FontSize = 12;


%%
col_B = [[0 0.4470 0.7410];[0.4940 0.1840 0.5560];]; 
figure
subplot(121)
hold on
plot(taskB_ini(1, :), taskCA_ini(1, :), "color", col_B(1, :), "Marker", 'o', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
plot(taskB_ini(2, :), taskCA_ini(2, :), "color", col_B(2, :), "Marker", '*', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
plot(0:50, 0:50, '--k')
legend("0 Cent", "400 Cent", 'FontSize',12)
xlabel("Task B Performance", 'FontSize',14)
ylabel("Task C Performance", 'FontSize',14)
title("Initialization 0-0.4s", 'FontSize',14)

subplot(122)
hold on
plot(taskB_stab(1, :), taskCA_stab(1, :), "color", col_B(1, :), "Marker", 'o', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
plot(taskB_stab(2, :), taskCA_stab(2, :), "color", col_B(2, :), "Marker", '*', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1);
plot(0:25, 0:25, '--k')
legend("0 Cent", "400 Cent", 'FontSize',12)
xlabel("Task B Performance", 'FontSize',14)
ylabel("Task C Performance", 'FontSize',14)
title("Stablization 0.4-1.2s", 'FontSize',14)
