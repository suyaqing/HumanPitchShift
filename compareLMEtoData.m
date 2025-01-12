% test whether individual betas can reconstruct the mean performance of
% each task. fitted and actual data are compared using two-sample KS test
tic
savepath = 'D:\Pitch adaptation\test\';
file_CENT = savepath + "PARAMETERS_CENT.mat";
file_var = savepath + "variables.mat";
file_GLMM = savepath + "IndividualMaxiGLMM.mat";
load(file_CENT)
load(file_var)
load(file_GLMM)
EXCLUDE = {'f12NC', 'f18NC','f19NC','f25NC','m11NC','m12NC','m18NC','m23','m27'};

%%% Set Resolution and Format for Plots
folder = savepath + "\Plots\GLMM\";
RESOLUTION  = 300; %Resolution in dots per inch (DPI), specified as a whole number that is greater than or equal to 1.
FORMAT = ".pdf"; %Possible formats are: '.jpg', '.png', '.tif', '.gif', '.pdf', '.emf' or '.eps'
CONTENT_TYPE = 'auto'; 

T = 120;
tt = (1:T)/100 - 1/100; tt = tt(:);
col_A = [0 0.4470 0.7410];
col_B = [[0 0.4470 0.7410];[0.4940 0.1840 0.5560];]; 
col_CB = [[0.4940 0.1840 0.5560];[0 0.4470 0.7410];[0.3010 0.7450 0.9330];[0.4660 0.6740 0.1880];[0.9290 0.6940 0.1250];[0.8500 0.3250 0.0980];[0.6350 0.0780 0.1840];]; %RGB color code, nice for 7 shifts.
trend_mean = zeros(T, 1); trend_SE = trend_mean; trend_p = trend_mean;
target_mean = zeros(T, 2); target_SE = target_mean; target_p = target_mean;
FB_mean = zeros(T, 2); FB_SE = FB_mean; FB_p = FB_mean;
shift_mean = zeros(T, 6); shift_SE = shift_mean; shift_p = shift_mean;
% construct input table for the GLM, one subject after another
varNames = {'pitch', 'T0', 'T400', 'F0', 'F400', 'S100n', 'S50n', 'S25n', ...
    'S25', 'S50', 'S100'};
nvar = length(varNames);
tA = 16;
tB = 16;
tCA = 70;
tCB = 10;
fmt = '%.3e';
%%
tic
sv = 0;
for s = 1:nsubj
    data = PARAMETERS_CENT(s); 
    subjname = stats_all_sub(s).sbjName;
    if ~isempty(subjname)
%         sv = sv+1; %valid subj
        % first get all betas into a matrix
        betas = zeros(T, 11);
        glmms = stats_all_sub(s).glmm;
        for t = 1:T
            str = dataset2struct(glmms(t).coeff);
            for ind = 1:11
                betas(t, ind) = str(ind).Estimate;
            end
        end
        % get all mean pitch traces
        mean_A = data.A_mean(1:T); mean_B = data.B_mean(1:T, :); 
        mean_CA = data.CA_mean(1:T, :); 
        mean_CB = squeeze(data.CB_mean(1:T, :, :, :)); % including 0 shift
        mean_CB = mean_CB(:, :, [1:3, 5:7]); % remove 0 shift
        sd_A = sqrt(data.A_variance(1:T));
        sd_B = sqrt(data.B_variance(1:T, :));
        sd_CA = sqrt(data.CA_variance(1:T, :));
        var_CB = squeeze(data.CB_variance(1:T, :, :, :)); % including 0 shift
        sd_CB = sqrt(var_CB(:, :, [1:3, 5:7])); 
        % reconstruct pitch traces from betas and KS test
        betas_B = zeros(T, 2); pB = zeros(1,2);
        betas_CA = zeros(T, 2); pCA = zeros(1,2);
        betas_CB = zeros(T, 2, 6); pCB = zeros(2,6);

        betas_A = betas(:, 1); [~, pA] = kstest2(betas_A, mean_A);
        for targ = 1:2
            betas_B(:, targ) = betas(:, 1) + betas(:, targ+1);
            if targ == 2
                betas_B(:, targ) = betas_B(:, targ) + 400;
            end
            [~, pB(targ)] = kstest2(betas_B(:, targ), mean_B(:, targ));
            betas_CA(:, targ) = betas(:, 1) + betas(:, targ+1) + betas(:, targ+3);
            if targ == 2
                betas_CA(:, targ) = betas_CA(:, targ) + 400;
            end
            [~, pCA(targ)] = kstest2(betas_CA(:, targ), mean_CA(:, targ));
            for shft = 1:6
                betas_CB(:, targ, shft) = betas(:, 1) + betas(:, targ+1) + ...
                    betas(:, targ+3) + betas(:, shft+5);
                if targ == 2
                    betas_CB(:, targ, shft) = betas_CB(:, targ, shft) + 400;
                end
                [~, pCB(targ, shft)] = kstest2(betas_CB(:, targ, shft), mean_CB(:, targ, shft));
            end
        end
        % plot individual figure
        fig = figure();
        sp = subplot(241);
        hold on
        plot(tt, mean_A, Color=col_A, LineWidth=2)
        patch([tt; flipud(tt)], [mean_A-sd_A;  flipud(mean_A+sd_A)], col_A,'FaceAlpha',0.2, 'EdgeColor','none')
        yline(0, LineStyle="--");
        yl = ylim;
        yl = [min(-100, yl(1)), max(100, yl(2))];
        ylim(sp, yl);
        legend("Mean","SD",'FontSize',10 );
        title("Task A", 'FontSize',15)
        xlabel('Time (s)', 'FontSize',15)
        ylabel('Magnitude (cent)', 'FontSize',15)
        clear ylim

        sp=subplot(245);
        hold on
        plot(tt, betas_A, Color=col_A, LineWidth=2)
        yline(0, LineStyle="--");
        yl = ylim;
        yl = [min(-100, yl(1)), max(100, yl(2))];
        ylim(sp, yl);
        title("Reconstructed Task A", 'FontSize',15)
        xlabel('Time (s)', 'FontSize',15)
        ylabel('Magnitude (cent)', 'FontSize',15)
        legend(['p=' num2str(pA, fmt)], 'FontSize',10)
        clear ylim

        % subplot 2: target
        sp=subplot(242);
        hold on
        for targ = 1:2
            plot(tt, mean_B(:, targ), Color=col_B(targ, :), LineWidth=2)
            patch([tt; flipud(tt)], [mean_B(:, targ)-sd_B(:, targ); ...
            flipud(mean_B(:, targ)+sd_B(:, targ))], ...
            col_B(targ, :),'FaceAlpha',0.2, 'EdgeColor','none')
        end
        yline(0, LineStyle="--");
        yline(400, LineStyle="--");
        yl = ylim;
        yl = [min(-200, yl(1)), max(600, yl(2))];
        ylim(sp, yl);
        title("Task B", 'FontSize',15)
        xlabel('Time (s)', 'FontSize',15)
        ylabel('Magnitude (cent)', 'FontSize',15)
        legend("Mean 0 Cent","SD 0 Cent", "Mean 400 Cent","SD 400 Cent",'FontSize',10,'NumColumns',2)
        hold off
        clear ylim

        sp=subplot(246);
        hold on
        for targ = 1:2
            plot(tt, betas_B(:, targ), Color=col_B(targ, :), LineWidth=2)
        end
        yline(0, LineStyle="--");
        yline(400, LineStyle="--");
        yl = ylim;
        yl = [min(-200, yl(1)), max(600, yl(2))];
        ylim(sp, yl);
        title("Reconstructed Task B", 'FontSize',15)
        xlabel('Time (s)', 'FontSize',15)
        ylabel('Magnitude (cent)', 'FontSize',15)
        legend(['p=' num2str(pB(1), fmt)], ['p=' num2str(pB(2), fmt)],'FontSize',10,'NumColumns',2)
        clear ylim

        % subplot 3: feedback
        sp=subplot(243);
        hold on
        for targ = 1:2
            plot(tt, mean_CA(:, targ), Color=col_B(targ, :), LineWidth=2)
            patch([tt; flipud(tt)], [mean_CA(:, targ)-sd_CA(:, targ); ...
            flipud(mean_CA(:, targ)+sd_CA(:, targ))], ...
            col_B(targ, :),'FaceAlpha',0.2, 'EdgeColor','none')
        end
        yline(0, LineStyle="--");
        yline(400, LineStyle="--");
        yl = ylim;
        yl = [min(-200, yl(1)), max(600, yl(2))];
        ylim(sp, yl);
        title("Task CA", 'FontSize',15)
        xlabel('Time (s)', 'FontSize',15)
        ylabel('Magnitude (cent)', 'FontSize',15)
        legend("Mean 0 Cent","SD 0 Cent", "Mean 400 Cent","SD 400 Cent",'FontSize',10,'NumColumns',2)
        clear ylim

        sp=subplot(247);
        hold on
        for targ = 1:2
            plot(tt, betas_CA(:, targ), Color=col_B(targ, :), LineWidth=2)
        end
        yline(0, LineStyle="--");
        yline(400, LineStyle="--");
        yl = ylim;
        yl = [min(-200, yl(1)), max(600, yl(2))];
        ylim(sp, yl);
        title("Reconstructed Task CA", 'FontSize',15)
        xlabel('Time (s)', 'FontSize',15)
        ylabel('Magnitude (cent)', 'FontSize',15)
        legend(['p=' num2str(pCA(1), fmt)], ['p=' num2str(pCA(2), fmt)],'FontSize',10,'NumColumns',2)
        clear ylim

        sp=subplot(244);
        hold on
        for targ = 1:2
            for sh = 1:6
                plot(tt, mean_CB(:, targ, sh), Color=col_CB(sh, :), LineWidth=2)
                patch([tt; flipud(tt)], [mean_CB(:, targ, sh)-sd_CB(:, targ, sh); ...
                flipud(mean_CB(:, targ, sh)+sd_CB(:, targ, sh))], ...
                col_CB(sh, :),'FaceAlpha',0.2, 'EdgeColor','none')
            end
        end
        yline(0, LineStyle="--");
        yline(400, LineStyle="--");
        yl = ylim;
        yl = [min(-200, yl(1)), max(600, yl(2))];
        ylim(sp, yl);
        title("Task CB", 'FontSize',15)
        xlabel('Time (s)', 'FontSize',15)
        ylabel('Magnitude (cent)', 'FontSize',15)
        legend("-100","", "-50", "", "-25","", "25","", "50","", "100","", 'FontSize',10,'Location', 'west');
        clear ylim

        sp=subplot(248);
        hold on
        for targ = 1:2
            for sh = 1:6
                plot(tt, betas_CB(:, targ, sh), Color=col_CB(sh, :), LineWidth=2)
            end
        end
        yline(0, LineStyle="--");
        yline(400, LineStyle="--");
        yl = ylim;
        yl = [min(-200, yl(1)), max(600, yl(2))];
        ylim = yl;
        title("Reconstructed Task CB", 'FontSize',15)
        xlabel('Time (s)', 'FontSize',15)
        ylabel('Magnitude (cent)', 'FontSize',15)
        legend(['p=' num2str(pCB(1,1), fmt), ', ' num2str(pCB(2,1), fmt)], ...
            ['p=' num2str(pCB(1,2), fmt), ', ' num2str(pCB(2,2), fmt)], ...
            ['p=' num2str(pCB(1,3), fmt), ', ' num2str(pCB(2,3), fmt)], ...
            ['p=' num2str(pCB(1,4), fmt), ', ' num2str(pCB(2,4), fmt)], ...
            ['p=' num2str(pCB(1,5), fmt), ', ' num2str(pCB(2,5), fmt)], ...
            ['p=' num2str(pCB(1,6), fmt), ', ' num2str(pCB(2,6), fmt)], ...
            'Location', 'west','FontSize',10)
        clear ylim
%         legend("-100", "-50", "-25", "25","50","100",'FontSize',10);

        sgtitle(data.Sbjnames)
        fig.WindowState = 'maximized';
        exportgraphics(fig, folder + data.Sbjnames + "_GLMMvsData" + FORMAT, "Resolution",RESOLUTION,'ContentType',CONTENT_TYPE)
        close(fig);
    end
end

toc
