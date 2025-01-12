
%%% Creates one plot per participant with showing the results of test A, B, CA, CB and D and saves them in the Folder "\Plots\ALL\" created at the defined "savepath" directory.
%%% NEEDS PARAMETERS: "FORMAT", "savepath", "TIME", "TIME_D", "RESOLUTION", 
%%% NEEDS FILE: PARAMETERS_CENTS.


% Author: Philipp Eugster <eugsteph@student.ethz.ch>
% Created: 10.12.2022




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Set Colors for different tests.

col_A = [0 0.4470 0.7410];
col_A_raw = [0 0.4470 0.7410 0.3];
col_B = [[0 0.4470 0.7410];[0.4940 0.1840 0.5560];];
col_B_raw = [[0 0.4470 0.7410 0.3];[0.4940 0.1840 0.5560 0.3];];
col_CA = col_B;
col_CA_raw =  [[0 0.4470 0.7410 0.1];[0.4940 0.1840 0.5560 0.1];];
col_CB = [[0.4940 0.1840 0.5560];[0 0.4470 0.7410];[0.3010 0.7450 0.9330];[0.4660 0.6740 0.1880];[0.9290 0.6940 0.1250];[0.8500 0.3250 0.0980];[0.6350 0.0780 0.1840];]; %RGB color code, nice for 7 shifts.
col_CB_raw = [[0.4940 0.1840 0.5560 0.2];[0 0.4470 0.7410 0.2];[0.3010 0.7450 0.9330 0.2];[0.4660 0.6740 0.1880 0.2];[0.9290 0.6940 0.1250 0.2];[0.8500 0.3250 0.0980 0.2];[0.6350 0.0780 0.1840 0.2];]; %RGB color code, nice for 7 shifts.
col_D = col_CB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Generate directory to save plots ("savepath" is defined in Main.m)

folder = savepath + "\Plots\ALL\";
if ~exist(folder, 'dir')
    mkdir(folder);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Create plots for the different tests from the PARAMETERS_CENT.m file.
%%
for i = 1:nsubj

    %%%%%%%%%%%%%%%%%%%%%%%% Plot Test A in subplot %%%%%%%%%%%%%%%%%%%%%%%
    fig = figure();
    subplot(3,3,1)
    keepIndex = ~isnan(PARAMETERS_CENT(i).A_mean);
    x = TIME(keepIndex);                                                       
    y = PARAMETERS_CENT(i).A_mean(keepIndex);                                 
    sd = sqrt(PARAMETERS_CENT(i).A_variance(keepIndex));
    plot(x(:), y, Color=col_A,LineWidth=2)
    hold on
    patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_A,'FaceAlpha',0.2, 'EdgeColor','none')
    idx = find(~PARAMETERS_CENT(i).A_outlr == 1);
    plot(TIME,PARAMETERS_CENT(i).A_cent(:,idx(1)), Color = col_A_raw);
    plot(TIME,PARAMETERS_CENT(i).A_cent(:,~PARAMETERS_CENT(i).A_outlr), 'HandleVisibility', 'off',Color = col_A_raw);
    yline(0,'-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
    
    if PARAMETERS_CENT(i).A_nOutlr > 0
        plot(TIME, PARAMETERS_CENT(i).A_cent(:,PARAMETERS_CENT(i).A_outlr), Color = col_A, LineStyle='--');
        legend("Mean","STDEV", "Rawdata", "Outlier",'FontSize',5 );
    else
        legend("Mean","STDEV", "Rawdata",'FontSize',5 );
    end

    title("Test A");
    xlabel('Time in seconds');
    ylabel('Cent');
    ylim([-400 400])
    xlim([0 2])
    hold off

    %%%%%%%%%%%%%%%%%%%%%%%% Plot Test B in subplot %%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,3,2)
    hold on
    keepIndex = ~isnan(PARAMETERS_CENT(i).B_mean);
    for id = 1:length(keepIndex(1,:))
        x = TIME(keepIndex(:,id));
        y = PARAMETERS_CENT(i).B_mean(keepIndex(:,id),id);
        sd = sqrt(PARAMETERS_CENT(i).B_variance(keepIndex(:,id),id));
        plot(x(:), y, Color=col_B(id,:),LineWidth=2)
        patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_B(id,:),'FaceAlpha',0.2, 'EdgeColor','none')
        idx = find(~PARAMETERS_CENT(i).B_outlr(:,:,id) == 1);
        plot(TIME,PARAMETERS_CENT(i).B_cent(:,idx(1),id), Color = col_B_raw(id,:));
        plot(TIME,PARAMETERS_CENT(i).B_cent(:,~PARAMETERS_CENT(i).B_outlr(:,:,id),id), 'HandleVisibility', 'off',Color = col_B_raw(id,:));
        if PARAMETERS_CENT(i).B_nOutlr(id) == 1
            plot(TIME, PARAMETERS_CENT(i).B_cent(:,PARAMETERS_CENT(i).B_outlr(:,:,id) ,id), LineStyle="--",Color = col_B(id,:) );
        elseif PARAMETERS_CENT(i).B_nOutlr(id) > 1
            idx = find(PARAMETERS_CENT(i).B_outlr(:,:,id) == 1);
            plot(TIME, PARAMETERS_CENT(i).B_cent(:,idx(1),id), LineStyle="--",Color = col_B(id,:) ); 
            plot(TIME, PARAMETERS_CENT(i).B_cent(:,PARAMETERS_CENT(i).B_outlr(:,:,id),id),'HandleVisibility', 'off', LineStyle="--",Color = col_B(id,:) );
        end
    end


    yline(0,'-.b','HandleVisibility','off',LineWidth=0.8, Color="k");

    if sum(PARAMETERS_CENT(i).B_nOutlr(:,:,1)) > 0 &&  sum(PARAMETERS_CENT(i).B_nOutlr(:,:,2)) > 0
        legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    elseif sum(PARAMETERS_CENT(i).B_nOutlr(:,:,1)) > 0
        legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    elseif sum(PARAMETERS_CENT(i).B_nOutlr(:,:,2)) > 0
         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    else
         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    end

    title("Test B");
    yline(0, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
    yline(400, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
    ylim([-400 800])

    xlabel('Time in seconds');
    ylabel('Cent');
    hold off
   %%%%%%%%%%%%%%%%%%%%%%%% Plot effect of sensory target in subplot %%%%%%%%%%%%%%%%%%%%%%

    subplot(3,3,4)
    hold on

    keepIndex = ~isnan(PARAMETERS_CENT(i).CA_mean);
    for id = 1:length(keepIndex(1,:))
        x = TIME(keepIndex(:,id));
        y = PARAMETERS_CENT(i).B_sensoryTargetEffect(keepIndex(:,id),id);
%         sd = sqrt(PARAMETERS_CENT(i).CA_variance(keepIndex(:,id),id));
        plot(x(:), y, Color=col_B(id,:),LineWidth=2)
%         patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_CA(id,:),'FaceAlpha',0.2, 'EdgeColor','none')
%         idx = find(~PARAMETERS_CENT(i).CA_outlr(:,:,id) == 1);
%         plot(TIME,PARAMETERS_CENT(i).CA_cent(:,idx(1),id), Color = col_CA_raw(id,:));
%         plot(TIME,PARAMETERS_CENT(i).CA_cent(:,~PARAMETERS_CENT(i).CA_outlr(:,:,id),id), 'HandleVisibility', 'off',Color = col_CA_raw(id,:));
%         if PARAMETERS_CENT(i).CA_nOutlr(id) == 1
%             plot(TIME, PARAMETERS_CENT(i).CA_cent(:,PARAMETERS_CENT(i).CA_outlr(:,:,id) ,id), LineStyle="--",Color = col_B(id,:) );
%         elseif PARAMETERS_CENT(i).CA_nOutlr(id) > 1
%             idx = find(PARAMETERS_CENT(i).CA_outlr(:,:,id) == 1);
%             plot(TIME, PARAMETERS_CENT(i).CA_cent(:,idx(1),id), LineStyle="--",Color = col_B(id,:) ); 
%             plot(TIME, PARAMETERS_CENT(i).CA_cent(:,PARAMETERS_CENT(i).CA_outlr(:,:,id),id),'HandleVisibility', 'off', LineStyle="--",Color = col_B(id,:) );
%         end
    end

    title("Test B - Sensory Target");
    yline(0, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
%     yline(400, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
    ylim([-400 800]);
    xlabel('Time in seconds');
    ylabel('Cent');
    legend("Effect 0Cent", "Effect 400Cent");
    
%     if sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,1)) > 0 &&  sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,2)) > 0
%         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
%     elseif sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,1)) > 0
%         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
%     elseif sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,2)) > 0
%          legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
%     else
%          legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
%     end
    hold off


    %%%%%%%%%%%%%%%%%%%%%%%% Plot Test CA in subplot %%%%%%%%%%%%%%%%%%%%%%

    subplot(3,3,3)
    hold on

    keepIndex = ~isnan(PARAMETERS_CENT(i).CA_mean);
    for id = 1:length(keepIndex(1,:))
        x = TIME(keepIndex(:,id));
        y = PARAMETERS_CENT(i).CA_mean(keepIndex(:,id),id);
        sd = sqrt(PARAMETERS_CENT(i).CA_variance(keepIndex(:,id),id));
        plot(x(:), y, Color=col_CA(id,:),LineWidth=2)
        patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_CA(id,:),'FaceAlpha',0.2, 'EdgeColor','none')
        idx = find(~PARAMETERS_CENT(i).CA_outlr(:,:,id) == 1);
        plot(TIME,PARAMETERS_CENT(i).CA_cent(:,idx(1),id), Color = col_CA_raw(id,:));
        plot(TIME,PARAMETERS_CENT(i).CA_cent(:,~PARAMETERS_CENT(i).CA_outlr(:,:,id),id), 'HandleVisibility', 'off',Color = col_CA_raw(id,:));
        if PARAMETERS_CENT(i).CA_nOutlr(id) == 1
            plot(TIME, PARAMETERS_CENT(i).CA_cent(:,PARAMETERS_CENT(i).CA_outlr(:,:,id) ,id), LineStyle="--",Color = col_B(id,:) );
        elseif PARAMETERS_CENT(i).CA_nOutlr(id) > 1
            idx = find(PARAMETERS_CENT(i).CA_outlr(:,:,id) == 1);
            plot(TIME, PARAMETERS_CENT(i).CA_cent(:,idx(1),id), LineStyle="--",Color = col_B(id,:) ); 
            plot(TIME, PARAMETERS_CENT(i).CA_cent(:,PARAMETERS_CENT(i).CA_outlr(:,:,id),id),'HandleVisibility', 'off', LineStyle="--",Color = col_B(id,:) );
        end
    end

    title("Test C - part A");
    yline(0, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
    yline(400, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
    ylim([-400 800]);
    xlabel('Time in seconds');
    ylabel('Cent');
    
    if sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,1)) > 0 &&  sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,2)) > 0
        legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    elseif sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,1)) > 0
        legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    elseif sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,2)) > 0
         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    else
         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    end
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%% Plot effect of vocal feedback in subplot %%%%%%%%%%%%%%%%%%%%%%

    subplot(3,3,5)
    hold on

    keepIndex = ~isnan(PARAMETERS_CENT(i).CA_mean);
    for id = 1:length(keepIndex(1,:))
        x = TIME(keepIndex(:,id));
        y = PARAMETERS_CENT(i).CA_vocalFeebackEffect(keepIndex(:,id),id);
%         sd = sqrt(PARAMETERS_CENT(i).CA_variance(keepIndex(:,id),id));
        plot(x(:), y, Color=col_CA(id,:),LineWidth=2)
%         patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_CA(id,:),'FaceAlpha',0.2, 'EdgeColor','none')
%         idx = find(~PARAMETERS_CENT(i).CA_outlr(:,:,id) == 1);
%         plot(TIME,PARAMETERS_CENT(i).CA_cent(:,idx(1),id), Color = col_CA_raw(id,:));
%         plot(TIME,PARAMETERS_CENT(i).CA_cent(:,~PARAMETERS_CENT(i).CA_outlr(:,:,id),id), 'HandleVisibility', 'off',Color = col_CA_raw(id,:));
%         if PARAMETERS_CENT(i).CA_nOutlr(id) == 1
%             plot(TIME, PARAMETERS_CENT(i).CA_cent(:,PARAMETERS_CENT(i).CA_outlr(:,:,id) ,id), LineStyle="--",Color = col_B(id,:) );
%         elseif PARAMETERS_CENT(i).CA_nOutlr(id) > 1
%             idx = find(PARAMETERS_CENT(i).CA_outlr(:,:,id) == 1);
%             plot(TIME, PARAMETERS_CENT(i).CA_cent(:,idx(1),id), LineStyle="--",Color = col_B(id,:) ); 
%             plot(TIME, PARAMETERS_CENT(i).CA_cent(:,PARAMETERS_CENT(i).CA_outlr(:,:,id),id),'HandleVisibility', 'off', LineStyle="--",Color = col_B(id,:) );
%         end
    end

    title("Test C - Vocal FB");
    yline(0, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
%     yline(400, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
    ylim([-400 800]);
    xlabel('Time in seconds');
    ylabel('Cent');
    legend("Effect 0Cent", "Effect 400Cent");
    
%     if sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,1)) > 0 &&  sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,2)) > 0
%         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
%     elseif sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,1)) > 0
%         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
%     elseif sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,2)) > 0
%          legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
%     else
%          legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
%     end
    hold off

    %%%%%%%%%%%%%%%%%%%%%%%% Plot cumulative effect of target and vocal feedback in subplot %%%%%%%%%%%%%%%%%%%%%%

    subplot(3,3,6)
    hold on

    keepIndex = ~isnan(PARAMETERS_CENT(i).CA_mean);
    for id = 1:length(keepIndex(1,:))
        x = TIME(keepIndex(:,id));
        y = PARAMETERS_CENT(i).CA_vocalFeebackEffect(keepIndex(:,id),id) + PARAMETERS_CENT(i).B_sensoryTargetEffect(keepIndex(:,id),id);
%         sd = sqrt(PARAMETERS_CENT(i).CA_variance(keepIndex(:,id),id));
        plot(x(:), y, Color=col_CA(id,:),LineWidth=2)
%         patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_CA(id,:),'FaceAlpha',0.2, 'EdgeColor','none')
%         idx = find(~PARAMETERS_CENT(i).CA_outlr(:,:,id) == 1);
%         plot(TIME,PARAMETERS_CENT(i).CA_cent(:,idx(1),id), Color = col_CA_raw(id,:));
%         plot(TIME,PARAMETERS_CENT(i).CA_cent(:,~PARAMETERS_CENT(i).CA_outlr(:,:,id),id), 'HandleVisibility', 'off',Color = col_CA_raw(id,:));
%         if PARAMETERS_CENT(i).CA_nOutlr(id) == 1
%             plot(TIME, PARAMETERS_CENT(i).CA_cent(:,PARAMETERS_CENT(i).CA_outlr(:,:,id) ,id), LineStyle="--",Color = col_B(id,:) );
%         elseif PARAMETERS_CENT(i).CA_nOutlr(id) > 1
%             idx = find(PARAMETERS_CENT(i).CA_outlr(:,:,id) == 1);
%             plot(TIME, PARAMETERS_CENT(i).CA_cent(:,idx(1),id), LineStyle="--",Color = col_B(id,:) ); 
%             plot(TIME, PARAMETERS_CENT(i).CA_cent(:,PARAMETERS_CENT(i).CA_outlr(:,:,id),id),'HandleVisibility', 'off', LineStyle="--",Color = col_B(id,:) );
%         end
    end

    title("Target + Vocal FB");
    yline(0, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
%     yline(400, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
    ylim([-400 800]);
    xlabel('Time in seconds');
    ylabel('Cent');
    legend("Effect 0Cent", "Effect 400Cent");
    
%     if sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,1)) > 0 &&  sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,2)) > 0
%         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
%     elseif sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,1)) > 0
%         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
%     elseif sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,2)) > 0
%          legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
%     else
%          legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
%     end
    hold off

    %%%%%%%%%%%%%%%%%%%%%%%% Plot Test CB in subplot %%%%%%%%%%%%%%%%%%%%%%
    
    subplot(3,3,7)
    hold on
    keepIndex = ~isnan(PARAMETERS_CENT(i).CB_mean);
    for targ = 1:length(TARGET)
        for shif = 1:length(SHIFT)
            x = TIME(keepIndex(:,:,targ,shif));
            y = PARAMETERS_CENT(i).CB_mean(keepIndex(:,:,targ,shif),:,targ,shif);
            sd = sqrt(PARAMETERS_CENT(i).CB_variance(keepIndex(:,:,targ,shif),:,targ,shif));
            plot(x(:), y, Color=col_CB(shif,:),LineWidth=2)
            patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_CB(shif,:),'FaceAlpha',0.2, 'EdgeColor','none')
        end
    end

%%% Uncomment this code fragment if the rawdata should be plotted for test CB as well.    
%     for targ = 1:length(TARGET)
%         for shif = 1:length(SHIFT)
%             plot(TIME,PARAMETERS_CENT(i).CB_cent(:,~PARAMETERS_CENT(i).CB_outlr(:,:,targ,shif),targ,shif), Color = col_CB_raw(shif,:));
%         end
%     end
%%%

    title("Test C - part B");
    yline(0, LineStyle="--");
    yline(400, LineStyle="--");
    ylim([-400 800])

%%% Uncomment this code fragment if the outliers should be plotted for test CB as well.    
%     for targ = 1:length(TARGET)
%         for shif = 1:length(SHIFT)
%             if sum(PARAMETERS_CENT(i).CB_nOutlr(:,:,targ,shif)) > 0
%                 plot(TIME,PARAMETERS_CENT(i).CB_cent(:,PARAMETERS_CENT(i).CB_outlr(:,:,targ,shif) ,targ,shif), LineStyle="--", Color = col_CB_raw(shif,:) )
%             end
%         end
%     end
%%%
    
    legend("-100","", "-50", "", "-25", "", "0", "", "25","","50","","100", 'Orientation', 'horizontal', 'Location', 'south','FontSize',5);
    xlabel('Time in seconds');
    ylabel('Cent');
    hold off
   
    %%%%%%%%%%%%%%%%%%%%%%%% Plot Test Effect of shifted feedback in subplot %%%%%%%%%%%%%%%%%%%%%%
    
    subplot(3,3,8)
    hold on
    keepIndex = ~isnan(PARAMETERS_CENT(i).CB_mean);
%     for targ = 1:length(TARGET)
    targ = 1;
    for shif = 1:length(SHIFT)
        x = TIME(keepIndex(:,:,targ,shif));
        y = PARAMETERS_CENT(i).CB_shiftedFeedbackEffect(keepIndex(:,:,targ,shif),:,targ,shif);
%             sd = sqrt(PARAMETERS_CENT(i).CB_shiftedFeedbackError(keepIndex(:,:,targ,shif),:,targ,shif));
        plot(x(:), y, Color=col_CB(shif,:),LineWidth=2)
%             patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_CB(shif,:),'FaceAlpha',0.2, 'EdgeColor','none')
    end
%     end

    title("Test C - Shifted FB 0Cent");
    yline(0, LineStyle="--");
%     yline(400, LineStyle="--");
    ylim([-400 800])


%     legend("-100","", "-50", "", "-25", "", "0", "", "25","","50","","100", 'Orientation', 'horizontal', 'Location', 'south','FontSize',5);
    xlabel('Time in seconds');
    ylabel('Cent');
    hold off

    subplot(3,3,9)
    hold on
    keepIndex = ~isnan(PARAMETERS_CENT(i).CB_mean);
    targ = 2;
%     for targ = 1:length(TARGET)
    for shif = 1:length(SHIFT)
        x = TIME(keepIndex(:,:,targ,shif));
        y = PARAMETERS_CENT(i).CB_shiftedFeedbackEffect(keepIndex(:,:,targ,shif),:,targ,shif);
%             sd = sqrt(PARAMETERS_CENT(i).CB_shiftedFeedbackError(keepIndex(:,:,targ,shif),:,targ,shif));
        plot(x(:), y, Color=col_CB(shif,:),LineWidth=2)
%             patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_CB(shif,:),'FaceAlpha',0.2, 'EdgeColor','none')
    end
%     end

    title("Test C - Shifted FB 400Cent");
    yline(0, LineStyle="--");
%     yline(400, LineStyle="--");
    ylim([-400 800])

%     legend("-100","", "-50", "", "-25", "", "0", "", "25","","50","","100", 'Orientation', 'horizontal', 'Location', 'south','FontSize',5);
    xlabel('Time in seconds');
    ylabel('Cent');
    hold off
    %%%%%%%%%%%%%%%%%%%%%%%%% Plot Test D in subplot %%%%%%%%%%%%%%%%%%%%%%

   
%     subplot(2,3,6)
%     hold on
%     keepIndex = ~isnan(PARAMETERS_CENT(i).D_mean);
% 
%     for targ = 1:length(keepIndex(1,:,:,1))
%         for shif = 1:length(keepIndex(1,:,1,:))
%             x = TIME_D(keepIndex(:,:,targ,shif));
%             y = PARAMETERS_CENT(i).D_mean(keepIndex(:,:,targ,shif),:,targ,shif);
%             sd = sqrt(PARAMETERS_CENT(i).D_variance(keepIndex(:,:,targ,shif),:,targ,shif));
%             plot(x(:), y, Color=col_D(shif,:),LineWidth=2)
%             patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_D(shif,:),'FaceAlpha',0.2, 'EdgeColor','none')
%         end
%     end

%%% Uncomment this code fragment if the rawdata and outliers should be plotted for test D as well.    
%     for targ = 1:length(TARGET)
%         for shif = 1:length(SHIFT)
%             plot(TIME_D,PARAMETERS_CENT(i).D_cent(:,~PARAMETERS_CENT(i).D_outlr(:,:,targ,shif),targ,shif),'HandleVisibility', 'off', Color = col_D_raw(shif,:));
%             if PARAMETERS_CENT(i).D_nOutlr(:,:,targ,shif) > 0
%                 plot(TIME_D,PARAMETERS_CENT(i).D_cent(:,PARAMETERS_CENT(i).D_outlr(:,:,targ,shif),targ,shif),'HandleVisibility', 'off', Color = col_D_raw(shif,:), LineStyle="--" );
%             end
%         end
%     end
%%%

% 
%     title("Test D");
%     yline(0, LineStyle="--");
%     yline(400, LineStyle="--");
%     ylim([-400 800])
%     leg = legend("-100","", "-50", "", "-25", "", "0", "", "25","","50","","100",'Orientation', 'horizontal', 'Location','south','FontSize',5);
%     leg.ItemTokenSize = [16,100];
%     xlabel('Time in seconds');
%     ylabel('Cent');

%     hold off
    sgtitle(sbjnames{i})
    fig.WindowState = 'maximized';
    exportgraphics(fig, folder + sbjnames{i} + "_CombinedPlotsNEW" + FORMAT, "Resolution",RESOLUTION,'ContentType',CONTENT_TYPE)
    close(fig);
end
disp("Plot of Rawdata in CENT are DONE & SAVED....")