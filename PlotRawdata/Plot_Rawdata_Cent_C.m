%%% Creates one plot per participant for test CA and CB and saves them in the Folder "\Plots\testCA\" and "\Plots\testCB\" created at the defined "savepath" directory.
%%% NEEDS PARAMETERS: "FORMAT", "savepath", "TIME", "RESOLUTION", 
%%% NEEDS FILE: PARAMETERS_CENTS.


% Author: Philipp Eugster <eugsteph@student.ethz.ch>
% Created: 10.12.2022


%%% Set color for test CA 
col_CA = [[0 0.4470 0.7410];[0.4940 0.1840 0.5560];];
col_CA_raw =  [[0 0.4470 0.7410 0.1];[0.4940 0.1840 0.5560 0.1];];

%%% Create Directory to save Plots for test CA
folder = savepath + "\Plots\testCA\";
if ~exist(folder, 'dir')
    mkdir(folder);
end


disp("START generating plots for Test C part A...")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test_CA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nsubj
    
    fig = figure(i);
    hold on
    
    %%% Using the patch function below only works with non-nan values. Save indexes of non-NaN values in "keepIndex"
    keepIndex = ~isnan(PARAMETERS_CENT(i).CA_mean);

    %%% for loop to plot data for both targets
    for id = 1:length(keepIndex(1,:))

        %%% Plots Mean +/- SD as a shaded area.
        x = TIME(keepIndex(:,id));
        y = PARAMETERS_CENT(i).CA_mean(keepIndex(:,id),id);
        sd = sqrt(PARAMETERS_CENT(i).CA_variance(keepIndex(:,id),id));
        plot(x(:), y, Color=col_CA(id,:),LineWidth=2)
        patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_CA(id,:),'FaceAlpha',0.2, 'EdgeColor','none')

        %%% Plot the rawdata as well (Only the first trial is plotted with "HandleVisibility" = On, so that it only shows up once in the legend.
        idx = find(~PARAMETERS_CENT(i).CA_outlr(:,:,id) == 1);
        plot(TIME,PARAMETERS_CENT(i).CA_cent(:,idx(1),id), Color = col_CA_raw(id,:));
        plot(TIME,PARAMETERS_CENT(i).CA_cent(:,~PARAMETERS_CENT(i).CA_outlr(:,:,id),id), 'HandleVisibility', 'off',Color = col_CA_raw(id,:));

        %%% Depending on the number of outliers per target the outliers need to be plotted with 'HandleVisibility' set to 'off' so the outliers don't show up in the legend mutlible times.
        if PARAMETERS_CENT(i).CA_nOutlr(id) == 1
            plot(TIME, PARAMETERS_CENT(i).CA_cent(:,PARAMETERS_CENT(i).CA_outlr(:,:,id) ,id), LineStyle="--",Color = col_B(id,:) );
        elseif PARAMETERS_CENT(i).CA_nOutlr(id) > 1
            idx = find(PARAMETERS_CENT(i).CA_outlr(:,:,id) == 1);
            plot(TIME, PARAMETERS_CENT(i).CA_cent(:,idx(1),id), LineStyle="--",Color = col_B(id,:) ); 
            plot(TIME, PARAMETERS_CENT(i).CA_cent(:,PARAMETERS_CENT(i).CA_outlr(:,:,id),id),'HandleVisibility', 'off', LineStyle="--",Color = col_B(id,:) );
        end
    end
    
    %%% Plot reference line for the two targets
    yline(0, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
    yline(400, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
    
    %%% Generate legend according if / how many outliers are present.
   if sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,1)) > 0 &&  sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,2)) > 0
        legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    elseif sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,1)) > 0
        legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    elseif sum(PARAMETERS_CENT(i).CA_nOutlr(:,:,2)) > 0
         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    else
         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    end

    %%% Generate axis labels and limits
    title(strcat(sbjnames(i), " - ", "Test CA"));
    ylim([-400 800]);
    xlabel('Time in seconds');
    ylabel('Cent');
    hold off

    %%% Export and close figure.
    exportgraphics(fig, folder + sbjnames{i} + "_TestCA" + FORMAT, "Resolution",RESOLUTION);
    close(figure(i));
end
disp("Plots of Rawdata in for Test C part A are DONE & SAVED....")





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Test_CB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Set color for test CB 
col_CB = [[0.4940 0.1840 0.5560];[0 0.4470 0.7410];[0.3010 0.7450 0.9330];[0.4660 0.6740 0.1880];[0.9290 0.6940 0.1250];[0.8500 0.3250 0.0980];[0.6350 0.0780 0.1840];]; %RGB color code, nice for 7 shifts.
col_CB_raw = [[0.4940 0.1840 0.5560 0.2];[0 0.4470 0.7410 0.2];[0.3010 0.7450 0.9330 0.2];[0.4660 0.6740 0.1880 0.2];[0.9290 0.6940 0.1250 0.2];[0.8500 0.3250 0.0980 0.2];[0.6350 0.0780 0.1840 0.2];]; %RGB color code, nice for 7 shifts.


%%% Create Directory to save Plots for test CB
folder = savepath + "\Plots\testCB\";
if ~exist(folder, 'dir')
    mkdir(folder);
end



disp("START generating plots for Test C part B...")


for i = 1:nsubj
    fig = figure(i);
    hold on

    %%% Using the patch function below only works with non-nan values. Save indexes of non-NaN values in "keepIndex"
    keepIndex = ~isnan(PARAMETERS_CENT(i).CB_mean);

    %%% for loop to plot data for both targets and the seven shift conditions
    for targ = 1:length(TARGET)
        for shif = 1:length(SHIFT)
            %%% Plots Mean +/- SD as a shaded area.
            x = TIME(keepIndex(:,:,targ,shif));
            y = PARAMETERS_CENT(i).CB_mean(keepIndex(:,:,targ,shif),:,targ,shif);
            sd = sqrt(PARAMETERS_CENT(i).CB_variance(keepIndex(:,:,targ,shif),:,targ,shif));
            plot(x(:), y, Color=col_CB(shif,:),LineWidth=2)
            patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_CB(shif,:),'FaceAlpha',0.2, 'EdgeColor','none')
        end
    end
    
    %%% Plot raw data and outliers
    for targ = 1:length(TARGET)
        for shif = 1:length(SHIFT)
            plot(TIME,PARAMETERS_CENT(i).CB_cent(:,~PARAMETERS_CENT(i).CB_outlr(:,:,targ,shif),targ,shif),'HandleVisibility', 'off', Color = col_CB_raw(shif,:));
            if PARAMETERS_CENT(i).CB_nOutlr(:,:,targ,shif) > 0
                plot(TIME,PARAMETERS_CENT(i).CB_cent(:,PARAMETERS_CENT(i).CB_outlr(:,:,targ,shif),targ,shif),'HandleVisibility', 'off', Color = col_CB_raw(shif,:), LineStyle="--" );
            end
        end
    end
    
    %%% Plot reference line for the two targets
    yline(0, LineStyle="--");
    yline(400, LineStyle="--");

    %%% Generate axis labels and limits and set the legend.
    title(strcat(sbjnames(i), " - ", "Test CB"));
    ylim([-400 800])
    legend("-100","", "-50", "", "-25", "", "0", "", "25","","50","","100", 'Orientation', 'horizontal', 'Location', 'south','FontSize',5);
    xlabel('Time in seconds');
    ylabel('Cent');
    hold off

    %%% Export and close figure.
    exportgraphics(fig, folder + sbjnames{i} + "_TestCB" + FORMAT, "Resolution",RESOLUTION,'ContentType',CONTENT_TYPE)
    close(figure(i));
end


disp("Plots of Rawdata in for Test C part B are DONE & SAVED....")



