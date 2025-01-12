
%%% Creates one plot per participant for test B and saves them in the Folder "\Plots\testB\" created at the defined "savepath" directory.
%%% NEEDS PARAMETERS: "FORMAT", "savepath", "TIME", "RESOLUTION", 
%%% NEEDS FILE: PARAMETERS_CENTS.


% Author: Philipp Eugster <eugsteph@student.ethz.ch>
% Created: 10.12.2022


%%% Specify Colours
col_B = [[0 0.4470 0.7410];[0.4940 0.1840 0.5560];];
col_B_raw = [[0 0.4470 0.7410 0.3];[0.4940 0.1840 0.5560 0.3];];

%%% Create Directory to save Plots
folder = savepath + "\Plots\testB\";
if ~exist(folder, 'dir')
    mkdir(folder);
end


disp("START generating plots for TEST_B...")


for i = 1:nsubj
    fig = figure(i);
    hold on

    %%% Using the patch function below only works with non-nan values. Save indexes of non-NaN values in "keepIndex"
    keepIndex = ~isnan(PARAMETERS_CENT(i).B_mean);

    %%% for loop to plot data for both targets
    for id = 1:length(keepIndex(1,:))

        %%% Plots Mean +/- SD as a shaded area.
        x = TIME(keepIndex(:,id));
        y = PARAMETERS_CENT(i).B_mean(keepIndex(:,id),id);
        sd = sqrt(PARAMETERS_CENT(i).B_variance(keepIndex(:,id),id));
        plot(x(:), y, Color=col_B(id,:),LineWidth=2)
        patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_B(id,:),'FaceAlpha',0.2, 'EdgeColor','none')
        
        %%% Plot the rawdata as well (Only the first trial is plotted with "HandleVisibility" = On, so that it only shows up once in the legend.
        idx = find(~PARAMETERS_CENT(i).B_outlr(:,:,id) == 1);
        plot(TIME,PARAMETERS_CENT(i).B_cent(:,idx(1),id), Color = col_B_raw(id,:));
        plot(TIME,PARAMETERS_CENT(i).B_cent(:,~PARAMETERS_CENT(i).B_outlr(:,:,id),id), 'HandleVisibility', 'off',Color = col_B_raw(id,:));

        %%% Depending on the number of outliers per target the outliers need to be plotted with 'HandleVisibility' set to 'off' so the outliers don't show up in the legend mutlible times.
        if PARAMETERS_CENT(i).B_nOutlr(id) == 1
            plot(TIME, PARAMETERS_CENT(i).B_cent(:,PARAMETERS_CENT(i).B_outlr(:,:,id) ,id), LineStyle="--",Color = col_B(id,:) );
        elseif PARAMETERS_CENT(i).B_nOutlr(id) > 1
            idx = find(PARAMETERS_CENT(i).B_outlr(:,:,id) == 1);
            plot(TIME, PARAMETERS_CENT(i).B_cent(:,idx(1),id), LineStyle="--",Color = col_B(id,:) ); 
            plot(TIME, PARAMETERS_CENT(i).B_cent(:,PARAMETERS_CENT(i).B_outlr(:,:,id),id),'HandleVisibility', 'off', LineStyle="--",Color = col_B(id,:) );
        end
    end

    %%% Plot reference line for the two targets
    yline(0, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
    yline(400, '-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
    
    %%% Generate legend according if / how many outliers are present.
    if sum(PARAMETERS_CENT(i).B_nOutlr(:,:,1)) > 0 &&  sum(PARAMETERS_CENT(i).B_nOutlr(:,:,2)) > 0
        legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    elseif sum(PARAMETERS_CENT(i).B_nOutlr(:,:,1)) > 0
        legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Outlier 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    elseif sum(PARAMETERS_CENT(i).B_nOutlr(:,:,2)) > 0
         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent","Outlier 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    else
         legend("Mean 0Cent","STDEV 0Cent","Rawdata 0Cent","Mean 400Cent","STDEV 400Cent","Rawdata 400Cent",'Location','best','FontSize',5,'NumColumns',2)
    end
    
    %%% Generate axis labels and limits
    title(strcat(sbjnames(i), " - ", "Test B"));
    ylim([-400 800])
    xlabel('Time in seconds');
    ylabel('Cent');
    hold off
    
    %%% Export and close figure.
    exportgraphics(fig, folder + sbjnames{i} + "_TestB" + FORMAT, "Resolution",RESOLUTION,'ContentType',CONTENT_TYPE)
    close(figure(i));
end

disp("Plot of Rawdata in CENTS are DONE & SAVED....")









