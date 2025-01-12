
%%% Creates one plot per participant for test A and saves them in the Folder "\Plots\testA\" created at the defined "savepath" directory.
%%% NEEDS PARAMETERS: "FORMAT", "savepath", "TIME", "RESOLUTION", 
%%% NEEDS FILE: PARAMETERS_CENTS.


% Author: Philipp Eugster <eugsteph@student.ethz.ch>
% Created: 10.12.2022


%%% Specify Colours
col_A = [0 0.4470 0.7410];
col_A_raw = [0 0.4470 0.7410 0.3];

%%% Create Directory to save Plots
folder = savepath + "\Plots\testA\";
if ~exist(folder, 'dir')
    mkdir(folder);
end

disp("START generating plots for Test A...")


for i = 1:nsubj
    fig = figure(i);

    %%% Using the patch function below only works with non-nan values. Save indexes of non-NaN values in "keepIndex"
    keepIndex = ~isnan(PARAMETERS_CENT(i).A_mean);

    %%% Plots Mean +/- SD as a shaded area.
    x = TIME(keepIndex);                                                       
    y = PARAMETERS_CENT(i).A_mean(keepIndex);                                 
    sd = sqrt(PARAMETERS_CENT(i).A_variance(keepIndex));
    plot(x(:), y, Color=col_A,LineWidth=2)
    hold on
    patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_A,'FaceAlpha',0.2, 'EdgeColor','none')
    

    %%% Plot the rawdata as well (Only the first trial is plotted with "HandleVisibility" = On, so that it only shows up once in the legend.
    idx = find(~PARAMETERS_CENT(i).A_outlr == 1);
    plot(TIME,PARAMETERS_CENT(i).A_cent(:,idx(1)), Color = col_A_raw);
    plot(TIME,PARAMETERS_CENT(i).A_cent(:,~PARAMETERS_CENT(i).A_outlr), 'HandleVisibility', 'off',Color = col_A_raw);

    %%% Reference line at 0 cent.
    yline(0,'-.b','HandleVisibility','off',LineWidth=0.8, Color="k");
    
    %%% Plot outliers if present (Incluences legend).
    if PARAMETERS_CENT(i).A_nOutlr > 0
        plot(TIME, PARAMETERS_CENT(i).A_cent(:,PARAMETERS_CENT(i).A_outlr), Color = col_A, LineStyle='--');
        legend("Mean","STDEV", "Rawdata", "Outlier",'FontSize',5 );
    else
        legend("Mean","STDEV", "Rawdata",'FontSize',5 );
    end
    
    %%% Set label names and limits
    title(strcat(sbjnames(i), " - ", "Test A"));
    xlabel('Time in seconds');
    ylabel('Cent');
    ylim([-400 400])
    xlim([0 2])
    hold off
    
    %%% Export and close figure
    exportgraphics(fig, folder + sbjnames{i} + "_TestA" + FORMAT, "Resolution",RESOLUTION,'ContentType',CONTENT_TYPE)
    close(figure(i));
end
disp("Plot of Rawdata in CENT are DONE & SAVED....")


