
%%% Creates one plot per participant for test D and saves them in the Folder "\Plots\testD\" created at the defined "savepath" directory.
%%% NEEDS PARAMETERS: "FORMAT", "savepath", "TIME_D", "RESOLUTION", 
%%% NEEDS FILE: PARAMETERS_CENTS.


% Author: Philipp Eugster <eugsteph@student.ethz.ch>
% Created: 10.12.2022


%%% Specify Colours
col_D = [[0.4940 0.1840 0.5560];[0 0.4470 0.7410];[0.3010 0.7450 0.9330];[0.4660 0.6740 0.1880];[0.9290 0.6940 0.1250];[0.8500 0.3250 0.0980];[0.6350 0.0780 0.1840];]; %RGB color code, nice for 7 shifts.
col_D_raw = [[0.4940 0.1840 0.5560 0.2];[0 0.4470 0.7410 0.2];[0.3010 0.7450 0.9330 0.2];[0.4660 0.6740 0.1880 0.2];[0.9290 0.6940 0.1250 0.2];[0.8500 0.3250 0.0980 0.2];[0.6350 0.0780 0.1840 0.2];]; %RGB color code, nice for 7 shifts.


%%% Create Directory to save Plots
folder = savepath + "\Plots\testD\";
if ~exist(folder, 'dir')
    mkdir(folder);
end



disp("START generating plots for Test D...")

for i = 1:nsubj
    fig = figure(i);
    hold on

    %%% Using the patch function below only works with non-nan values. Save indexes of non-NaN values in "keepIndex"
    keepIndex = ~isnan(PARAMETERS_CENT(i).D_mean);

    %%% Plots Mean +/- SD as a shaded area.
    for targ = 1:length(keepIndex(1,:,:,1))
        for shif = 1:length(keepIndex(1,:,1,:))
            x = TIME_D(keepIndex(:,:,targ,shif));
            y = PARAMETERS_CENT(i).D_mean(keepIndex(:,:,targ,shif),:,targ,shif);
            sd = sqrt(PARAMETERS_CENT(i).D_variance(keepIndex(:,:,targ,shif),:,targ,shif));
            plot(x(:), y, Color=col_D(shif,:),LineWidth=2)
            patch([x(:); flipud(x(:))], [y-sd(:);  flipud(y+sd(:))], col_D(shif,:),'FaceAlpha',0.2, 'EdgeColor','none')
        end
    end

    %%% Comment this code snipped out, if the Rawdata / Outliers should not be shown in the plot.
    for targ = 1:length(TARGET)
        for shif = 1:length(SHIFT)
            plot(TIME_D,PARAMETERS_CENT(i).D_cent(:,~PARAMETERS_CENT(i).D_outlr(:,:,targ,shif),targ,shif),'HandleVisibility', 'off', Color = col_D_raw(shif,:));
            if PARAMETERS_CENT(i).D_nOutlr(:,:,targ,shif) > 0
                plot(TIME_D,PARAMETERS_CENT(i).D_cent(:,PARAMETERS_CENT(i).D_outlr(:,:,targ,shif),targ,shif),'HandleVisibility', 'off', Color = col_D_raw(shif,:), LineStyle="--" );
            end
        end
    end
    %%%
    
    %%% Set axis names and limits.
    title(strcat(sbjnames(i), " - ", "Test D"));
    yline(0, LineStyle="--");
    yline(400, LineStyle="--");
    ylim([-400 800])
    leg = legend("-100","", "-50", "", "-25", "", "0", "", "25","","50","","100",'Orientation', 'horizontal', 'Location','south','FontSize',5);
    leg.ItemTokenSize = [16,100];
    xlabel('Time in seconds');
    ylabel('Cent');
    hold off

    %%% Save & close figure.
    exportgraphics(fig, folder + sbjnames{i} + "_TestD" + FORMAT, "Resolution",RESOLUTION,'ContentType',CONTENT_TYPE);
    close(figure(i));
end
disp("Plots of Rawdata in for Test D are DONE & SAVED....")

