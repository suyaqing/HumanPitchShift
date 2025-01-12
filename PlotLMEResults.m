% plot summary LME results
clear
savepath = 'D:\Pitch adaptation\test\';
fn = savepath + "MaxiGLMM_RandIntSubj.mat";
load(fn);
% first construct matrices for plotting
T = length(glms);
trend_mean = nan(T, 1); trend_SE = trend_mean;
target_mean = nan(T, 2); target_SE = target_mean;
FB_mean = nan(T, 2); FB_SE = FB_mean;
shift_mean = nan(T, 6); shift_SE = shift_mean;

for t = 1:T
    glme = glms(t).glme;
    if ~isempty(glme)
        str = dataset2struct(glme.Coefficients);
        trend_mean(t) = str(1).Estimate; trend_SE(t) = str(1).SE;
        for targ = 1:2
            target_mean(t, targ) = str(targ+1).Estimate;
            target_SE(t, targ) = str(targ+1).SE;
            FB_mean(t, targ) = str(targ+3).Estimate;
            FB_SE(t, targ) = str(targ+3).SE;
        end
        for shft = 1:6
            shift_mean(t, shft) = str(shft+5).Estimate;
            shift_SE(t, shft) = str(shft+5).SE;
        end
    end
end
offset = 100;
% save(fn, "trend_mean", "trend_SE", "target_mean", "target_SE", "FB_mean", ...
%     "FB_SE", "shift_mean", "shift_SE", '-append')
%%

tt = ([1:T]+(100-offset))/100 - 1/100; tt = tt(:);
col_A = [0 0.4470 0.7410];
col_B = [[0 0.4470 0.7410];[0.4940 0.1840 0.5560];]; 
col_CB = [[0.4940 0.1840 0.5560];[0 0.4470 0.7410];[0.3010 0.7450 0.9330];[0.4660 0.6740 0.1880];[0.9290 0.6940 0.1250];[0.8500 0.3250 0.0980];[0.6350 0.0780 0.1840];]; %RGB color code, nice for 7 shifts.

% subplot 1: pitch trend (interception)
figure
ax = subplot(221);
% subplot(221)
hold on
plot(tt, trend_mean, Color=col_A, LineWidth=2)
id = ~isnan(trend_mean);
patch([tt(id); flipud(tt(id))], [trend_mean(id)-trend_SE(id);  flipud(trend_mean(id)+trend_SE(id))], col_A,'FaceAlpha',0.2, 'EdgeColor','none')
yline(0, LineStyle="--");
ylim([-100 100])
title("Vocal Trend", 'FontSize',15)
xlabel('Time (s)', 'FontSize',15)
ylabel('Magnitude (cent)', 'FontSize',15)
legend("Mean","SE", 'FontSize',12);
% ax = gca;
% ax.FontSize = 15;
% ax.FontSize = 15;
set(ax, 'Fontsize', 12)
% subplot 2: target
ax = subplot(222);
hold on
for targ = 1:2
    plot(tt, target_mean(:, targ), Color=col_B(targ, :), LineWidth=2)
    patch([tt(id); flipud(tt(id))], [target_mean(id, targ)-target_SE(id, targ); ...
        flipud(target_mean(id, targ)+target_SE(id, targ))], ...
        col_B(targ, :),'FaceAlpha',0.2, 'EdgeColor','none')
end
ylim([-100 100])
yline(0, LineStyle="--");
% yline(400, LineStyle="--");
title("Auditory Target", 'FontSize',15)
xlabel('Time (s)', 'FontSize',15)
ylabel('Magnitude (cent)', 'FontSize',15)
legend("Mean 0 Cent","SE 0 Cent", "Mean 400 Cent","SE 400 Cent",'Location','best','FontSize',12,'NumColumns',2)
set(ax, 'Fontsize', 12)
% subplot 3: feedback
ax = subplot(223);
hold on
for targ = 1:2
    plot(tt, FB_mean(:, targ), Color=col_B(targ, :), LineWidth=2)
    patch([tt(id); flipud(tt(id))], [FB_mean(id, targ)-FB_SE(id, targ); flipud(FB_mean(id, targ)+FB_SE(id, targ))], ...
        col_B(targ, :),'FaceAlpha',0.2, 'EdgeColor','none')
end
yline(0, LineStyle="--");
ylim([-100 100])
title("Auditory Feedback", 'FontSize',15)
xlabel('Time (s)', 'FontSize',15)
ylabel('Magnitude (cent)', 'FontSize',15)
legend("Mean 0 Cent","SE 0 Cent", "Mean 400 Cent","SE 400 Cent",'Location','best','FontSize',12,'NumColumns',2)
set(ax, 'Fontsize', 12)
% subplot 4: shift
ax = subplot(224);
hold on
for s = 1:6
    plot(tt, shift_mean(:, s), Color=col_CB(s, :), LineWidth=2)
    patch([tt(id); flipud(tt(id))], [shift_mean(id, s)-shift_SE(id, s); ...
        flipud(shift_mean(id, s)+shift_SE(id, s))], ...
        col_CB(s, :),'FaceAlpha',0.2, 'EdgeColor','none')
end
yline(0, LineStyle="--");
ylim([-100 100])
title("Shifted Auditory Feedback", 'FontSize',15)
xlabel('Time (s)', 'FontSize',15)
ylabel('Magnitude (cent)', 'FontSize',15)
legend("-100","", "-50", "", "-25", "", "25","","50","","100", 'Orientation', 'horizontal', 'Location', 'south','FontSize',12);
set(ax, 'Fontsize', 12)
%% plot auxilliary stats, including random effect sizes, rsq and AIC
AIC = zeros(T, 1);
R2 = zeros(T, 1);
RandEffect = zeros(T, 1);
for t = 1:T
    crit = dataset2struct(glms(t).crit);
    re = dataset2struct(glms(t).random);
    AIC(t) = crit.AIC;
    RandEffect(t) = re.Estimate;
    R2(t) = glms(t).rsq.Adjusted;    
end
figure;
subplot(131)
plot(tt, AIC, LineWidth=2);
title("AIC", 'FontSize',15)
xlabel('Time (s)', 'FontSize',15)

subplot(132)
plot(tt, R2, LineWidth=2);
title("R Squared", 'FontSize',15)
xlabel('Time (s)', 'FontSize',15)

subplot(133)
plot(tt, RandEffect, LineWidth=2);
title("Random Effect STD", 'FontSize',15)
xlabel('Time (s)', 'FontSize',15)
