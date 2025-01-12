% Fit sigmoid function for task CB, separately for 0 and 400 cent targets
% x axis (true) shift in auditory input in cent, y axis actual shift in the
% vocalization relative to the target in cent
% one sigmoid model is fit for all data from all interval
% [b,dev,stats] = glmfit(X,y,'link','logit'), b includes an offset

savepath = 'D:\Pitch adaptation\test\';
file_CENT = savepath + "PARAMETERS_CENT.mat";
file_var = savepath + "variables.mat";
load(file_CENT)
load(file_var)

%%% Set Resolution and Format for Plots
folder = savepath + "\Plots\Sigmoid\";
RESOLUTION  = 300; %Resolution in dots per inch (DPI), specified as a whole number that is greater than or equal to 1.
FORMAT = ".pdf"; %Possible formats are: '.jpg', '.png', '.tif', '.gif', '.pdf', '.emf' or '.eps'
CONTENT_TYPE = 'auto'; % Type of content to store when saving as an EMF, EPS, or PDF file. Specify the value as one of these options:
                         % 'auto' — MATLAB controls whether the content is a vector graphic or an image.
                         % 'vector' — Stores the content as a vector graphic that can scale to any size. If you are saving a PDF file, embeddable fonts are included in the file.
                         % 'image' — Rasterizes the content into one or more images within the file.
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
shifts = [-100, -50, -25, 0, 25, 50, 100];
% tCA = 70;
% tCB = 10;
offset = 100;
%%
sigfit = fittype( @(a, b, c, x) c+(a./(1+exp(-b.*x))) );
tic
data_alltime = [];
% glms = struct([]);
a = nan(nsubj, 2); b = nan(nsubj, 2); c = nan(nsubj, 2); d = nan(nsubj, 2);
slope = nan(nsubj, 2); % slope at x=0
interc = nan(nsubj, 2); % intercept x=0
rmse = nan(nsubj, 2); % rms error of fit
rsq = nan(nsubj, 2); % Degree-of-freedom adjusted coefficient of determination
sse = nan(nsubj, 2);
stds = nan(nsubj, 2); % standard deviation averaged across all shift
dprime = nan(nsubj, 2); % fitted shift at two ends divided by STD
% mean0 = nan(nsubj, 7); var0 = nan(nsubj, 7);
est400 = nan(nsubj, 7);
est0 = nan(nsubj, 7); % estimated value for 0 target
stats_all_sub = struct([]);
for s = 1:nsubj
    data = PARAMETERS_CENT(s);       
    glms = struct([]);
   
    X0 = []; y0 = [];
    X400 = []; y400 = [];
    for t = (80:T)+offset
%         tt = t-offset;
        % zero shift for 0 cent and 400 cent target
%         y_CA0 = zeros(tCA-sum(data.CA_outlr(:, :, 1)), 1);
%         X_CA0 = y_CA0;
%         y_CA0 = data.CA_cent(t, ~data.CA_outlr(:, :, 1), 1)';
%         X_CA0 = zeros(size(y_CA0));
%         X_CA0(:) = 0; 

%         y_CA400 = zeros(tCA-sum(data.CA_outlr(:, :, 2)), 1);
%         X_CA400 = y_CA400;
%         y_CA400 = data.CA_cent(t, ~data.CA_outlr(:, :, 2), 2)'-400;
%         X_CA400 = zeros(size(y_CA400));
%         X_CA400(:) = 0; 
        
        y_CB0 = []; y_CB400 = [];
        X_CB0 = []; X_CB400 = [];
        for target = 1:2
            for shift = 1:7
%                 y_temp = zeros(tCB - sum(data.CB_outlr(:, :, target, shift)), 2);
%                 X_temp = y_temp;
                if target == 1
                    y_temp = data.CB_cent(t, ~data.CB_outlr(:, :, target, shift), target, shift)';
                    X_temp = shifts(shift)*ones(size(y_temp));                    
                    if shift==4 % zero shift
                        y_CA0 = data.CA_cent(t, ~data.CA_outlr(:, :, 1), 1)';
                        X_CA0 = zeros(size(y_CA0));
                        y_temp = [y_temp; y_CA0];
                        X_temp = [X_temp; X_CA0];
                    end 
%                     mean0(tt, shift) = mean(y_temp);
%                     var0(tt, shift) = sqrt(var(y_temp));
                    y_CB0 = [y_CB0; y_temp];
                    X_CB0 = [X_CB0; X_temp];
                else
                    y_temp = data.CB_cent(t, ~data.CB_outlr(:, :, target, shift), target, shift)'-400;
                    X_temp = shifts(shift)*ones(size(y_temp));                    
                    if shift==4 % zero shift
                        y_CA400 = data.CA_cent(t, ~data.CA_outlr(:, :, 2), 2)'-400;
                        X_CA400 = zeros(size(y_CA400));
                        y_temp = [y_temp; y_CA400];
                        X_temp = [X_temp; X_CA400];
                    end
%                     mean400(tt, shift) = mean(y_temp);
%                     var400(tt, shift) = sqrt(var(y_temp));
                    y_CB400 = [y_CB400; y_temp];
                    X_CB400 = [X_CB400; X_temp];
                end
            end
        end
        X0 = [X0; X_CB0]; y0 = [y0; y_CB0];
        X400 = [X400; X_CB400]; y400 = [y400; y_CB400];
    end

    nanid = find(isnan(y0));
    X0(nanid) = []; y0(nanid) = [];
    
    nanid = find(isnan(y400));
    X400(nanid) = []; y400(nanid) = [];

    std0 = 0; std400 = 0;
    for sh = 1:7 % compute variance per shift then average to get std
        shift = shifts(sh);
        idx0 = find(X0==shift);
        idx400 = find(X400==shift);
        std0 = std0 + std(y0(idx0));
        std400 = std400 + std(y400(idx400));
    end
    stds(s, 1) = std0/7; stds(s, 2) = std400/7;
    [f0,gof0,out0] = fit(X0,y0,sigfit, 'StartPoint', [-200, 0, 100]);
    [f400, gof400, out400] = fit(X400, y400, sigfit,  'StartPoint', [-200, 0, 100]);
    cv0 = coeffvalues(f0);
    a(s, 1) = cv0(1); b(s, 1) = cv0(2); c(s, 1) = cv0(3); %d(s, 1) = cv0(4);
    interc(s, 1) = c(s, 1) + a(s, 1)/2;
    slope(s, 1) = logitslope0(a(s, 1), b(s, 1), 0);
    rmse(s, 1) = gof0.rmse;
    sse(s, 1) = gof0.sse;
    rsq(s, 1) = gof0.adjrsquare;
    est0(s, :) = feval(f0, shifts)';
    dprime(s, 1) = abs((est0(s, 1)-est0(s, 7))/stds(s, 1));

    cv400 = coeffvalues(f400);
    a(s, 2) = cv400(1); b(s, 2) = cv400(2); c(s, 2) = cv400(3); %d(s, 2) = cv400(4);
    interc(s, 2) = c(s, 2) + a(s, 2)/2;
    slope(s, 2) = logitslope0(a(s, 2), b(s, 2), 0);
    rmse(s, 2) = gof400.rmse;
    sse(s, 2) = gof400.sse;
    rsq(s, 2) = gof400.adjrsquare;
    est400(s, :) = feval(f400, shifts)';
    dprime(s, 2) = abs((est400(s, 1)-est400(s, 7))/stds(s, 2));
    % plot for individual subject

%     fig = figure();
%     sp = subplot(121); % data and fit for target 0
%     hold on
%     plot(f0, X0, y0);
%     yline(0, LineStyle="--");
%     xlabel('Input shift (cent)', 'FontSize',15)
%     ylabel('Vocalization shift (cent)', 'FontSize',15)
%     title(['Target 0, a=' num2str(a(s, 1)) ', b=' num2str(b(s, 1)) ', bias=' num2str(bias(s, 1))], ...
%         'FontSize',15)
%     sp.FontSize = 15;
%     yl = [-200, 200];
%     ylim(sp, yl)
%     clear ylim
% 
%     sp = subplot(122); % data and fit for target 0
%     hold on
%     plot(f400, X400, y400)
%     yline(0, LineStyle="--");
%     xlabel('Input shift (cent)', 'FontSize',15)
%     ylabel('Vocalization shift (cent)', 'FontSize',15)
%     title(['Target 400, a=' num2str(a(s, 2)) ', b=' num2str(b(s, 2)) ', bias=' num2str(bias(s, 2))], ...
%          'FontSize',15)
%     sp.FontSize = 15;
%     yl = [-200, 200];
%     ylim(sp, yl);
%     clear ylim
% 
%     sgtitle(data.Sbjnames)
%     fig.WindowState = 'maximized';
%     exportgraphics(fig, folder + data.Sbjnames + "_SigmoidTaskC" + FORMAT, "Resolution",RESOLUTION,'ContentType',CONTENT_TYPE)
%     close(fig);
    
end
toc
%% plot summary
ol = [30, 36, 47];
nol = 1:nsubj;
nol(ol) = [];
figure
sp = subplot(221);
plot(slope(nol, 1), slope(nol, 2), "Marker", 'o', 'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1)
hold on
plot(slope(ol, 1), slope(ol, 2), "Marker", '*','Color', 'r', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1)
hold on
yl = [-1.5 0.5];
xlim(sp, yl);
ylim(sp, yl);
lim = -1.5:0.1:0.5;
plot(lim, lim, '--k')
xlabel('0 Cent Target')
ylabel('400 Cent Target')
title('Slope at 0 Shift')
sp.FontSize = 15;


sp = subplot(222);
plot(rmse(nol, 1), rmse(nol, 2), "Marker", 'o', 'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1)
hold on
plot(rmse(ol, 1), rmse(ol, 2), "Marker", '*','Color', 'r', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1)
yl = [0 200];
xlim(sp, yl);
ylim(sp, yl);
lim = 0:10:200;
plot(lim, lim, '--k')
xlabel('0 Cent Target')
ylabel('400 Cent Target')
title('RMS Error')
sp.FontSize = 15;

sp = subplot(234);
plot(a(nol, 1), a(nol, 2), "Marker", 'o', 'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1)
hold on
plot(a(ol, 1), a(ol, 2), "Marker", '*','Color', 'r', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1)
yl = [-7000 0];
xlim(sp, yl);
ylim(sp, yl);
lim = -7000:1000:0;
plot(lim, lim, '--k')
xlabel('0 Cent Target')
ylabel('400 Cent Target')
title('Parameter a')
sp.FontSize = 15;

sp = subplot(235);
plot(b(nol, 1), b(nol, 2), "Marker", 'o', 'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1)
hold on
plot(b(ol, 1), b(ol, 2), "Marker", '*','Color', 'r', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1)
% lim = -1.5:0.1:0.5;
% plot(lim, lim, '--k')
yl = [-0.2 0.101];
xlim(sp, yl);
ylim(sp, yl);
lim = -0.2:0.1:0.1;
plot(lim, lim, '--k')
xlabel('0 Cent Target')
ylabel('400 Cent Target')
title('Parameter b')
sp.FontSize = 15;

sp = subplot(236);
plot(interc(nol, 1), interc(nol, 2), "Marker", 'o', 'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1)
hold on
plot(interc(ol, 1), interc(ol, 2), "Marker", '*','Color', 'r', ...
    'LineStyle','none', 'MarkerSize', 8, 'LineWidth', 1)
yl = [-500 200];
xlim(sp, yl);
ylim(sp, yl);
lim = -500:10:200;
plot(lim, lim, '--k')
xlabel('0 Cent Target')
ylabel('400 Cent Target')
title('Bias')
sp.FontSize = 15;

toc
function s = logitslope0(a, b, x)
d = (1+exp(-b.*x)).^2;
s = a.*b.*exp(-b.*x)./d;
end


