% Author: Philipp Eugster
% Created: 21.12.2022

%%% NEEDS
% FILE: PARAMETERS_CENT.mat


%%% Suddend and unnatural jumps in pitch, wich almost only happen at the very end of a trial are removed manually here. Used same criteria as descripted in the follwoing paper:
%https://www.biorxiv.org/content/10.1101/2020.06.06.138263v2.full.pdf (Title: Spontaneous variability predicts adaptive motor response in vocal pitch control, Ryosuke O. Tachibana)






%%% Test CB (only checked if between 1second to 1.5 second no unnatrual outliers.

idx = find("f05" == sbjnames);
PARAMETERS_CENT(idx).CB_HZ(150:end,7,1,4) = nan;


idx = find("f12NC" == sbjnames);
PARAMETERS_CENT(idx).CB_HZ(120:end,8,2,3) = nan;
PARAMETERS_CENT(idx).CB_HZ(130:end,9,2,3) = nan;
PARAMETERS_CENT(idx).CB_HZ(130:end,10,2,4) = nan;



idx = find("f26NC" == sbjnames);
PARAMETERS_CENT(idx).CB_HZ(128:131,8,2,5) =  347.8319; %Unnatural jump is corrected with value befor the jump


idx = find("m21NC" == sbjnames);
PARAMETERS_CENT(idx).CB_HZ(140:end,2,1,4) = nan;


% t = 4;
% 
% for i = 1:10
%     A = ischange(PARAMETERS_CENT(idx).CB_cent(60:180,i,1,t), 'Threshold', 100000/20);
%     if PARAMETERS_CENT(idx).CB_outlr(:,i,1,t) == 1
%         continue
%     end
%     if sum(A) > 1
%         figure();
%         plot(TIME(60:180), PARAMETERS_CENT(idx).CB_cent(60:180,i,1,t));
%         xlim([0 2])
%         ylim([-400 800])
%         title(sprintf("Rapid change in Col %d and target %d \n", i,1))
%         fprintf("Rapid change in Col %d and target %d \n", i,1);
%     end
% end
% 
% 
% for i = 1:10
%     A = ischange(PARAMETERS_CENT(idx).CB_cent(60:180,i,2,t), 'Threshold', 10000);
%     if PARAMETERS_CENT(idx).CB_outlr(:,i,2,t) == 1
%         continue
%     end
%     if sum(A) > 1
%         figure();
%         plot(TIME(60:180), PARAMETERS_CENT(idx).CB_cent(60:180,i,2,t));
%         xlim([0 2])
%         ylim([-400 800])
%         title(sprintf("Rapid change in Col %d and target %d \n", i,2))
%         fprintf("Rapid change in Col %d and target %d \n", i,2);
%     end
% end