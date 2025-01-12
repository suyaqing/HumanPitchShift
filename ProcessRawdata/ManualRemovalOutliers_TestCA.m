% Author: Philipp Eugster
% Created: 21.12.2022

%%% NEEDS
% FILE: PARAMETERS_CENT.mat


%%% Suddend and unnatural jumps in pitch, wich almost only happen at the very end of a trial are removed manually here. Used same criteria as descripted in the follwoing paper:
%https://www.biorxiv.org/content/10.1101/2020.06.06.138263v2.full.pdf (Title: Spontaneous variability predicts adaptive motor response in vocal pitch control, Ryosuke O. Tachibana)






%%% Test CA

idx = find("f03NC" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(167:end,53,1) = nan;
PARAMETERS_CENT(idx).CA_HZ(149:end,61,1) = nan;

idx = find("f05" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(163:170,28,1) = nan;
PARAMETERS_CENT(idx).CA_HZ(166:170,59,1) = nan;

idx = find("f06" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(176:190,37,1) = nan;

idx = find("f07" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(185:190,12,1) = nan;

idx = find("f08" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(171:190,58,1) = nan;
PARAMETERS_CENT(idx).CA_HZ(181:190,65,1) = nan;

idx = find("f11NC" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(142:150,8,1) = nan;

idx = find("f15" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(72:79,16,1) =   214.4879;
PARAMETERS_CENT(idx).CA_HZ(160:end,25,1) = nan;

idx = find("f17" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(90:end,26,1) = nan;


idx = find("f21NC" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(140:end,9,1) = nan;
PARAMETERS_CENT(idx).CA_HZ(172:end,12,1) = nan;
PARAMETERS_CENT(idx).CA_HZ(151:end,56,2) = nan;


idx = find("f26NC" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(182:end,61,2) = nan;
PARAMETERS_CENT(idx).CA_HZ(173:end,62,2) = nan;

idx = find("m02NC" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(163:end,12,1) = nan;
PARAMETERS_CENT(idx).CA_HZ(175:end,20,1) = nan;


idx = find("m06" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(152:end,56,1) = nan;
PARAMETERS_CENT(idx).CA_HZ(164:end,33,1) = nan;
PARAMETERS_CENT(idx).CA_HZ(170:end,29,1) = nan;
PARAMETERS_CENT(idx).CA_HZ(150:end,28,1) = nan;
PARAMETERS_CENT(idx).CA_HZ(134:end,21,2) = nan;
PARAMETERS_CENT(idx).CA_HZ(150:end,22,2) = nan;
PARAMETERS_CENT(idx).CA_HZ(150:end,25,2) = nan;
PARAMETERS_CENT(idx).CA_HZ(150:end,53,2) = nan;
PARAMETERS_CENT(idx).CA_HZ(148:end,47,2) = nan;


idx = find("m08" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(160:end,52,2) = nan;


idx = find("m09NC" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(180:end,41,1) = nan;

idx = find("m13NC" == sbjnames);



idx = find("m15" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(120:end,66,1) = nan;


idx = find("m24NC" == sbjnames);
PARAMETERS_CENT(idx).CA_HZ(170:end,56,1) = nan;

%%% Used to detect ouliers... Was done semi automatic...

% for i = 1:70
%     A = ischange(PARAMETERS_CENT(idx).CA_cent(40:end,i,1), 'Threshold', 100000/5);
%     if PARAMETERS_CENT(idx).CA_outlr(:,i,1) == 1
%         continue
%     end
%     if sum(A) > 1
%         figure();
%         plot(TIME(40:end), PARAMETERS_CENT(idx).CA_cent(40:end,i,1));
%         xlim([0 2])
%         ylim([-400 800])
%         title(sprintf("Rapid change in Col %d and target %d \n", i,1))
%         fprintf("Rapid change in Col %d and target %d \n", i,1);
%     end
% end
% 
% 
% for i = 1:70
%     A = ischange(PARAMETERS_CENT(idx).CA_cent(40:end,i,2), 'Threshold', 10000/2);
%     if PARAMETERS_CENT(idx).CA_outlr(:,i,2) == 1
%         continue
%     end
%     if sum(A) > 1
%         figure();
%         plot(TIME(40:end), PARAMETERS_CENT(idx).CA_cent(40:end,i,2));
%         xlim([0 2])
%         ylim([-400 800])
%         title(sprintf("Rapid change in Col %d and target %d \n", i,2))
%         fprintf("Rapid change in Col %d and target %d \n", i,2);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






