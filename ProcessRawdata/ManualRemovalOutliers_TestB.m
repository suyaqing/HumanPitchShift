% Author: Philipp Eugster
% Created: 21.12.2022

%%% NEEDS
% FILE: PARAMETERS_CENT.mat


%%% Suddend and unnatural jumps in pitch, wich almost only happen at the very end of a trial are removed manually here. Used same criteria as descripted in the follwoing paper:
%https://www.biorxiv.org/content/10.1101/2020.06.06.138263v2.full.pdf (Title: Spontaneous variability predicts adaptive motor response in vocal pitch control, Ryosuke O. Tachibana)






%%% Test B
idx = find("f17" == sbjnames);
PARAMETERS_CENT(idx).B_HZ(104:108,12,1) = nan;

idx = find("f21NC" == sbjnames);
PARAMETERS_CENT(idx).B_HZ(94:96,3,2) = 203.0325;
idx = find("f22NC" == sbjnames);
PARAMETERS_CENT(idx).B_HZ(167:end,4,1)  = nan;
PARAMETERS_CENT(idx).B_HZ(177:end,5,1) = nan;

idx = find("f24" == sbjnames);
PARAMETERS_CENT(idx).B_HZ(125:end,12,2) = nan;

idx = find("m02NC" == sbjnames);
PARAMETERS_CENT(idx).B_HZ(150:160,16,2) = nan;

