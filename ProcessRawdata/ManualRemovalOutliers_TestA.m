% Author: Philipp Eugster
% Created: 21.12.2022

%%% NEEDS
% FILE: PARAMETERS_CENT.mat


%%% Suddend and unnatural jumps in pitch, wich almost only happen at the very end of a trial are removed manually here. Used same criteria as descripted in the follwoing paper:
%https://www.biorxiv.org/content/10.1101/2020.06.06.138263v2.full.pdf (Title: Spontaneous variability predicts adaptive motor response in vocal pitch control, Ryosuke O. Tachibana)



%%% Test A
idx = find("f14NC" == sbjnames);
PARAMETERS_CENT(idx).A_HZ(143:145,12) = nan;

idx = find("m04" == sbjnames);
PARAMETERS_CENT(idx).A_HZ(139:141,9) = nan;

idx = find("m07" == sbjnames);
PARAMETERS_CENT(idx).A_HZ(131:135,9) = nan;
PARAMETERS_CENT(idx).A_HZ(:,10) = nan;
