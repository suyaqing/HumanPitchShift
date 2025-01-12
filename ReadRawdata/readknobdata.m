%%% This Function reads in the data for the knob turning experiments in Test D

% Author: Dimitri Spicher
% Adapted by: Philipp Eugster <eugsteph@student.ethz.ch>
% Final changes: 07.12.2022


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [time,data,out] = readknobdata(filename)


fs = 100;
duration = 4;
l = load(filename);
time = l.trace(:,1);
data = l.trace(:,2);
data(1) = 0;
% remove overlapped time points
okid = find(diff(time)~=0);
time = [time(okid);time(end)];
data = [data(okid);data(end)];
% linear interporation
xi = (1:fs*duration)'/fs;
out = interp1(time,data,xi); 
