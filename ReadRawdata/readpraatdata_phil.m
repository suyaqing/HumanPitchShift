function [out,trace,ons] = readpraatdata_phil(filename,reffreq,REC_DURATION, TIMESTEPS_PER_SECOND)



% Usage: out = readpraatdata(filename,reffreq,REC_DURATION,TIMESTEPS_PER_SECOND)
%
%   reffreq ... for correcting octave confusion
%   REC_DURATION &
%   TIMESTEPS_PER_SECOND ... To make function work, if the pitch data is extracted using PRAAT with different sample rate.
%
% Author:  R.Tachibana 2018/03/08, adapted by Philipp Eugster on 12.10.2022 <eugsteph@student.ethz.ch>


% Adapted by PHILIPP EUGSTER, 12.10.2022
% Changes: The function now takes the two additional arguments
% "REC_DURATION" and "TIMESTEPS_PER_SECOND". With this change I make sure that if
% those parameters are changed in the main file, the function works with the
% right parameters. Befor those changes the parameters were hardcoded in
% both files (main analysis file and the function itself) which is prone for
% errors. Old hardcoded variable "duration" = REC_DURATION and old hardcoded
% variable "fs" = TIMESTEPS_PER_SECOND.


len = round(TIMESTEPS_PER_SECOND*REC_DURATION);

% read from file
fid = fopen(filename,'r'); % filename is outpath from praat_script
C = textscan(fid,'%s\t%f\t%f','headerlines',1);
fclose(fid);
ts = C{2};%time vector
f0 = C{3};%f0 traces
if isempty(find(ts>0, 1))
    disp(filename)
    ts = ones(400,1);
    f0 = ones(400,1);
end
trace = zeros(len+TIMESTEPS_PER_SECOND,1);
trace(1:size(ts,1)) = f0;


% remove short fragments and aligned by onset
mindur = 10; % remove fragments less than 100 ms

 % The following line
 % find positions where the trace increases in pitch
ons = find(diff([0;trace]>0)>0); %diff: differences between adjacent elements of X;
offs = find(diff([trace;0]>0)<0)+1;
dur = offs-ons;

% set small fragments to zeros
ngid = find(dur<mindur);
if ~isempty(ngid)
    disp(filename)
    for n=1:length(ngid)
        trace(ons(ngid(n)):offs(ngid(n))) = 0;
    end
end

ons = ons(dur>=mindur);
if isempty(ons)
    out = zeros(len,1);
else
    align = zeros(len,1);
    temp = trace(ons(1):end); % align onset to the beginning of the trace
    align(1:min(length(temp),len)) = temp(1:min(length(temp),len));
    align(align==0) = nan;
    
    % correct octave confusion
    out = align;
    if mean(align,'omitnan') < reffreq * 0.6
        out = align * 2;
    elseif mean(align, 'omitnan') > reffreq * 1.8
        out = align / 2;
    end
end
