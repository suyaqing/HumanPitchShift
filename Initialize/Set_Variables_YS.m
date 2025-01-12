
%%% In this file, parameters are set that should not be changed since they are specific to the 
%%% Experiments performed by Mengli. Saves all the set parameters in the file variables.mat at the specified savepath (specified in Main.m).


% Author: Philipp Eugster <eugsteph@student.ethz.ch>
% Created: 18.12.2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%% DONT CHANGE!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REC_DURATION = 2; %duration of voice recording in seconds for test A, B & C
REC_DURATION_D = 4; %duration of voice recording in seconds for test D

%LENGTH = REC_DURATION * TIMESTEPS_PER_SECOND; % Number of possible data points for a recording, since it includes ZERO, it goes from 0 to REC-DURATION - 1x (1/TIMESTEPS_PER_SECOND).
%(2 seconds with a intervall size of 0.01 seconds == 200 datapoints).

TIME = (-1:1/TIMESTEPS_PER_SECOND:REC_DURATION-(1/TIMESTEPS_PER_SECOND))'; % TIME stores all the time intervalls, used for time-axis in plots.
TIME_D = (0:1/TIMESTEPS_PER_SECOND:REC_DURATION_D-(1/TIMESTEPS_PER_SECOND))'; % TIME stores all the time intervalls, used for time-axis in plots.
LENGTH = length(TIME);
TARGET = [0 400]; % The two targets used in the experiment per participant in cent.
SHIFT = [-100 -50 -25 0 25 50 100]; % The different shifts in cent that are used in Test CB and also in D.

SESSIONS_C = 5; % Number of sessions in task C in Experiment perfomred by Mengli
SESSIONS_D = 5; % Number of sessions in task D in Experiment perfomred by Mengli

FS = 44100;
%%%%%%%%%%%%%%%% Read in reference frequency, do not change %%%%%%%%%%%%%%%

%%% Reads in the reference frequency from the rawdata (if READ_RAWDATA is set to "true" in Main.m) or from the PARAMETERS_CENT.mat file if set to "false".

if READ_RAWDATA == true
    for b=1:nsubj
        a = load([datapath sbjnames{b} filesep 'ref_freq.mat']);
        ref_freq(b) = a.ref_freq;
        if ref_freq(b) < 1
            disp("ERROR, NO Reference frequency found for:");
            disp(ref_freq(b));
        end
    end

elseif READ_RAWDATA == false
    for b=1:nsubj
        ref_freq(b) = PARAMETERS_CENT(b).ref_freq;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = [savepath 'variables' '.mat'];
delete([savepath 'variables' '.mat']);
save(filename,'REC_DURATION','TIMESTEPS_PER_SECOND','RESOLUTION','FORMAT','SINGLE_PLOTS','CONTENT_TYPE','WINDOW_OUTLIER_1','WINDOW_OUTLIER_2','WINDOW_OUTLIER_3','WINDOW_OUTLIER_4','WINDOW_OUTLIER_5', 'WINDOW_OUTLIER_D','LENGTH','TIME','nsubj','sbjnames',"SHIFT","TARGET","TIME_D","REC_DURATION_D","SESSIONS_D","ref_freq","SESSIONS_C","READ_RAWDATA","PRODUCE_PLOTS");


disp('ALL PARAMETERS SET & SAVED!')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%