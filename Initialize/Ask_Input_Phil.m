%CHECK IF YOU WANT TO READ IN RAWDATA OR DIRECT
quest = 'Read in from Rawdata or from PARAMETERS.m file?';
dlgtitle = 'Read in from Rawdata?';
btn1 = 'Raw data';
btn2 ='PARAMETERS file';

answer = questdlg(quest, dlgtitle, btn1, btn2, btn1);
% Handle response
switch answer
    case btn1
        READ_RAWDATA = true;
    case btn2
        READ_RAWDATA = false;
end

if READ_RAWDATA == false
    datapath = uigetdir("C:\Users\Philipp\Desktop\Semesterproject","Select Folder where Parameters file is stored");
    datapath = [datapath '\' ];
    filename = [datapath '\PARAMETERS_CENT' '.mat'];
    load(filename);
    quest = 'Plot rawdata & save plots for each participant?';
    dlgtitle = 'Plot rawdata?';
    btn1 = 'Yes';
    btn2 ='No';
    answer = questdlg(quest, dlgtitle, btn1, btn2, btn2);
    % Handle response
    switch answer
        case btn1
            PRODUCE_PLOTS = true;
        case btn2
            PRODUCE_PLOTS = false;
    end
end

if READ_RAWDATA == true
    datapath = uigetdir("C:\Users\Philipp\Desktop\Semesterproject","Select Folder where Rawdata is stored");
    datapath = [datapath '\'];
    quest = 'Plot rawdata & save plots for each participant?';
    dlgtitle = 'Plot rawdata?';
    btn1 = 'Yes';
    btn2 ='No';
    answer = questdlg(quest, dlgtitle, btn1, btn2, btn2);
    % Handle response
    switch answer
        case btn1
            PRODUCE_PLOTS = true;
        case btn2
            PRODUCE_PLOTS = false;
    end
end


savepath = uigetdir('C:\Users\Philipp\Desktop\Semesterproject\MYANALYSIS',"Select location where to save Parameters & Plots");
savepath = [savepath '\' ];

fid = fopen( 'getPath.m', 'wt' );
fprintf(fid, 'savepath = ''%s'';\n', savepath);
fprintf(fid, 'datapath = ''%s'';\n', datapath);
fclose(fid);




