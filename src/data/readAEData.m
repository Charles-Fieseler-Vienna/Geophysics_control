function [AE_Data, varargout ] = readAEData(fistr,interactive_select, tw_type)
%% Function to Import AE data from csv or text files
%% Parameters 
%  Example: 
%  Default -  [AE_Data ] = readAEData() 
%  Options - [AE_Data, s_int, ch, hit, time, np, time_ax ] = readAEData()
% Inputs: select folder with AE data
%                 data_cond: y for file with headers / n for file without
%                 headers
% Possible Outputs: 
% 1                    AE_Data - Array of AE waveforms, the columns represent a new waverform
% 2                    samp_int - sampling time 
% 3                    channel_data - channel number extracted from the waveform data
% 4                    hit_num_data - hit number extracted from the waveform data, correlates with the time the event occurred, i.e. 1 beign early - 10 being late in testing. 
% 5                    test_time_data - time the AE event/hit occurred could be in seconds (format: ss.ssss) or date and time (format: dd:hh:mm:ss.ssss) there is an option if this occurs
% 6                    np - number of data points in a waveform
% 7                    t_ax - time axis
% -----------------------------------------------------------
% When interactive select is activated the outputs are shifteds as follows:
% Possible outputs: 
%              1 = AE_Data
%             2 = samp_int;
%             3 = channel_data;
%             4 = hit_num_data;
%             5 = t_wind1;
%             6 = t_wind2;
%             7 = test_time_data;
%             8 = np;
%             9 = time_axis;
%         
%%
% Created by: Chven Mitchell 
% Contact: mitch240@purdue.edu | chven.mitchell@gmail.com
% Institution: Purdue University
% Last Updated: 2020/0420
% Advisor: Laura Pyrak-Nolte (ljpn@purdue.edu)
%% Last Updates Made
% Converted to function to accept both text or csv files for learning control signals  toolbox
% Updateds to Make
% Output text file should have header rows 
%% 
% clear
% clear all
clc

% Check version/ Matlab Year
ver_i = ver('Matlab');
ver_rel=strsplit(ver_i.Release,[")","(","a","R"]);
ver_rel = str2double(ver_rel{1,2});

if ver_rel >= 2019
    c_opt =1;
else
    c_opt = 2;
end
%%
if nargin==1
    interactive_select = 'n';
    tw_type = 'none';
end 
if nargin==2
    tw_type = 'none';
end 
%%
% Get the Directory and Files in the folders
% For mac weirdly: users will need to select csv or txt files using the option button on the bottom left of the file load dialog window or if you
% only work with text files remove *.csv or vice versa
disp('Select AE Files for analysis (*.csv | *.txt) : ')
% [baseName, dirName]=uigetfile({'*.txt';'*.csv'},'Select the INPUT DATA FILE(s)','MultiSelect','on');
% dirName = uigetdir();
dirName =  '/Volumes/vermicullite/CAM_Data/Res_AcousticEmissions/20201108/DTA/20-075_MOR/AEdata_Output/20-075_MOR/';
baseName = dir(fullfile(dirName, '*.csv'));
% disp(extractfield(baseName, 'name'))
baseName = {baseName.name}';

% [numfiles, ~] = size(baseName);
% baseName = [dirName, '/', baseName.name];
% Handling for if one file is selected
% if ischar(baseName)==1 || numfiles ==1
%     baseName = {baseName};
% else
%     baseName = baseName';
% end

% disp(numfiles)
%%
% Create a new folder in the top folder where the data is located to place
% all outputs from the script
new_folder_name = strcat(fistr, '_LControl_Signal_Outputs');
source_dir = dirName;
ParentFolder = source_dir;  % Defines the user selection as the parent folder i.e choose the data directory
NewSubFolder = fullfile(source_dir,new_folder_name); %defines location of the new folder
   if ~exist(NewSubFolder, 'dir')
       mkdir(ParentFolder,new_folder_name);
   end
disp('Moving into output folder')
cd(NewSubFolder)
%%
% Number of files loaded
numfiles = length(baseName);
f_chk_id = 1;           % File check ID for while loop
%% If interactive select is on and the window is fixed get window length from users
if strcmp(interactive_select,'y')==1 && tw_type=='f'
        % Get user input for the fixed time window
        fiprompt = {'Please enter fixed window length in seconds: '};
        dlgtitle = 'Fixed Window Input';
        dims = [1 50];
        definput = {'0.00005'};
        answer = inputdlg(fiprompt,dlgtitle,dims,definput);
        fixed_window = str2double(answer{1});
end 
%% Setup the Import Options
delimiter = ',\t';
startRow = 1;
formatSpec = '%s%s%s';
%% Checking data to see if it is two colum [ time, AE data] or one column [AE Data] 
filename = [dirName,'/',baseName{1,1}];
disp(filename)
% Open the text file.
fileID = fopen(filename,'r');

textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false);

% Close the text file.
fclose(fileID);

% Extract data from files
extract_Data = dataArray{1,1};
data_col2 = str2double(dataArray{1,2});
% disp(length(extract_Data));

% Checking the data for Nans. Assuming that the presence of Nans relate to
% the presence of header. Otherwise the file has no headers
index=isnan(str2double(extract_Data));
if sum(index)==0
    disp('This file has no headers')
    data_cond = 'n';
     d_line_id = 1;         % Start of AE Data in file
else
    data_cond ='y';
    d_line_id = 12;         % Start of AE Data in file
end
%%
% If the file has headers get this information
        if strcmp(data_cond,'y')==1
            % Get AE Testing information
            % Extract the sample interval in seconds: Can be used to build the time axis
            samp_int = strsplit(extract_Data(4,1),':');
            samp_int = str2double(samp_int(1,2));

            % Extract the number of data points for the dataset
            np = strsplit(extract_Data(8,1),':');
            np = str2double(np(1,2));
        else
            % get sampling dt from user
            fiprompt = {'Please provide the dt in (secs): '};
            dlgtitle = 'User Input';
            dims = [1 50];       
            definput = {'5.0e-7'};
            answer = inputdlg(fiprompt,dlgtitle,dims,definput);
            samp_int = str2double(answer{1});
            % samp_int = input('Please provide the dt in (secs): ');
           [ np,~]=size(extract_Data);
        end
 
        % Preallocate for AE output with headers + Data
        AE_Data_output = zeros(np+3, numfiles);      % for writing out text files
        
        % Preallocate for AE Data
        AE_Data = zeros(np, numfiles);
        AE_Time_Data = zeros(np, numfiles);
        
            if isnan(data_col2(d_line_id,1))==1
                ae_ax_id = 1;
                t_ax_id = 0; 
                time_axis = [0:np-1].*samp_int;
            else
                t_ax_id = 1;     % Time axis column id
                ae_ax_id = 2;        % AE axis column id
                
                % Extract Time Array
                % time_axis = str2double(extract_Data(d_line_id:end,t_ax_id));
                % time_axis = time_axis + abs(time_axis(1,1));
            end

             %% Preallocate Arrays and Extract Data
                channel_data = zeros(numfiles,1);
                hit_num_data = zeros(numfiles, 1);
                test_time_data = zeros(numfiles,1);
%                             
%                 if c_opt==1
%                     writematrix(time_axis, append(fistr,'_AE_Data_Time_Axis.txt'),'Delimiter','tab');
%                 else
%                     dlmwrite(strcat(fistr,'_AE_Data_Time_Axis.txt'),time_axis,'precision','%.9f', 'delimiter', '\t' );
%                 end
            if strcmp(interactive_select,'y')~=1
%                 prompt = {'Enter Start Time:','Enter End Time:'};
%                 dlgtitle = 'Enter Fixed Time Windows for Data';
%                 dims = [1 35];
%                 definput = {' 1.500e-05','8.000e-04'};
%                 answer_twind = inputdlg(prompt,dlgtitle,dims,definput);
%                 t_wind1 = str2double(answer_twind{1});
%                 t_wind2 = str2double(answer_twind{2});
                     t_wind1 =  '4.500e-05';
                     t_wind2 = '1.500e-04';
            end
%%  Loop over remaining files and extract the data
while f_chk_id <=(numfiles )
    filename = [dirName,'/', baseName{f_chk_id,1}];
    fileID = fopen(filename,'r');
    %  disp(fileID)
    textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false);
    fclose(fileID);
    
    if ae_ax_id==2
        extract_Data = dataArray{1, 2};
        wav_data = dataArray{1, 1};
         time_axis = str2double(wav_data(d_line_id:end));
         time_axis = time_axis + abs(time_axis(1,1));
    else
        extract_Data = dataArray{1,1};
        wav_data = dataArray{1, 1};
    end
  %% 
  
    if strcmp(data_cond,'y')==1
        % Get the Channel number per waveform
        ch_num = strsplit(wav_data(9,1),':');
        ch_num = str2double(ch_num(1,2));
        
        % Get the Hit number per waveform
        
        hit_num = strsplit(wav_data(10,1),':');
        hit_num = str2double(hit_num(1,2));
   
        
        % Time of Test in (seconds)
        % Note: This might be specific to my data be care if your output is in dd:hh:mm:ss.sssss
        ttime_num = strsplit(wav_data(11,1),':');
        tlen = length(ttime_num);
        
        if tlen>2
            ttime_num = join([strtrim(ttime_num(2)),ttime_num(3:end)]); % Join the strings, first removing the leading space in front of the first number
            ttime_num = strrep(ttime_num,' ',':');          % replace the white space with : for time
        else
            ttime_num = str2double(ttime_num(1,2));    % Else use this for time output in seconds (ss.ssssss) from start of testing i.e. 12.678559452
        end
        
         if f_chk_id > 960
            ttime_num = ttime_num + 84777.802153;
        end
        
        % Store data to arays
        channel_data(f_chk_id,1) = ch_num;
        hit_num_data(f_chk_id,1) = hit_num;
        test_time_data(f_chk_id,1) = ttime_num;
    else
        % Use the file name to build the channel, hit and time information
        filename_parts = split(baseName{f_chk_id,1},["_",".","-"]);
        f_len = length(filename_parts);
        % Get channel from file name
        f_ch = str2double(filename_parts{2,1});
        if isnan(f_ch)
            f_ch = 0;
        end 
        
        % Get the hit number from the file name and check to be sure its
        % numeric not string
        if f_len>2 
            f_hit = str2double(filename_parts{end-2,1}); 
            if isnan(f_hit)
                f_hit=0;
            end
            f_time = 0;
        end
        if f_len>3
            f_time = str2double(filename_parts{end-1,1}); 
            if isnan(f_time)
                f_time = 0;
            end
        end
        channel_data(f_chk_id,1) = f_ch;
        hit_num_data(f_chk_id,1) = f_hit;
        test_time_data(f_chk_id,1) = f_time;
    end  
    
    AE_Data(1:np,f_chk_id) = str2double(extract_Data(d_line_id:end,1));    
    AE_Time_Data(1:np,f_chk_id) = str2double(time_axis);    
%% -------------------------------------------------------------------------------------------------
%      Interactively selecting t_wind1 and twind_2
%      interactive_select = input('Do you wish to select the time windows (y/n): ','s');
if strcmp(interactive_select,'y')==1 && tw_type=='v'
    disp('Click to SET (start) t_wind1 and (end) t_wind2 time windows: ')
    plot(time_axis, AE_Data(:,f_chk_id));
    %          xlim([0 (np/2)]*samp_int);
    % Pick one set of two-dimensional
    %                 set(gcf, 'Position', get(0, 'Screensize'));
    enableDefaultInteractivity(gca);
    [pk_x, pk_y] =ginput(2);
    t_wind1(f_chk_id,1) = pk_x(1,1);
    t_wind2(f_chk_id,1) = pk_x(2,1);
    
else if strcmp(interactive_select,'y')==1 && tw_type=='f'
        disp('Click to SET (start) t_wind1: ')
        plot(time_axis, AE_Data(:,f_chk_id));
        %  set(gcf, 'Position', get(0, 'Screensize'));  % To  make the waveform full the screen
        enableDefaultInteractivity(gca);
        [pk_x, ~] =ginput(1);
        t_wind1(f_chk_id,1) = pk_x(1,1);
        t_wind2 = t_wind1 + fixed_window;
    end
end
    f_chk_id =f_chk_id + 1;
end
%% Set varargout Outputs
    if  nargout==2
        varargout{1} = samp_int;   
    end

        if  nargout==3
            varargout{1} = samp_int;
            varargout{2} = channel_data;
        end
        
        if  nargout==4
            varargout{1} = samp_int;
            varargout{2} = channel_data;
            varargout{3} = hit_num_data;
        end
        
        if  nargout==5
            varargout{1} = samp_int;
            varargout{2} = channel_data;
            varargout{3} = hit_num_data;
            varargout{4} = baseName;
        end
        
        
    if interactive_select=='y'
        if  nargout==7
            varargout{1} = samp_int;
            varargout{2} = channel_data;
            varargout{3} = hit_num_data;
            varargout{4} = baseName;            
            varargout{5} = t_wind1;
            varargout{6} = t_wind2;
        end

        if  nargout==8
            varargout{1} = samp_int;
            varargout{2} = channel_data;
            varargout{3} = hit_num_data;
            varargout{4} = baseName;            
            varargout{5} = t_wind1;
            varargout{6} = t_wind2;
            varargout{7} = test_time_data;
        end

        if  nargout==9
            varargout{1} = samp_int;
            varargout{2} = channel_data;
            varargout{3} = hit_num_data;
            varargout{4} = baseName;                   
            varargout{5} = t_wind1;
            varargout{6} = t_wind2;
            varargout{7} = test_time_data;
            varargout{8} = np;
        end
        
         if  nargout==10
            varargout{1} = samp_int;
            varargout{2} = channel_data;
            varargout{3} = hit_num_data;
            varargout{4} = baseName;                   
            varargout{5} = t_wind1;
            varargout{6} = t_wind2;
            varargout{7} = test_time_data;
            varargout{8} = np;
            varargout{9} = time_axis;
        end
        
    else
        if  nargout==6
            varargout{1} = samp_int;
            varargout{2} = channel_data;
            varargout{3} = hit_num_data;
            varargout{4} = baseName;                   
            varargout{5} = test_time_data;
        end

        if  nargout==7
            varargout{1} = samp_int;
            varargout{2} = channel_data;
            varargout{3} = hit_num_data;
            varargout{4} = baseName;                   
            varargout{5} = test_time_data;
            varargout{6} = np;
        end

        if  nargout==8
            varargout{1} = samp_int;
            varargout{2} = channel_data;
            varargout{3} = hit_num_data;
            varargout{4} = baseName;                   
            varargout{5} = test_time_data;
            varargout{6} = np;
            varargout{7} = time_axis;
        end

        if  nargout==10
            varargout{1} = samp_int;
            varargout{2} = channel_data;
            varargout{3} = hit_num_data;
            varargout{4} = baseName;                   
            varargout{5} = test_time_data;
            varargout{6} = np;
            varargout{7} = time_axis;
            varargout{8} = t_wind1;
            varargout{9} = t_wind2;
        end
    end
%%  Build output file with all informaiton 
        AE_Data_output(1,:) = channel_data';
        AE_Data_output(2,:) = hit_num_data';
        AE_Data_output(3,:) = test_time_data';
        AE_Data_output(4:end,:) = AE_Data;
        
        new_wind = ones(numfiles,1); 
        t_wind1 = new_wind*t_wind1;
        t_wind2 = new_wind*t_wind2;

        % Write data to text files
        if c_opt ==1
            writematrix(channel_data, append(fistr,'_AE_channel_Data.txt'),'Delimiter','tab');
            writematrix(hit_num_data, append(fistr, '_AE_HitNumber_Data.txt'),'Delimiter','tab');
            writematrix(test_time_data, append(fistr,'_AE_hit_Time_Data.txt'),'Delimiter','tab');
            writematrix([channel_data, t_wind1, t_wind2], append(fistr,'_AE_picked_Times_Data.txt'),'Delimiter','tab');
            writematrix(AE_Data_output, append(fistr,'_AE_Data_ch_hit_testime.txt'),'Delimiter','tab');
        else
            dlmwrite(strcat(fistr,'_AE_channel_Data.txt'),channel_data,'precision','%.0f', 'delimiter', '\t');
            dlmwrite(strcat(fistr,'_AE_HitNumber_Data.txt'),hit_num_data,'precision','%.0f', 'delimiter', '\t' );
            dlmwrite(strcat(fistr,'_AE_hit_Time_Data.txt'),test_time_data,'precision','%.9f', 'delimiter', '\t');
            dlmwrite(strcat(fistr,'_AE_picked_Times_Data.txt'),[channel_data, t_wind1, t_wind2],'precision','%.9f', 'delimiter', '\t');
            dlmwrite(strcat(fistr,'_AE_Data_ch_hit_testime.txt'),AE_Data_output,'precision','%.12f', 'delimiter', '\t' );
        end

%% Clear temporary variables

clear fileID dataArray
close all
% clearvars -except AE_Data
end
