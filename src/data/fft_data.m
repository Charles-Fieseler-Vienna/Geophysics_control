%% Moving window DFT to get frequencies
%---------------------------------------------
%% Analysis settings

% Which input files
% which_dataset = 'mortar_fnames';
% which_dataset = 'localized_fnames';
% which_dataset = 'distributed_fnames';
which_dataset = 'unknown_fnames2';
num_files = 250; % 250;
all_dft_len = 63;% 30;

% Intermediate results (output)
intermediate_fname = sprintf('spectrogram_%s.mat', which_dataset);

% Class with preprocessing functions
data_foldername = '';
pp = PurdueProject(data_foldername);
%% Get filenames
%---------------------------------------------
fnames = pp.(which_dataset);
if length(fnames) < num_files
    num_files = length(fnames);
    fprintf("Only found %d files\n", num_files);
end
if strcmp(which_dataset, 'localized_fnames')
    [all_dat, kept_ind] = pp.filter_by_activity(fnames);
    kept_ind = find(kept_ind); % Note that the mortar data is much cleaner
else
    kept_ind = 1:length(fnames);
end
%---------------------------------------------
%% Get spectrogram
%---------------------------------------------
% num_files = length(kept_ind);
all_dft = zeros(1024, all_dft_len, num_files);
all_dat = cell(num_files, 1);

for i = 1:num_files
    fprintf('File %d/%d\n', i, num_files);
    dat_with_time = readtable(fnames{kept_ind(i)});
    % Remove rise, and remove time column
    ind = dat_with_time{:,1} > 0;
    dat = dat_with_time{ind,2}';
    
    all_dat{i} = dat;
    
    t = dat_with_time{end-length(dat)+1:end,1};
    dat_table = timetable(dat', 'RowTimes', seconds(t));
    all_dft(:,:,i) = pspectrum(dat_table, 'spectrogram');
%     all_dft(:,:,i) = pspectrum(dat, ...
%         'spectrogram', 'FrequencyLimits',[0,1.5]);
end
[~, f] = pspectrum(dat_table, 'spectrogram');
%---------------------------------------------
%% Get and plot initial peaks
%---------------------------------------------

all_peaks = zeros(length(num_files));
all_locs = zeros(length(num_files));

for i = 1:num_files
    fprintf('File %d/%d\n', i, num_files);
    % Only look at the first window
    tmp = all_dft(:,1,i);
    [all_peaks(i), all_locs(i)] = max(tmp);
end

% Remove outliers, especially peaks at 0 frequency
ind_to_keep = ~isoutlier(all_locs);
all_locs_clean = all_locs(ind_to_keep);
%% Plot
figure;
histogram(f(all_locs))%, 'BinWidth',0.025)
% histogram(f(all_locs_clean), 'BinWidth',0.025)
xlabel('Frequency')
title(sprintf('Frequency Peaks in initial window for dataset %s',which_dataset))

figure;
scatter(rand(size(all_locs_clean)),f(all_locs_clean))
ylabel('Frequency')
xlabel('Visualization variable')
title(sprintf('Non-outlier Frequency Peaks in initial window for dataset %s',which_dataset))

%% Look at traces with frequencies
to_look_at_all_figures = false;
if to_look_at_all_figures
    for i = 1:num_files
        figure;

        plot(all_dat{i})
        title(sprintf("Time series %d, peak at: %.2f", i, f(all_locs_clean(i))))

        pause
        close all
    end
end
%---------------------------------------------
%% SAVE INTERMEDIATES
%---------------------------------------------
fname = [pp.intermediate_foldername, intermediate_fname];

% save(fname, 'kept_ind', 'all_locs', 'all_peaks', 'all_dft', '-v7.3');
% save(fname, 'f', 'all_dat', 'kept_ind', 'all_locs', 'all_peaks', 't_of_dat', '-v7.3');
save(fname, 'f', 'all_dat', 'kept_ind', 'all_locs', 'all_peaks', '-v7.3');


disp('Intermediate data saved')
% error('The following code is for visualization; only run if interested')

%---------------------------------------------
