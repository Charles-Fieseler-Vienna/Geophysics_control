%% Moving window DFT to get frequencies
%---------------------------------------------
%% Analysis settings

% Which input files
% which_dataset = 'mortar_fnames';
which_dataset = 'localized_fnames';
% which_dataset = 'distributed_fnames';
num_files = 10;

% Intermediate results
intermediate_fname = 'test.mat';

% Class with preprocessing functions
data_foldername = '';
pp = PurdueProject(data_foldername);

%% Get filenames
%---------------------------------------------
fnames = pp.(which_dataset);
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
all_dft = zeros(1024, 30, num_files);

for i = 1:num_files
    fprintf('File %d/%d\n', i, num_files);
    dat = readtable(fnames{kept_ind(i)});
    % Remove rise
    ind = dat{:,1} > 0;
    dat = dat{ind,2}';
    
    all_dft(:,:,i) = pspectrum(dat, ...
        'spectrogram', 'FrequencyLimits',[0,1.5]);
end
[~, f, t] = pspectrum(dat, 'spectrogram', 'FrequencyLimits',[0,1.5]);
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
histogram(f(all_locs), 'BinWidth',0.025)
% histogram(f(all_locs_clean), 'BinWidth',0.025)
xlabel('Frequency')
title(sprintf('Frequency Peaks in initial window for dataset %s',which_dataset))

figure;
scatter(rand(size(all_locs_clean)),f(all_locs_clean))
ylabel('Frequency')
xlabel('Visualization variable')
title(sprintf('Non-outlier Frequency Peaks in initial window for dataset %s',which_dataset))

%---------------------------------------------
%% SAVE INTERMEDIATES
%---------------------------------------------
fname = [pp.intermediate_foldername,intermediate_fname];

% save(fname, 'kept_ind', 'all_locs', 'all_peaks', 'all_dft', '-v7.3');
save(fname, 'f', 'kept_ind', 'all_locs', 'all_peaks', '-v7.3');


disp('Intermediate data saved')
error('The following code is for visualization; only run if interested')

%---------------------------------------------
%% Find recurrence
%---------------------------------------------
mean_location = round(mean(all_locs_clean));
% dat = mean(all_dft((mean_location-5):(mean_location+5),:,:),1);
dat = all_dft(mean_location,:,:);
dat = reshape(dat, size(all_dft,2), size(all_dft,3));
% Normalize the power spectra
all_maxes = max(dat,[],1);
dat = dat ./ all_maxes;

% Plot two different subsets, split by signal amplitude
split_max = median(all_maxes);
ind_large = all_maxes > 200*split_max;
ind_small = ~ind_large;

fig = plot_std_fill(dat(:,ind_large), 2);
title(sprintf('Large Amplitude %s Datasets (%d/%d)', ...
    which_dataset, nnz(ind_large), length(ind_large)))
ylim([0,1])

plot_std_fill(dat(:,ind_small), 2);
title(sprintf('Small Amplitude %s Datasets (%d/%d)', ...
    which_dataset, nnz(ind_small), length(ind_large)))
ylim([0,1])
%==========================================================================




%% 3d surface visualization (spectrogram)
which_dataset = 'mortar_fnames';
% which_dataset = 'localized_fnames';
fnames = pp.(which_dataset);

% i = 11; % Good complicated-ish example
% i = 15; % Unique peak
% i = 21; % Good recurrence
i = 32;
dat = readtable(fnames{kept_ind(i)});
% Remove rise
ind = dat{:,1} > 0;
dat = dat{ind,2}';

title_str = sprintf('Dataset %s; file %d', which_dataset, i);

[p,f,t] = pspectrum(dat, 'spectrogram', 'FrequencyLimits',[0,1.5]);
figure;
surf(f,t,p')
shading interp
xlabel('Frequency')
ylabel('Time (frames)')
zlabel('Power')
title(title_str)

figure;
plot(dat)
title(title_str)
%==========================================================================



%% Spectrogram animation
pp = PurdueProject();
% which_dataset = 'mortar_fnames';
% which_dataset = 'localized_fnames';
which_dataset = 'distributed_fnames';
fnames = pp.(which_dataset);
if strcmp(which_dataset, 'localized_fnames')
    [all_dat, kept_ind] = pp.filter_by_activity(fnames);
    kept_ind = find(kept_ind); % Note that the mortar data is much cleaner
else
    kept_ind = 1:length(fnames);
end
num_files = length(kept_ind);

%% Sort by real time
% Plot DMD values vs. (real) time
crack_times = zeros(1, num_files);
for i = 1:num_files
    fprintf('File %d/%d\n', i, num_files);
    this_file = fnames{kept_ind(i)};
    crack_times(i) = pp.read_time_of_test(this_file);
end

h = crack_times;
% Sort
[sort_times, sort_ind] = sort(crack_times, 'ascend');

kept_ind = kept_ind(sort_ind);

%% Open video object
vidObj = VideoWriter('../figures/spectrogram_localized_raw.avi');
vidObj.FrameRate = 1;
open(vidObj);

% Plot!
for i = 1:num_files
    dat = readtable(fnames{kept_ind(i)});
    % Remove rise
    ind = dat{:,1} > 0;
    dat = dat{ind,2}';
    title_str = sprintf('Dataset %s; file %d', ...
        which_dataset, kept_ind(i));

    [p,f,t] = pspectrum(dat, 'spectrogram', 'FrequencyLimits',[0,1.5]);
    fig = figure;
    subplot(2,1,1)
    surf(f,t,p')
    shading interp
    xlabel('Frequency')
    ylabel('Time (frames)')
    zlabel('Power')
    zlim([0, 2e-5])
    title(title_str)
    
    subplot(2,1,2)
    plot(dat)
    xlim([0,1500])
    ylim([-.1, .1])
    
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    close(fig)
end

close(vidObj);
%==========================================================================

