

%% Get largest DMD modes (time-delay embed)
%---------------------------------------------
%% Get filenames
%---------------------------------------------
pp = PurdueProject();
% which_dataset = 'mortar_fnames';
which_dataset = 'localized_fnames';
% which_dataset = 'distributed_fnames';
fnames = pp.(which_dataset);
if strcmp(which_dataset, 'localized_fnames')
    [all_dat, kept_ind] = pp.filter_by_activity(fnames);
    kept_ind = find(kept_ind); % Note that the mortar data is much cleaner
else
    kept_ind = 1:length(fnames);
end

%% DMD (takes a while)
num_files = length(kept_ind);
aug = 50;

num_modes = 6;
all_coef = zeros(num_modes, num_files);
all_eigs = zeros(num_modes, num_files);
all_phi = zeros(aug, num_modes, num_files);


for i = 1:num_files
% for i = 1:10
    fprintf('File %d/%d\n', i, num_files);
    dat = readtable(fnames{kept_ind(i)});
    % Remove rise time
    ind = dat{:,1} > 0;
    dat = dat{ind,2}';
    dat_embed = time_delay_embed(dat, aug);

    % Strongly truncated SVD-based DMD
    [all_coef(:,i), all_eigs(:,i), all_phi(:,:,i)] = ...
        dmd(dat_embed, [], num_modes);
end

%% Use AIC to see how many modes to keep
all_aic = zeros(num_modes/2, num_files);
all_aic_mins = zeros(1, num_files);

for i = 1:num_files
    fprintf('File %d/%d\n', i, num_files);
    dat = readtable(fnames{kept_ind(i)});
    % Remove rise time
    ind = dat{:,1} > 0;
    dat = dat{ind,2}';
%     figure;plot(dat)
    dat_embed = time_delay_embed(dat, aug);

    % Calculate reconstruction to see how many modes to keep
    for i2 = 1:(num_modes/2)
        D = all_eigs(1:2*i2,i);
        V = all_phi(:,1:2*i2,i);
        A = V*diag(D)/V;
        all_aic(i2, i) = aic_dmd(dat_embed, A, i2, 100);
    end
    [~, all_aic_mins(i)] = min(all_aic(:,i));
end

%% SAVE INTERMEDIATE
% fname = pp.intermediate_foldername+"test.mat";
% fname = pp.intermediate_foldername+"distributed_dmd.mat";
% fname = pp.intermediate_foldername+"mortar_dmd.mat";
fname = pp.intermediate_foldername+"localized_dmd_50aug.mat";

% save(fname, 'all_coef', 'all_eigs', 'all_phi', 'all_aic', 'all_aic_mins',...
%     '-v7.3');
save(fname, 'all_coef', 'all_eigs', 'aug', 'num_modes',...
    'all_aic', 'all_aic_mins',...
    '-v7.3');

%% Plot maxes
figure
plot(all_aic_mins, 'o')
title("Number of modes for each dataset")

% figure;
% plot(all_aic(:,11))

%% Calculate some features to use as colors or sizes below
h = zeros(1,num_files);
% size_func = @(x) max(abs(x(1,:)));
color_func = @(x) max(svd(x))/sum(svd(x));
for i = 1:num_files
    fprintf('File %d/%d\n', i, num_files);
    dat = readtable(fnames{kept_ind(i)});
    % Remove rise
    ind = dat{:,1} > 0;
    dat = dat{ind,2}';
    dat_embed = time_delay_embed(dat, aug);
    
    h(i) = size_func(dat_embed);
end

h2 = h;

%% VERSION 2: USE COEFS
h = abs(all_coef(1,:));
h2 = abs(all_coef(3,:));
%% VERSION 3: Use metadata!
% mdata = pp.import_metadata_file('mortar');
crack_times = zeros(1, num_files);

for i = 1:num_files
    this_file = fnames{kept_ind(i)};
%     m = pp.fname2metadata(this_file, mdata);
    crack_times(i) = pp.read_time_of_test(this_file);
end

h = crack_times;
% Remove the extremely late events (days later)
% thresh = median(h);
% h(h>thresh) = 1e-5;
h = log(1 + h);

%% Process imaginary eigenvalues to get primary and secondary

all_e_real = [real(all_eigs(1,:)); imag(all_eigs(1,:));...
    real(all_eigs(3,:)); imag(all_eigs(3,:))]';

%% Plot all frequencies
figure;
% scatter(all_e_real(:,1),all_e_real(:,2), 100*(h./max(h)), all_aic_mins)
% scatter(all_e_real(:,3),all_e_real(:,4), 100*(h./max(h)), all_aic_mins)
s = scatter(all_e_real(:,1),all_e_real(:,2), 100*(h./max(h)), all_aic_mins);
% colormap(hsv)
colorbar
% hold on
% scatter(all_e_real(:,3),all_e_real(:,4), 100*(h2./max(h)))
title(sprintf('Dataset %s', which_dataset))

%% Simple GUI for looking at the data that produced the point
x = ginput(1);
% all_e_real = [real(all_e(1,:)); imag(all_e(1,:))]';
idx = knnsearch(all_e_real(:,1:2), x);

% Import
dat = readtable(fnames{kept_ind(idx)});
% Remove rise
ind = dat{:,1} > 0;
dat = dat{ind,2}';
figure;
plot(dat)
title(sprintf('File %d; eigenvalue %.4f,%.4f', kept_ind(idx), x(1), x(2)))


%% Make a video
% vidObj = VideoWriter('../figures/time_of_test_mortar.avi');
vidObj = VideoWriter('../figures/time_of_test_localized.avi');
open(vidObj);

figure;
hold on
title(sprintf('Dataset %s', which_dataset))
h_norm = 100*(h./max(h));
[~, ind] = sort(h_norm);
for i = 1:length(h)
    i2 = ind(i);
    plot(all_e_real(i2,1),all_e_real(i2,2),'bo')
%     drawnow
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end
close(vidObj);
%==========================================================================


