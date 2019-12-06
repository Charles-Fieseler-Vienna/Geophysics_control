

%% "De-noise" via SVD truncation
pp = PurdueProject();

% i = 9;
% fnames = pp.mortar_fnames;
%i = 8;
i = 12;
fnames = pp.localized_fnames;
dat = readtable(fnames{i});
to_include_rise = false;
if ~to_include_rise
    ind = dat{:,1} > 0;
    dat = dat{ind,2}';
else
    dat = dat{:,2}';
end

aug = 30;
dat_embed = time_delay_embed(dat, aug);
m = size(dat_embed, 2);

%---------------------------------------------
% Preprocess the td-embeded data
%---------------------------------------------
[r, S_vec, dat_signal1, U, S, V] = optimal_truncation(dat_embed);

r_custom = 2;
[dat_signal2] = svd_truncate(dat_embed, r_custom);

% Plot the SVD
plotSVD(dat_embed', struct('sigma_modes',[1,3,5]));
% plotSVD(dat_signal2');

%% Compare
figure;
plot(dat_embed(1,:));
hold on
plot(dat_signal1(1,:))
plot(dat_signal2(1,:))
legend('Data', sprintf('SVD reconstruction %d modes',r), ...
    sprintf('SVD reconstruction %d modes',r_custom))

% For signal processing
d0 = dat_embed(1,:)';
d1 = dat_signal1(1,:)';
d2 = dat_signal2(1,:)';
%==========================================================================



%% 2-d DMD truncation for many datasets

which_dataset = 'mortar_fnames';
% which_dataset = 'localized_fnames';
fnames = pp.(which_dataset);
[all_dat, kept_ind] = pp.filter_by_activity(fnames);
kept_ind = find(kept_ind); % Note that the mortar data is much cleaner
%% (above takes a while)
num_files = length(kept_ind);
all_e = zeros(2, num_files);

aug = 30;

for i = 1:num_files
    fprintf('File %d\n', i);
    dat = readtable(fnames{kept_ind(i)});
    % Remove rise
    ind = dat{:,1} > 0;
    dat = dat{ind,2}';
    dat_embed = time_delay_embed(dat, aug);

    % Strongly truncated Optimized DMD
    r = 2;
    imode = 1;
    t = 1:size(dat_embed, 2);
    [~,all_e(:,i)] = optdmd(dat_embed, t, r, imode);
end

%% Calculate some features to use as colors or sizes below
h = zeros(1,num_files);
% color_func = @(x) max(abs(x(1,:)));
color_func = @(x) max(svd(x))/sum(svd(x));
for i = 1:num_files
    dat = readtable(fnames{kept_ind(i)});
    % Remove rise
    ind = dat{:,1} > 0;
    dat = dat{ind,2}';
    dat_embed = time_delay_embed(dat, aug);
    
    h(i) = color_func(dat_embed);
end

%% Plot all frequencies
figure;
all_e_real = [real(all_e(1,:)); imag(all_e(1,:))]';
scatter(all_e_real(:,1),all_e_real(:,2), 1000*h)
title(sprintf('Dataset %s', which_dataset))

%% Simple GUI for looking at the data that produced the point
x = ginput(1);
% all_e_real = [real(all_e(1,:)); imag(all_e(1,:))]';
idx = knnsearch(all_e_real, x);

% Import
dat = readtable(fnames{kept_ind(idx)});
% Remove rise
ind = dat{:,1} > 0;
dat = dat{ind,2}';
figure;
plot(dat)
title(sprintf('Data for eigenvalue %.4f,%.4f', x(1), x(2)))

%==========================================================================



