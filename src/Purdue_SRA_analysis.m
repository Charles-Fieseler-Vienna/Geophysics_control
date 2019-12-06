%% Do SRA and compare eigenvalues and control signals
%---------------------------------------------
%% Get filenames
%---------------------------------------------
which_dataset = 'mortar_fnames';
% which_dataset = 'localized_fnames';
fnames = pp.(which_dataset);
if strcmp(which_dataset, 'localized_fnames')
    [all_dat, kept_ind] = pp.filter_by_activity(fnames);
    kept_ind = find(kept_ind); % Note that the mortar data is much cleaner
else
    kept_ind = 1:length(fnames);
end

%% Time-delay embed settings
which_dataset = 12;
dat = readtable(fnames{kept_ind(which_dataset)});
dat = dat{1:1500,2}';
aug = 20;
dat_embed_raw = time_delay_embed(dat, aug);

% Preprocess data: keep only top 2 SVD modes
% for i = [2, 4, 6]
%     dat_embed = svd_truncate(dat_embed_raw, i);
%     figure;plot(dat_embed(1,:));
%     title(sprintf('%d', i))
% end
% dat_embed = svd_truncate(dat_embed_raw, 2);
dat_embed = dat_embed_raw;
%---------------------------------------------
%% Get Signals using time delay
%---------------------------------------------
s = struct();
s.r_ctr = 1;
s.num_iter = 100;
max_func = 'aic';

[all_U, all_A, all_B, my_model_base] = ...
    sparse_residual_analysis(dat_embed, s);
all_nnz = zeros(size(all_U));
all_objective = all_nnz;
for i = 1:length(all_U)
    all_nnz(i) = nnz(all_U{i});
    all_objective(i) = acf(all_U{i}', 1, false);
end
all_objective(end) = all_objective(end-1);

%% Get "best" signals and plot them
figure;
subplot(2,1,1)
plot(20:length(all_objective), all_objective(20:end))
title('Objective function')

[~, i] = max(all_objective(20:end));
i = i + 19;
U = all_U{i};
subplot(2,1,2);
plot(U); title(sprintf('max %s=%d', max_func, i))

%% Calculate A and B matrices with these ideal signals
X2 = dat_embed(:, 2:end);
n = size(X2, 1);

AB = X2/[dat_embed(:, 1:end-1); U];
A = AB(:, 1:n);
B = AB(:, (n+1):end);

% Reconstruction
dat_dmd = zeros(size(dat_embed));
dat_dmd(:,1) = dat_embed(:, 1);
for i = 2:size(dat_embed, 2)
    dat_dmd(:, i) = A*dat_dmd(:, i-1) + B*U(:, i-1);
end

figure;
subplot(2,1,1)
plot(dat_embed(1, :));
hold on; plot(dat_dmd(1, :))
title(sprintf('Reconstruction of dataset %d', which_dataset))

subplot(2,1,2)
plot(U);title('Control signal')

%% Look at eigenvalues of A matrix and compare to regular dmd

% DMDc
eigs = log(diag(eig(A)));
dmdc_eigs_real = [real(eigs(end)); imag(eigs(end))];

% DMD
[ ~, Omega] = dmd( dat_embed, 1, 2);
dmd_eigs_real = [real(Omega(end)); abs(imag(Omega(end)))];

% FFT
[p,f,t] = pspectrum(dat, 'spectrogram', 'FrequencyLimits',[0,1.5]);
[max_f, ind_f] = max(p);
[max_of_max_f, ind_max_f] = max(max_f);
actual_max_frequency = f(ind_f(ind_max_f));

%Plot
figure
plot(dmdc_eigs_real(1,:), dmdc_eigs_real(2,:), 'bo')
hold on
plot(dmd_eigs_real(1,:), dmd_eigs_real(2,:), 'ro')
plot(0, actual_max_frequency, 'ko')
legend('DMDc', 'DMD', 'FFT')
%==========================================================================










%% Do SRA over many files
%---------------------------------------------
%% Get filenames
%---------------------------------------------
% which_dataset = 'mortar_fnames';
which_dataset = 'localized_fnames';
fnames = pp.(which_dataset);
if strcmp(which_dataset, 'localized_fnames')
    [all_dat, kept_ind] = pp.filter_by_activity(fnames);
    kept_ind = find(kept_ind); % Note that the mortar data is much cleaner
else
    kept_ind = 1:length(fnames);
end

%% Time-delay embed settings
aug = 20;
s = struct();
s.r_ctr = 1;
s.num_iter = 100;
max_func = 'aic';

U_across_files = zeros(length(kept_ind), 1479);
dmd_across_files = zeros(length(kept_ind), 2);
dmdc_across_files = zeros(length(kept_ind), 2);
fft_across_files = zeros(length(kept_ind), 1);

for i_dat = 1:length(kept_ind)
    dat = readtable(fnames{kept_ind(i_dat)});
    dat = dat{1:1500,2}';
    dat_embed = time_delay_embed(dat, aug);
    %---------------------------------------------
    % Get Signals using time delay
    %---------------------------------------------
    [all_U, all_A, all_B, my_model_base] = ...
        sparse_residual_analysis(dat_embed, s);
    all_nnz = zeros(size(all_U));
    all_objective = all_nnz;
    for i = 1:length(all_U)
        all_nnz(i) = nnz(all_U{i});
        all_objective(i) = acf(all_U{i}', 1, false);
    end
    all_objective(end) = all_objective(end-1);

    % Get "best" signals and plot them
    [~, i] = max(all_objective(20:end));
    i = i + 19;
    U = all_U{i};
    U_across_files(i_dat, :) = U;

    % Calculate A and B matrices with these ideal signals
    X1 = dat_embed(:, 1:end-1);
    X2 = dat_embed(:, 2:end);
    Lambda = func_DMDc(X1, X2, U, 2, 3);

    % Look at eigenvalues of A matrix and compare to regular dmd
    % DMDc
    eigs = log(Lambda);
    dmdc_across_files(i_dat, :) = ...
        [real(eigs(end)), imag(eigs(end))];

    % DMD
    [ ~, Omega] = dmd( dat_embed, 1, 2);
    dmd_across_files(i_dat, :) = ...
        [real(Omega(end)), abs(imag(Omega(end)))];

    % FFT
    [p,f,t] = pspectrum(dat, 'spectrogram', 'FrequencyLimits',[0,1.5]);
    [max_f, ind_f] = max(p);
    [max_of_max_f, ind_max_f] = max(max_f);
    fft_across_files(i_dat) = f(ind_f(ind_max_f));
end

%% Plot eigenvalues
figure
plot(dmdc_across_files(:,1), abs(dmdc_across_files(:,2)), 'ro')
hold on
plot(dmd_across_files(:,1), dmd_across_files(:,2), 'bo')
plot(0, fft_across_files, 'ko')
legend('DMDc', 'DMD', 'FFT')

title(sprintf('Eigenvalues for %s datasets', which_dataset))

%% SAVE INTERMEDIATE
% fname = pp.intermediate_foldername+"distributed_sra.mat";
% fname = pp.intermediate_foldername+"mortar_sra.mat";
fname = pp.intermediate_foldername+"localized_sra.mat";

save(fname, 'dmdc_across_files', 'U_across_files');

%% Plot histograms of controllers
hist_dat = [];
for i = 1:length(kept_ind)
    hist_dat = [hist_dat ...
        find(U_across_files(i,:))];
end

figure;
histogram(hist_dat, 100)
title(sprintf('Timing of control signals for %s datasets', which_dataset))
%==========================================================================
