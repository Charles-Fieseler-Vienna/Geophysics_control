

%% Reviewing whole-error metrics
%---------------------------------------------
% % Get data, possibly filtered
%---------------------------------------------
pp = PurdueProject();
trim = 500;

aug = 20;
num_datasets = 50;
min_activity_thresh = 0.015;
kept_dat_ind = false(num_datasets, 1);
fnames = pp.localized_fnames;
% fnames = pp.mortar_fnames;
dat = cell(num_datasets, 1);
all_times = dat;
for i = 1:num_datasets
    dat{i} = readtable(fnames{i});
    all_times{i} = dat{i}{1:end-trim, 1}';
    dat{i} = dat{i}{1:end-trim,2}';
    kept_dat_ind(i) = (max(dat{i}) > min_activity_thresh);
end

dat = dat(kept_dat_ind);
num_datasets = length(dat);

% Time-delay embed
for i = 1:num_datasets
    dat{i} = time_delay_embed(dat{i}, aug);
end

%---------------------------------------------
%% Look at the errors over the entire reconstruction
%---------------------------------------------
this_dat = dat{1};
s.r_ctr = 1;
s.verbose = false;
s.only_positive_U = false;
s.num_iter = 100;
max_func = 'acf';

[all_U, all_A, all_B] = sparse_residual_analysis(this_dat, s);

all_err = zeros(s.num_iter, 1);
all_err_large = all_err;
all_recon = zeros(s.num_iter, size(this_dat, 2));
all_BUnorm = all_err;
all_p = all_err;
all_acf = all_err;
all_nnz = all_err;
x0 = this_dat(:, 1);
for i = 1:s.num_iter
    U = all_U{i};
    tmp = calc_reconstruction_dmd(x0, [], all_A{i}, all_B{i}, U);
    all_recon(i,:) = tmp(1,:);
    all_BUnorm(i) = norm(all_B{i}*U, 'fro');
    all_nnz(i) = nnz(U);
    this_err = tmp - this_dat;
    all_err(i) = norm(this_err, 'fro') / length(x0);
    all_err_large(i) = norm(this_err(abs(this_err)>noise_level), 'fro')...
         / length(x0);
%     [~, all_p(i)] = lillietest(this_err(1,:));
    all_acf(i) = acf(U', 1, false);
end

% figure;plot(all_err); title('All errors')
% figure;plot(all_err_large); title('All errors above the noise threshold')
% figure;plot(all_BUnorm); title('Magnitude of the control signal')
% 
% figure;plot(all_err + all_BUnorm); title('Errors + control signal norm')
% figure;plot(all_acf); title('Autocorrelation')

%---------------------------------------------
% Get an estimate for the noise level
%---------------------------------------------
num_iter = 100;
X = this_dat;
X1 = X(:, 1:end-1);
X2 = X(:, 2:end);
raw_err = X2 - (X2/X1)*X1;
raw_err = raw_err(1,:);
err = raw_err;
all_p = zeros(num_iter, 1);
for i = 1:num_iter
    [~, ind] = max(abs(err));
    val = err(ind);
    err(ind) = 0;
    [h, all_p(i)] = lillietest(err);
end
p_thresh = 0.05;
good_ind_max = find(all_p > p_thresh, 1, 'last');
[~, good_ind_all] = maxk(abs(raw_err), good_ind_max);
noise = zeros(size(err));
bad_ind_all = 1:length(err);
bad_ind_all(good_ind_all) = [];
noise(bad_ind_all) = raw_err(bad_ind_all);

noise_level = max(noise);

%% Look at individual controllers and reconstructions
figure;
i = 60;
subplot(2,1,1);plot(this_dat(1,:)); hold on; plot(all_recon(i,:))
subplot(2,1,2);plot(all_U{i})

%% Post-processing loop
% motivation: the initial hard thresholding is a very blunt tool that
% should be refined as the signal selection has to decide between less
% obviously noise terms

%---------------------------------------------
% Get the starting point
%---------------------------------------------
% First, try to automatically find a good starting point

min_starting_entries = 20;
max_starting_entries = 200;
% start_mode = 'error';
% start_mode = 'acf';
start_mode = 'early';
burnin = [30, 5];
if strcmp(start_mode, 'error')
    tmp = all_err + all_BUnorm;
    [~, min_ind] = min(tmp(burnin(1):end-burnin(2)));
    start_ind = min_ind + burnin(1) - 1;
elseif strcmp(start_mode, 'acf')
    [~, start_ind] = max(all_acf);
elseif strcmp(start_mode, 'early')
    start_ind = burnin(1);
end
U0 = all_U{start_ind};

% Catch if above logic produced a bad start
if ...true ||...
        start_ind==burnin(1) || start_ind==length(tmp)-burnin(2) ||...
        nnz(U0)<=min_starting_entries || nnz(U0) <= nnz(all_U{end})
    warning('Automatic detection of starting point failed; using default')
    start_ind = find(all_nnz <= max_starting_entries, 1);
    U0 = all_U{start_ind};
end

fprintf('Starting at control signal %d with %d entries\n',...
    start_ind, nnz(U0));

% start_ind = 91;

% figure;plot(burnin:length(all_err), all_err(burnin:end)); hold on;
% plot(start_ind, all_err(start_ind), '*')
% title(sprintf('Starting at minimum value of raw error (%d)',...
%     start_ind))

n = nnz(U0) - 1;
new_U = cell(n, 1);
new_U{1} = U0;
all_new_err = cell(n,1);
x0 = this_dat(:,1);
A = all_A{start_ind};
B = all_B{start_ind};
tmp = calc_reconstruction_dmd(x0, [], A, B, U0);
base_err_mat = tmp - this_dat;

%---------------------------------------------
% Settings for how to calculate the improvement
%---------------------------------------------
% Particularly for time-delay embedded data, should we count all rows?
only_use_first_row = true;
use_mean_of_rows = false;
% Should we count errors below the noise level?
hard_threshold_noise_level = false;
% an improvement of 0 is not black and white
improvement_threshold = 0;
% A false kick won't change the error in the entire time series, but only
% as far as the largest eigenvalue
calc_kick_range = false;
% Instead of comparing changes in global error and implying that each
% individual kick is informed by all the data, maybe we can use a small
% window to calculate local aicc values instead
%   Note: this means the base error won't make sense and will have to be
%   calculated for each entry separately (i.e. just keep the matrix)
use_local_aicc = false;
% fast_descent_threshold = 2;
% err_mode = 'raw';
% err_mode = 'BU';
% err_mode = 'aicc';
err_mode = 'dyn_to_ctr'; % Ratio of dynamic reconstruction to ctr signals

if true %calc_kick_range || use_local_aicc
    % Need to remove kicks that are below the noise level
    for i2 = 1:5
        ctr_ind = find(U0);
        to_remove = false(size(ctr_ind));
        for i = 1:length(ctr_ind)
            if abs(B*U0(ctr_ind(i))) < noise_level
                to_remove(i) = true;
            end
        end
        U0(ctr_ind(to_remove)) = 0;
        [A, B] = exact_dmdc(this_dat, U0);
        if isempty(to_remove)
            break
        end
    end
    fprintf('Removal of noise-level kicks leaves %d entries\n', ...
        nnz(U0));
    new_U{1} = U0;
    n = nnz(U0) - 1;
end
if only_use_first_row
    base_err_mat = base_err_mat(1,:);
end
if use_mean_of_rows
    err_mat_factor_base = 1/length(x0);
else
    err_mat_factor_base = 1;
end
if hard_threshold_noise_level
    base_err_mat(abs(base_err_mat) < noise_level) = 0;
end
if strcmp(err_mode, 'raw')
    base_err = norm(base_err_mat, 'fro');
elseif strcmp(err_mode, 'BU')
    base_err = norm(base_err_mat, 'fro') + norm(B*U0, 'fro');
elseif strcmp(err_mode, 'aic')
    base_err = log(norm(base_err_mat, 'fro'))*size(U0, 2)*err_mat_factor_base...
        + 2*nnz(U0);
elseif strcmp(err_mode, 'aicc')
    B(abs(B)<1e-6) = 0;
    A(abs(A)<1e-6) = 0;
    k_t = nnz(U0) + nnz(A) + nnz(B); % Total number of parameters
    correction = (2*k_t.^2 + 2*k_t) / abs(size(U0, 2) - k_t - 1);
%     aic = log(norm(err_mat, 'fro'))*(size(U0, 2)-nnz(U0))...
%         + 2*nnz(U0);
    if calc_kick_range
        % i.e. how many steps will a kick influence, above the
        % noise level (usually ~200)
        %   Note that many kicks are themselves nearly at the noise
        %   level, i.e. B*U ~ noise_level
        %   Note that this is a an approximation, averaged over all control
        %   signals
        BU = B*U0;
        BU = abs(BU(1,:));
        BU = mean(BU(BU>0));
        err_mat_factor = max( log(noise_level/BU) / ...
            log(abs(eigs(A, 1))) * err_mat_factor_base, 0);
    else
        err_mat_factor = err_mat_factor_base*size(U0, 2);
    end
    aic = log(norm(base_err_mat, 'fro'))*err_mat_factor...
        + 2*nnz(U0);
    base_err = aic + correction;
elseif strcmp(err_mode, 'dyn_to_ctr')
    base_err = norm(base_err_mat(1,:), 'fro') + norm(B(1)*U0, 'fro');
end

%---------------------------------------------
% Iterative algorithm
%---------------------------------------------
% Throw out the signal that improves the full reconstruction the least
for i_U = 1:n
    U = new_U{i_U};
    num_ctr = nnz(U);
    this_improvement = zeros(num_ctr, 1);
    ctr_ind = find(U);
    for i_ctr = 1:num_ctr
        % Remove the test control signal
        test_U = U;
        test_U(ctr_ind(i_ctr)) = 0;
        
        % Get dynamics matrices and reconstruction
        [A, B] = exact_dmdc(this_dat, test_U);
        tmp = calc_reconstruction_dmd(x0, [], A, B, test_U);
        err_mat = tmp - this_dat;
        if only_use_first_row
            err_mat = err_mat(1,:);
            base_err_mat = base_err_mat(1,:);
        end
        if hard_threshold_noise_level
            err_mat(abs(err_mat) < noise_level) = 0;
        end
        if use_local_aicc
            BU = abs(B*U(ctr_ind(i_ctr)));
            range = round(max( log(noise_level/BU(1)) / ...
                log(abs(eigs(A, 1))), 0));
            range_ind = ctr_ind(i_ctr):ctr_ind(i_ctr) + range;
            base_err = log(norm(base_err_mat(:, range_ind), 'fro'))*range + 2;
            this_err = log(norm(err_mat(:, range_ind), 'fro'))*range;
            this_improvement(i_ctr) = base_err - this_err;
            continue
        end
        if strcmp(err_mode, 'raw')
            this_improvement(i_ctr) = base_err - norm(err_mat, 'fro');
        elseif strcmp(err_mode, 'BU')
            this_BU_magnitude = norm(B*U(ctr_ind(i_ctr)), 'fro');
    %         this_improvement(i_ctr) = (base_err - norm(err_mat, 'fro')) ...
    %             / length(x0) + this_BU_magnitude;
            this_improvement(i_ctr) = (base_err - norm(err_mat, 'fro')) ...
                + this_BU_magnitude;
        elseif strcmp(err_mode, 'aic')
            aic = log(norm(err_mat, 'fro'))*size(test_U, 2)*err_mat_factor_base...
                + 2*nnz(test_U);
            this_improvement(i_ctr) = base_err - aic;
        elseif strcmp(err_mode, 'aicc')
            % 2nd order correction
            B(abs(B)<1e-6) = 0;
            A(abs(A)<1e-6) = 0;
            k_t = nnz(test_U) + nnz(A) + nnz(B); % Total number of parameters
            correction = (2*k_t.^2 + 2*k_t) / abs(size(U0, 2) - k_t - 1);
%             aic = log(norm(err_mat, 'fro'))*(size(test_U, 2)-nnz(test_U)) ...
%                 + 2*nnz(test_U);
            if calc_kick_range
                % i.e. how many steps will a kick influence, above the
                % noise level (usually ~200)
                %   Note that many kicks are themselves nearly at the noise
                %   level, i.e. B*U ~ noise_level
                BU = abs(B*U(ctr_ind(i_ctr)));
                err_mat_factor = max( log(noise_level/BU(1)) / ...
                    log(abs(eigs(A, 1))) * err_mat_factor_base, 0);
            else
                err_mat_factor = err_mat_factor_base*size(test_U, 2);
            end
            aic = log(norm(err_mat, 'fro'))*err_mat_factor ...
                + 2*nnz(test_U);
            aicc = aic + correction;
            this_improvement(i_ctr) = base_err - aicc;
        elseif strcmp(err_mode, 'dyn_to_ctr')
            this_quality = norm(err_mat(1,:), 'fro') + ...
                norm(B(1)*test_U, 'fro');
            this_improvement(i_ctr) = base_err - this_quality;
        end
            
    end
    % Did any removals improve the error?
    [val, ind] = max(this_improvement);
    % Sanity check
    if max(this_improvement) - min(this_improvement) < 1e-6
        warning('Almost no difference between different removals')
    end
    fprintf('Iteration %d/%d:\n', i_U, n)
    fprintf('The best improvement is %.5f via removal of signal at %d\n',...
        val, ctr_ind(ind))
    if val < improvement_threshold
        disp('Improvement is below threshold, aborting')
        U_best = U;
        new_U{i_U+1} = U;
        break
    end
    if calc_kick_range
        %   Note that this is a an approximation, averaged over all 
        % control signals
        %   NOTE: only AICc
        [A, B] = exact_dmdc(this_dat, U);
        BU = B*U;
        BU = abs(BU(1,:));
        BU = mean(BU(BU>0));
        err_mat_factor = max( log(noise_level/BU) / ...
            log(abs(eigs(A, 1))) * err_mat_factor_base, 0);

        B(abs(B)<1e-6) = 0;
        A(abs(A)<1e-6) = 0;
        k_t = nnz(U) + nnz(A) + nnz(B); % Total number of parameters
        correction = (2*k_t.^2 + 2*k_t) / abs(size(U, 2) - k_t - 1);

        aic = log(norm(err_mat, 'fro'))*err_mat_factor ...
            + 2*nnz(U);
        base_err = aic + correction;
    else
        base_err = base_err - val;
    end
    if use_local_aicc
        [A, B] = exact_dmdc(this_dat, U);
        tmp = calc_reconstruction_dmd(x0, [], A, B, U);
        base_err_mat = tmp - this_dat;
    end
    U(ctr_ind(ind)) = 0;
    new_U{i_U+1} = U;
end

%---------------------------------------------
% Look at individual controllers and reconstructions
%---------------------------------------------
figure;
i = i_U;
try
    U = new_U{i};
catch
    U = new_U{1};
end
[A, B] = exact_dmdc(this_dat, U);
tmp = calc_reconstruction_dmd(x0, [], A, B, U);
subplot(3,1,1);plot(this_dat(1,:)); hold on; plot(tmp(1,:))
dyn_to_ctr = norm(tmp(1,:), 'fro') / norm(B(1)*U, 'fro');
title(sprintf('Reconstruction (quality is %.2f)', dyn_to_ctr))

subplot(3,1,2);plot(U)
title('Processed control signal')
subplot(3,1,3);plot(U0)
title(sprintf('Initial control signal (%d/%d entries removed)', i_U-1, n))
%==========================================================================

%% Using Lillietest to get error estimate
% A completely different way to do signal learning:
%   Start with the raw residual
%   Iteratively remove large values, doing a p-test to see how Gaussian the rest of the residual is
%   Keep the control signal that, when removed from the residual, makes the remaining residual the most Gaussian
n = 100;
X = this_dat;
X1 = X(:, 1:end-1);
X2 = X(:, 2:end);
raw_err = X2 - (X2/X1)*X1;
raw_err = raw_err(1,:);
err = raw_err;
all_p = zeros(n, 1);
for i = 1:n
    [~, ind] = max(abs(err));
    val = err(ind);
    err(ind) = 0;
    [h, all_p(i)] = lillietest(err);
end
p_thresh = 0.05;
good_ind_max = find(all_p > p_thresh, 1, 'last');
[~, good_ind_all] = maxk(abs(raw_err), good_ind_max);
learned_U = zeros(size(err));
learned_U(good_ind_all) = raw_err(good_ind_all);
noise = zeros(size(err));
bad_ind_all = 1:length(err);
bad_ind_all(good_ind_all) = [];
noise(bad_ind_all) = raw_err(bad_ind_all);

noise_level = max(noise);
% figure;plot(all_p)
% figure;plot(learned_U)

%---------------------------------------------
%% Calc control 'realness' for these test signals
%---------------------------------------------
settings = define_ideal_settings();
settings.dmd_mode = 'naive';
settings.augment_data = 0;
settings.custom_control_signal = [learned_U 0];
my_model = CElegansModel(this_dat, settings);

[all_err, out] = calc_control_realness(my_model);

my_model.plot_reconstruction_interactive(true, 1);
%==========================================================================


%% Scratch for post processing

%     if fast_descent_threshold > 0
%         all_removals = this_improvement > fast_descent_threshold;
%         if ~isempty(all_removals)
%             fprintf('Fast removal of %d indices\n',...
%                 length(find(all_removals)));
%             U(ctr_ind(all_removals)) = 0;
%             new_U{i_U+1} = U;
%             % TODO: update base_err
%             B(abs(B)<1e-6) = 0;
%             A(abs(A)<1e-6) = 0;
%             k_t = nnz(U) + nnz(A) + nnz(B); % Total number of parameters
%             correction = (2*k_t.^2 + 2*k_t) / abs(size(U0, 2) - k_t - 1);
%             aic = log(norm(base_err_mat, 'fro'))*size(U0, 2)*err_mat_factor...
%                 + 2*nnz(U);
%             base_err = aic + correction;
%             continue
%         end
%     end





%==========================================================================


%% Get signals and postprocess them for an entire dataset

% Import data
disp('Importing and preprocessing...')
pp = PurdueProject();
fnames = pp.mortar_fnames;
% fnames = pp.localized_fnames;

% Preprocess, i.e. filter and time-delay embed
% dat = pp.filter_by_activity(fnames(1:100));
dat = pp.filter_by_activity(fnames(1:2));

aug = 20;
n = length(dat);
for i = 1:n
    dat{i} = time_delay_embed(dat{i}, aug);
end

% Get initial control signals, and then post-process
all_U_best = cell(n,1);

sra_settings.r_ctr = 1;
sra_settings.verbose = false;
sra_settings.only_positive_U = false;
sra_settings.num_iter = 70;

post_settings.row_ind = 1;
post_settings.verbose = 1;

for i = 1:n
    fprintf('Getting control signals for dataset %d...\n', i)
    this_dat = dat{i};
    [all_U, all_A, all_B] = sparse_residual_analysis(this_dat, sra_settings);
    fprintf('Postprocessing control signals for dataset %d...\n', i)
    [U_best] = postprocess_control_signals(this_dat, all_U, all_A, all_B,...
        post_settings);
    all_U_best{i} = U_best;
end
disp('====================================================')
disp('Finished')
disp('====================================================')

plot_cell_array(all_U_best);

%==========================================================================








