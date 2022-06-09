%% Purdue AGU plots
% Scripts in different sections to generate various plots used in the AGU
% presentation
%
%
% INPUTS - Previously generated intermediate data files
%
% OUTPUTS - plots!
%
% EXAMPLES - go to a section and press ctr-enter
%
%
% Dependencies
%   .m files, .mat files, and MATLAB products required:(updated on 05-Dec-2019)
%         MATLAB (version 9.4)
%         ControlProject.m
%         time_delay_embed.m
%
%   See also: OTHER_FUNCTION_NAME
%
%
%
% Author: Charles Fieseler
% University of Washington, Dept. of Physics
% Email address: charles.fieseler@gmail.com
% Website: coming soon
% Created: 05-Dec-2019
%========================================

%% 1: Data intro
% Example datasets that display clear features
%   Note: All from Distributed dataset
pp = ControlProject();
% which_dataset = 8;
% fnames = pp.get_borehole_fnames(which_dataset);
%% Explore many waveforms
i_vec = 1:5000;
figure('WindowState', 'maximized')
for i = i_vec
    % i = 1;
    dat = readtable(fnames{i});
    % ind = dat{:,1} > 0; % Remove rise time
    % dat = dat{ind,2}';
    dat = dat{:,1}';
    plot(dat)
    ylim([-1, 1])
    title(sprintf("file %d (%s)", i, fnames{i}))
%     pause(0.01)
    drawnow
end
%% Basic trials for SRA

% Pick one dataset
pp = ControlProject();
which_dataset = 7;
which_file = 5;

fnames = pp.get_borehole_fnames(which_dataset);
dat = readtable(fnames{which_file});
X = dat{:,1}';

% Augment with delays
num_delays = 20;
% P = 5; Q = 10;
% X = resample(X, P, Q);
X = time_delay_embed(X, num_delays);
%% Control signals
% Returns a data class with control signals of increasing sparsity
settings = struct('r_ctr',1, ... % Number of control signals to search for
                  'verbose',false, ...
                  'only_positive_U', true,...
                   'num_iter',100); % Number of iterations
control_signal_path = learn_control_signals(X, settings);
%% Calculate best control signal
% objective_function = 'aic_window';
% this_window = {20};
% control_signal_path.calc_best_control_signal(...
%     objective_function, this_window);
objective_function = 'acf';
control_signal_path.calc_best_control_signal(objective_function);

best_U = control_signal_path.U;

% Plot
objective_vals = control_signal_path.objective_values;
% 
figure;
plot(objective_vals)
title('Objective function')
xlabel('Iteration')
%% Get reconstruction and plot
% Exact DMDc; uses L2 solution
%   Note: known to be biased for large noise
[A, B, X_reconstruction] = simple_dmdc(X, best_U);

% For spectrogram viewing
tmp = X_reconstruction(1,:);

% TODO: use actual amplitude of controller
tol = mean(abs(X(1,:)));
disp("Calculated noise level:")
disp(tol)

e = eigs(A);
disp("Calculated window:")
disp( log10(tol) / log10(max(abs(e))) )

% Plot Reconstruction along with data
figure('DefaultAxesFontSize', 24);
subplot(211)
plot(X(1,:), 'linewidth', 1)
hold on
plot(X_reconstruction(1,:))
legend("Data", "Reconstruction")
title("DMDc Model; known control signal")

% figure('DefaultAxesFontSize', 24);
subplot(212)
plot(best_U, 'linewidth', 2)
title("Learned Control signal")
%==========================================================================
%% Plot the top several control signals

num_top_vals = 30;


% Plot
objective_vals = control_signal_path.objective_values;
top_U = zeros(length(best_U), 1);
for i = 1:num_top_vals
    [~, ind] = max(objective_vals);
    disp(ind)
    top_U(:,i) = control_signal_path.all_U{ind} + 0.1*i;
    objective_vals(ind) = -Inf;
end


% 
figure;
plot(top_U)
title('Objective function')
xlabel('Iteration')
%%
% plot_std_fill(top_U, 2);


%==========================================================================










%% Plot as a function of downsampling

% Same base dataset
pp = ControlProject();
which_dataset = 7;
which_file = 5;

fnames = pp.get_borehole_fnames(which_dataset);
dat = readtable(fnames{which_file});
dat_raw = dat{:,1}';

% Augment with delays
%% Settings
num_delays = 8;

% Resampling
P_vec = [10, 8, 6, 4, 2]; 
Q = 10;

settings = struct('r_ctr',1, ... % Number of control signals to search for
                  'verbose',false, ...
                  'only_positive_U', true,...
                   'num_iter',100); % Number of iterations

dat_func = @(p) time_delay_embed(resample(dat_raw, p, Q), num_delays);

all_paths = {};
%% Calculate the controllers for each data set

for i = 1:length(P_vec)
    
    fprintf("Controller for iteration %d / %d \n", i, length(P_vec))
    
    X = dat_func(P_vec(i));
    all_paths{i} = learn_control_signals(X, settings);
end
%% Calculate best path
% objective_function = 'aic_window';
objective_function = 'acf';
% this_window = {};
for i = 1:length(P_vec)
    
%     this_window = {10 * P_vec(i)/Q};
    this_window = {2};
    fprintf("Controller for iteration %d / %d \n", i, length(P_vec))
    all_paths{i}.calc_best_control_signal(...
        objective_function, this_window);
end
%% Plot best control signals
num_top_vals = 10;

% Get values
ts = 1:9978;
top_U = zeros(length(ts), length(all_paths)*num_top_vals);
for i = 1:length(all_paths)
    this_path = all_paths(i);
    objective_vals = all_paths{i}.objective_values;
    
    for i2 = 1:num_top_vals
        % Get single max val
        [~, max_ind] = max(objective_vals);
    %         disp(ind)
        tmp = all_paths{i}.all_U{max_ind};
        objective_vals(max_ind) = -Inf;
    
        % Get corresponding controller
        U_ind = (i-1)*num_top_vals+i2;
        disp(U_ind)
        tmp = resample(tmp, length(ts), length(tmp)) + 0.1*U_ind + 0.5*i;
        
        top_U(:,U_ind) = tmp;
    end
end
%% Plot
figure;
plot(top_U)
title('Objective function')
xlabel('Iteration')


%==========================================================================






%% Plot for several data files

% Same base dataset
pp = ControlProject();
which_dataset = 7;

fnames = pp.get_borehole_fnames(which_dataset);
% which_file = 1:length(fnames)/2;
which_file = [2];

% Augment with delays
%% Settings
num_delays = 8;

settings = struct('r_ctr',1, ... % Number of control signals to search for
                  'verbose',false, ...
                  'only_positive_U', true,...
                   'num_iter',140); % Number of iterations
               
my_f1 = @(i) readtable(fnames{i});
% my_f2 = @(dat) dat{:,1}';
my_f2 = @(dat) lowpass(dat{:,1}', 0.5);
my_f3 = @(X) time_delay_embed(X, num_delays);
%% Calculate the controllers for each data set

all_paths = {};
all_X = [];
for i = 1:length(which_file)
    
    fprintf("Controller for iteration %d / %d \n", i, length(which_file))
    
    X = my_f1(which_file(i));
    X = my_f2(X);
    all_X = [all_X; X / max(max(X)) + i];
    X = my_f3(X);
    all_paths{i} = learn_control_signals(X, settings);
end
%% Calculate best path
objective_function = 'aic_window';
% objective_function = 'acf';
% this_window = {};
for i = 1:length(which_file)
    
%     this_window = {10 * P_vec(i)/Q};
    this_window = {20};
    fprintf("Controller for iteration %d / %d \n", i, length(which_file))
    all_paths{i}.calc_best_control_signal(...
        objective_function, this_window);
end
%% Plot best control signals
num_top_vals = 30;

% Get values
ts = 1:9990;
top_U = zeros(length(ts), length(all_paths)*num_top_vals);
for i = 1:length(all_paths)
    objective_vals = all_paths{i}.objective_values;
    
%     figure;plot(objective_vals)
    
    for i2 = 1:num_top_vals
        % Get single max val
        [~, max_ind] = max(objective_vals);
    %         disp(ind)
        tmp = all_paths{i}.all_U{max_ind};
        objective_vals(max_ind) = -Inf;
    
        % Get corresponding controller
        U_ind = (i-1)*num_top_vals+i2;
        disp(U_ind)
        tmp = tmp + 0.1*U_ind + 0.5*i;
%         tmp = resample(tmp, length(ts), length(tmp)) + 0.1*U_ind + 0.5*i;
        
        top_U(:,U_ind) = tmp;
    end
end
%% Plot
f1 = figure('DefaultAxesFontSize', 24);
plot(top_U)
t1 = sprintf('Top %d control signals (dataset %d)', num_top_vals, which_dataset);
title(t1)
xlabel('Iteration')


f2 = figure('DefaultAxesFontSize', 24);
plot(all_X')
t2 = sprintf('Time series (dataset %d)', which_dataset);
title(t2)
% ylim([0, 1.2])

axes = [f1.Children, f2.Children];
linkaxes(axes, 'x')
%% Save
% savefig(f1, 
fname = string(pp.intermediate_foldername) + strrep(t1, ' ', '-') + ".png";
saveas(f1, fname);

fname = string(pp.intermediate_foldername) + strrep(t2, ' ', '-') + ".png";
saveas(f2, fname);

%==========================================================================



%% Plot for one data file, different low-pass

% Same base dataset
pp = ControlProject();
which_dataset = 7;

fnames = pp.get_borehole_fnames(which_dataset);
% which_file = 1:length(fnames)/2;
which_file = [7];

% Augment with delays
%% Settings
num_delays = 8;

settings = struct('r_ctr',1, ... % Number of control signals to search for
                  'verbose',false, ...
                  'only_positive_U', true,...
                  ...'to_smooth_initialization', true,...
                   'num_iter',140); % Number of iterations
               
               
% lowp_vec = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
lowp_vec = [2,3,4,5];

my_f1 = @(i) readtable(fnames{i});
% my_f2 = @(dat) dat{:,1}';
% my_f2 = @(dat, lowp) lowpass(dat{:,1}', lowp);
% my_f2 = @(dat, lowp) movmean(lowpass(dat{:,1}',0.5), lowp);
my_f2 = @(dat, lowp) lowpass(movmean(dat{:,1}', lowp),0.8);
% my_f2 = @(dat, lowp) padarray(conv(dat{:,1}', gausswin(lowp)), 10-lowp);
my_f3 = @(X) time_delay_embed(X, num_delays);
%% Calculate the controllers for each data set

all_paths = {};
all_X = [];
for i = 1:length(lowp_vec)
    
    fprintf("Controller for iteration %d / %d \n", i, length(lowp_vec))
    
    X = my_f1(which_file(1));
    X = my_f2(X, lowp_vec(i));
    all_X = [all_X; X / max(max(X)) + i];
    X = my_f3(X);
    all_paths{i} = learn_control_signals(X, settings);
end
%% Calculate best path
objective_function = 'aic_window';
% objective_function = 'acf';
% this_window = {};
for i = 1:length(lowp_vec)
    
%     this_window = {10 * P_vec(i)/Q};
    this_window = {50};
    fprintf("Controller for iteration %d / %d \n", i, length(lowp_vec))
    all_paths{i}.calc_best_control_signal(...
        objective_function, this_window);
end
%% Plot best control signals
num_top_vals = 30;

% Get values
ts = 1:9990;
top_U = zeros(length(ts), length(all_paths)*num_top_vals);
for i = 1:length(all_paths)
    objective_vals = all_paths{i}.objective_values;
    
%     figure;plot(objective_vals)
    
    for i2 = 1:num_top_vals
        % Get single max val
        [~, max_ind] = max(objective_vals);
    %         disp(ind)
        tmp = all_paths{i}.all_U{max_ind};
        objective_vals(max_ind) = -Inf;
    
        % Get corresponding controller
        U_ind = (i-1)*num_top_vals+i2;
        disp(U_ind)
        tmp = tmp + 0.1*U_ind + 0.5*i;
%         tmp = resample(tmp, length(ts), length(tmp)) + 0.1*U_ind + 0.5*i;
        
        top_U(:,U_ind) = tmp;
    end
end
%% Plot
f1 = figure('DefaultAxesFontSize', 24);
plot(top_U)
t1 = sprintf('Top %d control signals (dataset %d)', num_top_vals, which_dataset);
title(t1)
xlabel('Iteration')


f2 = figure('DefaultAxesFontSize', 24);
plot(all_X')
t2 = sprintf('Time series (dataset %d)', which_dataset);
title(t2)
% ylim([0, 1.2])

axes = [f1.Children, f2.Children];
linkaxes(axes, 'x')
%% Save
% savefig(f1, 
fname = string(pp.intermediate_foldername) + strrep(t1, ' ', '-') + ".png";
saveas(f1, fname);

fname = string(pp.intermediate_foldername) + strrep(t2, ' ', '-') + ".png";
saveas(f2, fname);

%==========================================================================







%% RETURN TO CHVEN


%% Plot for several data files

% Same base dataset
pp = ControlProject();
% which_dataset = 'mortar_fnames';
which_dataset = 'distributed_fnames';
% which_dataset = 'localized_fnames';
fnames = pp.(which_dataset);
% [all_dat, kept_ind] = pp.filter_by_activity(fnames);
% fnames = fnames(kept_ind);

% which_file = 1:50;
which_file = 101:251;
% which_file = [6];
%% Settings
num_delays = 8;

settings = struct('r_ctr',1, ... % Number of control signals to search for
                  'verbose',false, ...
                  'only_positive_U', true,...
                   'num_iter',140); % Number of iterations
               
my_f1 = @(i) readtable(fnames{i});
ind = 100:500;
my_f2 = @(dat) dat{ind,2}';
% my_f2 = @(dat) lowpass(dat{:,2}', 0.8);
my_f3 = @(X) time_delay_embed(X, num_delays);
%% Calculate the controllers for each data set

all_paths = {};
all_X = [];
for i = 1:length(which_file)
    
    fprintf("Controller for iteration %d / %d \n", i, length(which_file))
    
    X = my_f1(which_file(i));
    X = my_f2(X);
    all_X = [all_X; X];
    X = my_f3(X);
    all_paths{i} = learn_control_signals(X, settings);
end
%% Calculate best path
objective_function = 'aic_window';
% objective_function = 'acf';

for i = 1:length(which_file)
    
    this_window = {3};
    fprintf("Controller for iteration %d / %d \n", i, length(which_file))
    all_paths{i}.calc_best_control_signal(...
        objective_function, this_window);
end

%% Summary function for how well SRA does
all_var_explained = zeros(length(which_file),1);
all_X_reconstructed = [];
all_num_events = all_var_explained;

for i = 1:length(which_file)
    % Calculate variance explained using the best controller
    best_U = all_paths{i}.U;
    best_U = remove_isolated_spikes(best_U);
    X = my_f1(which_file(i));
    X = my_f2(X);
    this_X = my_f3(X);
    % Get accuracy
    [this_var, X_reconstruction] = var_explained_by_dmdc(this_X, best_U);
    all_var_explained(i) = this_var;
    all_X_reconstructed = [all_X_reconstructed; X_reconstruction(end,:)];
    % Get number of control signals
    all_num_events(i) = length(calc_contiguous_blocks(...
        logical(best_U)));
end
%% Save best control signals
num_top_vals = 1;

% Get values
ts = 1:(size(all_X,2)-num_delays-1);
top_U = zeros(length(ts), length(all_paths)*num_top_vals);
offsets = zeros(1,length(all_paths)*num_top_vals);
for i = 1:length(all_paths)
    objective_vals = all_paths{i}.objective_values;
    
%     figure;plot(objective_vals)
    
    for i2 = 1:num_top_vals
        % Get single max val
        [~, max_ind] = max(objective_vals);
    %         disp(ind)
        tmp = all_paths{i}.all_U{max_ind};
        objective_vals(max_ind) = -Inf;
    
        % Get corresponding controller
        U_ind = (i-1)*num_top_vals+i2;
        disp(U_ind)
        if num_top_vals == 1
            offsets(U_ind) = U_ind;
        else
            offsets(U_ind) = i2/(num_top_vals + 10) + i;
        end

        top_U(:,U_ind) = tmp;
    end
end
%% Plot
f1 = figure('DefaultAxesFontSize', 24);
% plot(remove_isolated_spikes(top_U')'+offsets)
plot(top_U + offsets)
t1 = sprintf('Top %d control signals (dataset %s)', num_top_vals, which_dataset);
title(t1)
xlabel('Iteration')


f2 = figure('DefaultAxesFontSize', 24);
plot((0.5*all_X ./ max(all_X, [], 2))' + offsets)
t2 = sprintf('Time series (dataset %s)', which_dataset);
title(t2)
xlabel('Iteration')
% ylim([0, 1.2])


% f3 = figure('DefaultAxesFontSize', 24);
% plot(all_X' + offsets)
% hold on
% plot(all_X_reconstructed' + offsets)
% legend(num2str(all_var_explained))

figure;
violinplot(all_var_explained);

figure;
scatter(all_var_explained, all_num_events)
xlabel('Correlation')
ylabel('NNZ')

axes = [f1.Children, f2.Children];
linkaxes(axes, 'xy')
%% Look at the reconstructions one-by-one
figure
for i = 1:length(all_var_explained)
    plot(all_X(i,:))
    hold on
    plot(all_X_reconstructed(i,:))
    hold off
    title(sprintf("Percent %.2f; nnz %d; file %d", ...
        all_var_explained(i), all_num_events(i), which_file(i)))
    pause
end
%% Look at residuals one-by-one
figure
for i = 1:length(all_var_explained)
    subplot(211)
    plot(all_X(i,:))
    hold on
    plot(all_X_reconstructed(i,:))
    hold off
    xlim([1,400])
    title(sprintf("Percent %.2f; nnz %d", ...
        all_var_explained(i), all_num_events(i)))
    subplot(212)
    plot(all_X(i,1:size(all_X_reconstructed,2)) - all_X_reconstructed(i,:))
    title('Residual')
    xlim([1,400])
    pause
end
%% Save
% savefig(f1, 
fname = string(pp.intermediate_foldername) + strrep(t1, ' ', '-') + ".png";
saveas(f1, fname);

fname = string(pp.intermediate_foldername) + strrep(t2, ' ', '-') + ".png";
saveas(f2, fname);

%==========================================================================




%% Side test: correcting eigenvalues
% Raw, biased eigenvalues (noisy data)
i = 1;
best_U = all_paths{i}.U;
this_X = my_f3(all_X(i,:));
[A, B, xx1] = simple_dmdc(this_X, best_U);
e1 = abs(eig(A));

% eigenvalues with an offset
[A2, B2, xx2] = simple_dmdc(this_X+100, best_U);
e2 = abs(eig(A2));

% Plot
figure;
plot(this_X(1,:)+100)
hold on
plot(xx2(1,:))

figure;
plot(this_X(1,:))
hold on
plot(xx1(1,:))
%% Side test: using _data script
%% Test: density plot
X = [accuracy_d, num_events_d];
X2 = [accuracy_m, num_events_m];
X = [X; X2];
% From: https://www.mathworks.com/matlabcentral/fileexchange/65166-densityplot-x-y-varargin
% And: https://www.mathworks.com/matlabcentral/answers/225934-using-matlab-to-plot-density-contour-for-scatter-plot
[N, C] = hist3(X);
% [N, C] = hist3([accuracy_m, num_events_m]);
wx=C{1}(:);
wy=C{2}(:);
% display
figure
H = pcolor(wx, wy, N');
% box on
hold on
contour(C{1}, C{2}, N', 'k');
shading interp
colormap hot
ylim([1.6,5])
%% Test: GMM model
X = [accuracy_d, num_events_d/ 10];
X2 = [accuracy_m, num_events_m/ 10];
X = [X; X2];

figure
% scatter(X(:,1), X(:,2));
[N, C] = hist3(X);
wx=C{1}(:);
wy=C{2}(:);
H = pcolor(wx, wy, N');
box on
hold on
% shading interp

% Add gmm
gmm = fitgmdist(X,6, 'SharedCov',true);
% gmm = fitgmdist(X,10, 'RegularizationValue', 0.01);
% To converge: https://www.mathworks.com/matlabcentral/answers/100210-why-do-i-receive-an-error-while-trying-to-generate-the-gaussian-mixture-parameter-estimates-from-a-d
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gmm,[x0 y0]),x,y);
g = gca;
fcontour(gmPDF,[g.XLim g.YLim], 'r')
%%
% [best_gmm,all_AIC,all_gmm] = choose_cluster_num(X, 7, {'SharedCov',true});
opt = {'CovarianceType','diagonal', ...
    'RegularizationValue', 0.01,...
    'SharedCov',true};
[gmm,all_AIC,all_gmm] = choose_cluster_num(X, 10, opt);







%% NEW: 10/11/2021

%% Return to power spectrum: can you tell events apart?
%% Data
pp = ControlProject();
% which_dataset = 'mortar_fnames';
which_dataset = 'distributed_fnames';
% which_dataset = 'localized_fnames';
fnames = pp.(which_dataset);
% [all_dat, kept_ind] = pp.filter_by_activity(fnames);
% fnames = fnames(kept_ind);

which_file = 101:251;

read_from_file = @(i) readtable(fnames{i});
ind = 100:500;
extract_time_series = @(dat) dat{ind,2}';
% my_f2 = @(dat) lowpass(dat{:,2}', 0.8);
% my_f3 = @(X) time_delay_embed(X, num_delays);
%% Looping
all_paths = {};
all_X = [];
for i = 1:length(which_file)
    
    fprintf("Iteration %d / %d \n", i, length(which_file))
    
    X = read_from_file(which_file(i));
    X = extract_time_series(X);
    all_X = [all_X; X];
%     X = my_f3(X);
%     all_paths{i} = learn_control_signals(X, settings);
end

%% 



%%