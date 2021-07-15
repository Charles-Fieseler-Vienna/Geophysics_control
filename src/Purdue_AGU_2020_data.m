%% SHARED
%% Settings
pp = PurdueProject();

num_delays = 8;

% Building control signals
settings = struct('r_ctr',1, ... % Number of control signals to search for
                  'verbose',false, ...
                  'only_positive_U', true,...
                   'num_iter',130); % Number of iterations
objective_function = 'aic_window';
this_window = {20};
force_nonempty = true;
% objective_function = 'acf';

% Preprocessing functions
ind = 100:500;
my_f2 = @(dat) dat{ind,2}';
% my_f2 = @(dat) lowpass(dat{:,2}', 0.8);
my_f3 = @(X) time_delay_embed(X, num_delays);

error("Please run sections individually")


%% MORTAR
%%

%% Analyze several data files

% Same base dataset
which_dataset = 'mortar_fnames';
% which_dataset = 'distributed_fnames';
% which_dataset = 'localized_fnames';
fnames = pp.(which_dataset);
% [all_dat, kept_ind] = pp.filter_by_activity(fnames);
% fnames = fnames(kept_ind);
my_f1 = @(i) readtable(fnames{i});

this_window = {20};

% which_file = 1:length(fnames);
which_file = 1:250;
% which_file = 101:251;
% which_file = [6];
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

for i = 1:length(which_file)
    
    fprintf("Controller for iteration %d / %d \n", i, length(which_file))
    all_paths{i}.calc_best_control_signal(...
        objective_function, this_window, force_nonempty);
end
%% Summary function for how well SRA does
all_var_explained = zeros(length(which_file),1);
all_X_reconstructed = [];
all_num_events = all_var_explained;

for i = 1:length(which_file)
    % Calculate variance explained using the best controller
    best_U = all_paths{i}.U;
%     best_U = remove_isolated_spikes(best_U);
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
%% Data 1: Example good datasets

ex1 = 3;
fname = pp.intermediate_foldername + "mortar_ex1.mat";

X = all_X(ex1,:);
X_recon = all_X_reconstructed(ex1,:);
U = all_paths{ex1}.U;

save(fname, 'X', 'X_recon', 'U');
%% Data 2: Another example
ex2 = 5;
fname = pp.intermediate_foldername + "mortar_ex2.mat";

X = all_X(ex2,:);
X_recon = all_X_reconstructed(ex2,:);
U = all_paths{ex2}.U;

save(fname, 'X', 'X_recon', 'U');
%% Data 3: Accuracy
fname = pp.intermediate_foldername + "mortar_acc.mat";

accuracy = all_var_explained;
save(fname, 'accuracy');
%% Data 4: Scatterplot
fname = pp.intermediate_foldername + "mortar_scatter.mat";

accuracy = all_var_explained;
num_events = all_num_events;
save(fname, 'accuracy', 'num_events');
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
    subplot(211)
    plot(all_X(i,:))
    hold on
    plot(all_X_reconstructed(i,:))
    hold off
    title(sprintf("Percent %.2f; nnz %d; file %d", ...
        all_var_explained(i), all_num_events(i), which_file(i)))
    subplot(212)
    plot(all_paths{i}.U)
    pause
end
%% Save
% savefig(f1, 
fname = string(pp.intermediate_foldername) + strrep(t1, ' ', '-') + ".png";
saveas(f1, fname);

fname = string(pp.intermediate_foldername) + strrep(t2, ' ', '-') + ".png";
saveas(f2, fname);

%==========================================================================



%% DISTRIBUTED
%%

%% Analyze several data files

% Same base dataset
which_dataset = 'distributed_fnames';
% which_dataset = 'localized_fnames';
fnames = pp.(which_dataset);
% [all_dat, kept_ind] = pp.filter_by_activity(fnames);
% fnames = fnames(kept_ind);
my_f1 = @(i) readtable(fnames{i});

which_file = 1:length(fnames);
fprintf('%d files remaining\n', length(fnames))
% which_file = 1:50;
% which_file = 101:251;
% which_file = [6];
%% Test: fewer datasets
which_file = 1:250;

% this_window = {10};
%% Calculate the controllers for each data set

all_paths = {};
all_X = [];
for i = 1:length(which_file)
    
    fprintf("Controller for file %d / %d \n", i, length(which_file))
    
    X = my_f1(which_file(i));
    X = my_f2(X);
    all_X = [all_X; X]; % Save before embedding
    X = my_f3(X);
    all_paths{i} = learn_control_signals(X, settings);
end
%% Calculate best path

for i = 1:length(which_file)
    
    fprintf("Controller for file %d / %d \n", i, length(which_file))
    all_paths{i}.calc_best_control_signal(...
        objective_function, this_window, force_nonempty);
end
%% Summary function for how well SRA does
all_var_explained = zeros(length(which_file),1);
all_X_reconstructed = [];
all_num_events = all_var_explained;

for i = 1:length(which_file)
    fprintf("Accuracy for file %d / %d \n", i, length(which_file))

    % Calculate variance explained using the best controller
    best_U = all_paths{i}.U;
%     best_U = remove_isolated_spikes(best_U);
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
%         disp(U_ind)
        if num_top_vals == 1
            offsets(U_ind) = U_ind;
        else
            offsets(U_ind) = i2/(num_top_vals + 10) + i;
        end

        top_U(:,U_ind) = tmp;
    end
end
%% Data 1: Example mortar-like

ex1 = 7;
fname = pp.intermediate_foldername + "distributed_ex1.mat";

X = all_X(ex1,:);
X_recon = all_X_reconstructed(ex1,:);
U = all_paths{ex1}.U;

save(fname, 'X', 'X_recon', 'U');
%% Data 2: Example "clean packet" dataset
ex2 = 19;
fname = pp.intermediate_foldername + "distributed_ex2.mat";

X = all_X(ex2,:);
X_recon = all_X_reconstructed(ex2,:);
U = all_paths{ex2}.U;

save(fname, 'X', 'X_recon', 'U');
%% Data 3: Example "messy packet" dataset
ex2 = 43;
fname = pp.intermediate_foldername + "distributed_ex3.mat";

X = all_X(ex2,:);
X_recon = all_X_reconstructed(ex2,:);
U = all_paths{ex2}.U;

save(fname, 'X', 'X_recon', 'U');
%% Data 4: Accuracy
fname = pp.intermediate_foldername + "distributed_acc.mat";

accuracy = all_var_explained;
save(fname, 'accuracy');
%% Data 5: Scatterplot
fname = pp.intermediate_foldername + "distributed_scatter.mat";

accuracy = all_var_explained;
num_events = all_num_events;
save(fname, 'accuracy', 'num_events');
%% Data 6: Inset example 1
[val, ex1] = max(all_var_explained);
X = all_X(ex1,:);
X_recon = all_X_reconstructed(ex1,:);
U = all_paths{ex1}.U;

fname = pp.intermediate_foldername + "distributed_inset1_max.mat";

save(fname, 'X', 'X_recon', 'U');

% figure;
% subplot(211)
% plot(X)
% hold on
% plot(X_recon)
% title(sprintf("Best recontruction (%.2f)", val))
% 
% subplot(212)
% plot(U)
%% Data 7: Inset example 2
% [val, ex2] = max(all_var_explained .* (all_num_events==1));
% X = all_X(ex2,:);
% X_recon = all_X_reconstructed(ex2,:);
% U = all_paths{ex2}.U;
% 
% fname = pp.intermediate_foldername + "distributed_inset2_max_single_event.mat";
% 
% save(fname, 'X', 'X_recon', 'U');

% figure;
% subplot(211)
% plot(X)
% hold on
% plot(X_recon)
% title(sprintf("Best recontruction (%.2f)", val))
% 
% subplot(212)
% plot(U)
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
    subplot(211)
    plot(all_X(i,:))
    hold on
    plot(all_X_reconstructed(i,:))
    hold off
    title(sprintf("Percent %.2f; nnz %d; file %d", ...
        all_var_explained(i), all_num_events(i), which_file(i)))
    subplot(212)
    plot(all_paths{i}.U)
    pause
end



%% DISTRIBUTED
%%

%% Analyze several data files
which_dataset = 'localized_fnames';
fnames = pp.(which_dataset);
[all_dat, kept_ind] = pp.filter_by_activity(fnames);
fnames = fnames(kept_ind);
my_f1 = @(i) readtable(fnames{i});

% which_file = 1:length(fnames);
fprintf('%d files remaining\n', length(fnames))
%% Test: fewer datasets
% which_file = 1:250;
which_file = 1:length(fnames);

% this_window = {10};
%% Calculate the controllers for each data set

all_paths = {};
all_X = [];
for i = 1:length(which_file)
    
    fprintf("Controller for file %d / %d \n", i, length(which_file))
    
    X = my_f1(which_file(i));
    X = my_f2(X);
    all_X = [all_X; X];
    X = my_f3(X);
    all_paths{i} = learn_control_signals(X, settings);
end
%% Calculate best path

for i = 1:length(which_file)
    
    fprintf("Controller for file %d / %d \n", i, length(which_file))
    all_paths{i}.calc_best_control_signal(...
        objective_function, this_window, force_nonempty);
end
%% Summary function for how well SRA does
all_var_explained = zeros(length(which_file),1);
all_X_reconstructed = [];
all_num_events = all_var_explained;

for i = 1:length(which_file)
    fprintf("Accuracy for file %d / %d \n", i, length(which_file))

    % Calculate variance explained using the best controller
    best_U = all_paths{i}.U;
%     best_U = remove_isolated_spikes(best_U);
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
%         disp(U_ind)
        if num_top_vals == 1
            offsets(U_ind) = U_ind;
        else
            offsets(U_ind) = i2/(num_top_vals + 10) + i;
        end

        top_U(:,U_ind) = tmp;
    end
end
%% Data 4: Accuracy
fname = pp.intermediate_foldername + "localized_acc.mat";

accuracy = all_var_explained;
save(fname, 'accuracy');
%% Data 5: Scatterplot
fname = pp.intermediate_foldername + "localized_scatter.mat";

accuracy = all_var_explained;
num_events = all_num_events;
save(fname, 'accuracy', 'num_events');
%% Data 6: Insets
ex1 = 18;
X = all_X(ex1,:);
X_recon = all_X_reconstructed(ex1,:);
U = all_paths{ex1}.U;
fname = pp.intermediate_foldername + "localized_inset1_good.mat";

save(fname, 'X', 'X_recon', 'U');


ex1 = 34;
X = all_X(ex1,:);
X_recon = all_X_reconstructed(ex1,:);
U = all_paths{ex1}.U;
fname = pp.intermediate_foldername + "localized_inset2_good.mat";

save(fname, 'X', 'X_recon', 'U');

ex1 = 224;
X = all_X(ex1,:);
X_recon = all_X_reconstructed(ex1,:);
U = all_paths{ex1}.U;
fname = pp.intermediate_foldername + "localized_inset1_multi.mat";

save(fname, 'X', 'X_recon', 'U');

ex1 = 311;
X = all_X(ex1,:);
X_recon = all_X_reconstructed(ex1,:);
U = all_paths{ex1}.U;
fname = pp.intermediate_foldername + "localized_inset2_multi.mat";

save(fname, 'X', 'X_recon', 'U');
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
%% Plot as a function of experimental time
good_ind = abs(all_var_explained) < 1.0;
 
figure;
scatter(all_var_explained(good_ind), all_num_events(good_ind), ...
    [],1:length(all_num_events(good_ind)))
colorbar
%% Look at the reconstructions one-by-one
figure
for i = 1:length(all_var_explained)
    subplot(211)
    plot(all_X(i,:))
    hold on
    plot(all_X_reconstructed(i,:))
    hold off
    title(sprintf("Percent %.2f; nnz %d; file %d", ...
        all_var_explained(i), all_num_events(i), which_file(i)))
    subplot(212)
    plot(all_paths{i}.U)
    pause
end