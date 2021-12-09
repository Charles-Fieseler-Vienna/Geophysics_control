%% Shared setup
%%
pp = PurdueProject();

offsets = normrnd(0, 0.05, 250, 1);
ts = 0.5*(1:401);
ts2 = 0.5*(1:393);
ts3 = 0.5*(1:392);

error("Please run individual sections")

%% Intro
%%

%% Figure 1: All about mortar
tikz_opt = {'width', '0.4\textwidth', 'height', '0.2\textheight'};
fig_opt = {'DefaultAxesFontSize', 20, 'WindowState', 'Maximized'};
line_opt = {'LineWidth',1};

%% Panel 1/4

f1 = figure(fig_opt{:});

fname = "mortar_ex1";
load(pp.intermediate_foldername + fname + ".mat")

subplot(211)
plot(ts, X, line_opt{:})
hold on
plot(ts2, X_recon, line_opt{:})
% title("Example Mortar Event (Successful Reconstruction)")
leg = legend({'Data', 'Reconstruction'}, 'Location', 'Southeast');
% set(leg,'color','none');
xlim([1,200])
ylabel('Amp (Volts)')

subplot(212)
plot(ts3, U, line_opt{:})
legend('Learned Control Signal')
xlim([1,200])
xlabel('Time (micro seconds)')
ylabel('Amp (A.U.)')

%
out_fname = pp.paper_foldername + "fig1/" + fname;
saveas(f1, out_fname + ".png");

cleanfigure();
matlab2tikz(char(out_fname + ".tex"), tikz_opt{:});

close(f1);
%% Panel 2: another example
f2 = figure(fig_opt{:});

fname = "mortar_ex2";
load(pp.intermediate_foldername + fname + ".mat")

subplot(211)
plot(ts, X, line_opt{:})
hold on
plot(ts2, X_recon, line_opt{:})
% title("Example Mortar Event (Successful Reconstruction)")
% legend({'Data', 'Reconstruction'}, 'Location', 'Southeast')
xlim([1,200])
ylabel('Amp (Volts)')

subplot(212)
plot(ts3, U, line_opt{:})
% legend('Learned Control Signal')
xlim([1,200])
xlabel('Time (micro seconds)')
ylabel('Amp (A.U.)')

%
out_fname = pp.paper_foldername + "fig1/" + fname;
saveas(f2, out_fname + ".png");

cleanfigure();
matlab2tikz(char(out_fname + ".tex"), tikz_opt{:});

close(f2);
%% Panel 3: Reconstruction accuracy
f3 = figure(fig_opt{:});

fname = "mortar_acc";
load(pp.intermediate_foldername + fname + ".mat")

violinplot(accuracy)
% title("Consistent Reconstruction")
ylabel("Reconstruction Accuracy")
xticklabels([])
ylim([0,1])

%
out_fname = pp.paper_foldername + "fig1/" + fname;
saveas(f3, out_fname + ".png");

cleanfigure();
matlab2tikz(char(out_fname + ".tex"), tikz_opt{:});
close(f3);
%% Plot 4: Scatterplot
f4 = figure(fig_opt{:});

fname = "mortar_scatter";
load(pp.intermediate_foldername + fname + ".mat")

scatter(num_events, accuracy, '*', 'LineWidth',2)
ylabel('Reconstruction Accuracy')
xlabel('Number of Control Signals')
% title("Mortar Reconstructions Require Few Control Signals")
ylim([0,1])
xlim([0,5])

%
out_fname = pp.paper_foldername + "fig1/" + fname;
saveas(f4, out_fname + ".png");

cleanfigure();
matlab2tikz(char(out_fname + ".tex"), tikz_opt{:});

close(f4);


%% Figure 1 ALT: violin + examples

f1a = figure(fig_opt{:});
fname = "mortar_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;
subplot(121)
violinplot([acc_m], {"Mortar"})
ylabel("Reconstruction Accuracy")
ylim([0, 0.6])
title("Variance Explained")

% Load all raw data
fname = "mortar_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")

example_ind = [6, 2, 1, 18, 99];
num_rows = length(example_ind);
cmap = colormap(parula(6));

for i = 1:num_rows
    i_subplot = 2*i;
    subplot(num_rows, 2, i_subplot);
    
    i_data = example_ind(i);
    
    X = all_X(i_data,:);
    X_recon = all_X_reconstructed(i_data,:);
    U = all_paths{i_data}.U;
    
    plot(ts, X, line_opt{:}, 'color', cmap(i,:))
    hold on
    plot(ts2, X_recon, line_opt{:})
    xlim([1,200])
    plot(ts3, 0.01* U - max(X), 'k', 'linewidth', 2)
    
    yticks([])
    
    this_accuracy = acc_m(i_data);
    y = this_accuracy;
    
    if i==1
        title("Example time series")
    end
    if i<num_rows
        xticks([])
    else
        xlabel("Time (\mu s)")
        legend({"Time series", "Reconstruction", "Control signal"})
    end
    
    subplot(121)
    plot([0.5, 1.5], [y, y], 'color', cmap(i, :), 'linewidth', 3)
end

%
out_fname = pp.paper_foldername + "fig1/" + "alternate_violin";
saveas(f1a, out_fname + ".png");

%% Distributed
%%

%% Plot 1: Example good dataset
f1 = figure(fig_opt{:});

fname = "distributed_ex1";
load(pp.intermediate_foldername + fname + ".mat")

subplot(211)
plot(ts, X, 'Linewidth', 2)
hold on
plot(ts2, X_recon, 'LineWidth',2)
title("Example Distributed Event (Successful Reconstruction)")
legend({'Data', 'Reconstruction'}, 'Location', 'Southeast')
xlim([1,200])
ylabel('Amplitude (Volts)')

subplot(212)
plot(ts3, U, 'LineWidth',2)
legend('Learned Control Signal')
xlim([1,200])
xlabel('Time (micro seconds)')
ylabel('Amplitude (A.U.)')

saveas(f1, pp.paper_foldername + fname + ".png");
%
%% Plot 2: Example "clean packet" dataset
f2 = figure(fig_opt{:});

fname = "distributed_ex2";
load(pp.intermediate_foldername + fname + ".mat")

subplot(211)
plot(ts, X, 'Linewidth', 2)
hold on
plot(ts2, X_recon, 'LineWidth',2)
title("Example Distributed Event (Poor Reconstruction)")
legend({'Data', 'Reconstruction'}, 'Location', 'Southeast')
xlim([1,200])
ylabel('Amplitude (Volts)')

subplot(212)
plot(ts3, U, 'LineWidth',2)
legend('Learned Control Signal')
xlim([1,200])
xlabel('Time (micro seconds)')
ylabel('Amplitude (A.U.)')

saveas(f2, pp.paper_foldername + fname + ".png");
%
%% Plot 3: Example "messy packet" dataset
f3 = figure(fig_opt{:});

fname = "distributed_ex3";
load(pp.intermediate_foldername + fname + ".mat")

subplot(211)
plot(ts, X, 'Linewidth', 2)
hold on
plot(ts2, X_recon, 'LineWidth',2)
title("Example Distributed Event (Poor Reconstruction)")
xlim([1,200])
ylabel('Amplitude (Volts)')

subplot(212)
plot(ts3, U, 'LineWidth',2)
legend('Learned Control Signal')
xlim([1,200])
xlabel('Time (micro seconds)')
ylabel('Amplitude (A.U.)')

saveas(f3, pp.paper_foldername + fname + ".png");
%
%% Plot 4: Reconstruction accuracy
f3 = figure(fig_opt{:});

fname = "distributed_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_d = accuracy;
fname = "mortar_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;

violinplot([acc_m, acc_d], {"Mortar", "Distributed"})
ylabel("Reconstruction Accuracy")
title("Distributed Events Are Poorly Reconstructed")

saveas(f3, pp.paper_foldername + "distributed_acc" + ".png");
%
%% Plot 5: Scatterplot (with mortar as well)

% Load two datasets
fname = "distributed_scatter";
load(pp.intermediate_foldername + fname + ".mat")
accuracy_d = accuracy;
num_events_d = num_events;

fname = "mortar_scatter";
load(pp.intermediate_foldername + fname + ".mat")
accuracy_m = accuracy;
num_events_m = num_events;

% Finally, plot
f4 = figure(fig_opt{:});
scatter(accuracy_d, num_events_d, 'LineWidth',2, ...
    'markeredgecolor', [0.8500 0.3250 0.0980]	)
ylim([0,5])
xlim([0,1])
xlabel('Reconstruction accuracy')
ylabel('Number of Events')
title('Distributed Data Usually Require Many Control Signals')

f5 = figure(fig_opt{:});
scatter(accuracy_m, num_events_m, '*', 'LineWidth',2)
hold on
scatter(accuracy_d, num_events_d, 'LineWidth',2)
ylim([0,5])
xlim([0,1])
title('There Exist Mortar-Like Events within Distributed')

xlabel('Reconstruction accuracy')
ylabel('Number of Events')
legend('Mortar', 'Distributed')

saveas(f4, pp.paper_foldername + "distributed_scatter" + ".png");
saveas(f5, pp.paper_foldername + "distributed_mortar_scatter" + ".png");
%
%% Plot 6: Inset good examples

% Load two datasets
fname = "distributed_inset1_max";
load(pp.intermediate_foldername + fname + ".mat")

title_str = "";
plot_reconstruction_helper

saveas(f1, pp.paper_foldername + fname + ".png");


% Load two datasets
% fname = "distributed_inset2_max_single_event";
% load(pp.intermediate_foldername + fname + ".mat")
% 
% title_str = "";
% plot_reconstruction_helper
% 
% saveas(f1, pp.paper_foldername + fname + ".png");


%% Figure 2 ALT: violin + examples

f1a = figure(fig_opt{:});
fname = "distributed_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;
subplot(121)
violinplot([acc_m], {"Distributed"})
ylabel("Reconstruction Accuracy")
ylim([-0.1, 0.6])
title("Variance Explained")

% Load all raw data
% fname = "distributed_all_raw_data";
% load(pp.intermediate_foldername + fname + ".mat")

% example_ind = [1822, 85, 7, 5, 19, 22];
example_ind = [85, 7, 5, 22];
num_rows = length(example_ind);
num_delays = 10; % TODO: not hardcode
cmap = colormap(parula(6));

for i = 1:num_rows
    i_subplot = 2*i;
    subplot(num_rows, 2, i_subplot);
    
    i_data = example_ind(i);
    
    X = all_X(i_data,:);
    X_recon = all_X_reconstructed(i_data,:);
    U = all_paths{i_data}.U;
    
    plot(ts, X, line_opt{:}, 'color', cmap(i,:))
    hold on
    plot(ts2, X_recon, line_opt{:})
    
    plot(ts3, 0.025*U - 0.05, 'k', line_opt{:})
    xlim([1,200])
    xticks([])
    yticks([])
    
    this_accuracy = acc_m(i_data);
    y = this_accuracy;
    
%     title_str = sprintf("%.0f %s", 100*this_accuracy, "%");
%     title(title_str)
    
    if i==1
        title("Example time series")
    end
    if i<num_rows
        xticks([])
    else
        xlabel("Time (\mu s)")
        legend({"Time series", "Reconstruction", "Control signal"})
    end
    
    subplot(121)
    plot([0.5, 1.5], [y, y], 'color', cmap(i, :), 'linewidth', 3)
end

%
out_fname = pp.paper_foldername + "fig2/" + "alternate_violin_distributed";
saveas(f1a, out_fname + ".png");


%% Figure 2 ALT ALT: cannot disentangle events using FFT
f1 = figure(fig_opt{:});

% freq_threshold = 0.1;

fname = "\spectrogram_mortar_fnames";
load(pp.intermediate_foldername + fname + ".mat")
freq_m = f(all_locs);
% freq_m = freq_m(freq_m>freq_threshold);

%
fname = "\spectrogram_distributed_fnames";
load(pp.intermediate_foldername + fname + ".mat")
freq_d = f(all_locs);
% freq_d = freq_d(freq_d>freq_threshold);

%
fname = "\spectrogram_localized_fnames";
load(pp.intermediate_foldername + fname + ".mat")
freq_l = f(all_locs);
% freq_l = freq_l(freq_l>freq_threshold);

to_plot = struct(...
    'Mortar', freq_m,...
    'Distributed', freq_d,...
    'Localized', freq_l...
);

subplot(121)
violinplot(to_plot);
hold on
ylabel("Hz")
% title("Dominant frequency")
ylim([9e4, 2.75e5])
xlim([0.5, 3.5])

num_cols = 6;
cmap = colormap(parula(6));

%
fname = "mortar_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
freq_targets = [2.4e5, 2.3e5, 2.2e5, 2e5];
example_ind = get_indices_at_y_level(freq_targets, freq_m);
inset_subplots = [4, 10, 16, 22];
x_for_violin_plot = 1;
plot_time_series_column(all_X, [], [],...
    example_ind, ts, ts2, ts3, freq_m, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Mortar")

%
fname = "distributed_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
freq_targets = [2.4e5, 2.2e5, 1.4e5, 1e5];
example_ind = get_indices_at_y_level(freq_targets, freq_d);
inset_subplots = [5, 11, 17, 23];
x_for_violin_plot = 2;
plot_time_series_column(all_X, [], [],...
    example_ind, ts, ts2, ts3, freq_d, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Distributed")

%
fname = "localized_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
freq_targets = [2.4e5, 2.2e5, 1.4e5, 1e5];
example_ind = get_indices_at_y_level(freq_targets, freq_l);
inset_subplots = [6, 12, 18, 24];
x_for_violin_plot = 3;
plot_time_series_column(all_X, [], [],...
    example_ind, ts, ts2, ts3, freq_l, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Localized")

%

set_times_new_roman
out_basename = pp.paper_foldername + "fig3\triple_frequencies";
saveas(f1, out_basename + ".png");
saveas(f1, out_basename + ".pdf");
% cleanfigure();
% matlab2tikz(char(out_basename + ".tex"), tikz_opt{:});
% set(gcf, 'color', 'none');
% export_fig(out_basename + ".pdf");
% export_fig(out_basename + ".png");
%% Localized
%%

%% Plot 4: Reconstruction accuracy
f3 = figure(fig_opt{:});

fname = "localized_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_l = accuracy;
fname = "distributed_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_d = accuracy;
fname = "mortar_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;

violinplot([acc_m, acc_d, acc_l],...
    {"Mortar", "Distributed", "Localized"})
ylabel("Reconstruction Accuracy")
title("Localized Events Are Poorly Reconstructed")

saveas(f3, pp.paper_foldername + "localized_accuracy" + ".png");
%% Plot 5: Scatterplot (with mortar as well)
% Load all three datasets
fname = "localized_scatter";
load(pp.intermediate_foldername + fname + ".mat")
accuracy_l = accuracy;
num_events_l = num_events;

fname = "distributed_scatter";
load(pp.intermediate_foldername + fname + ".mat")
accuracy_d = accuracy;
num_events_d = num_events;

fname = "mortar_scatter";
load(pp.intermediate_foldername + fname + ".mat")
accuracy_m = accuracy;
num_events_m = num_events;

% Finally, plot
f4 = figure(fig_opt{:});
scatter(accuracy_l, num_events_l, 'LineWidth',2, ...
    'markeredgecolor', [0.9290, 0.6940, 0.1250]	)
ylim([0,5])
xlim([0,1])
xlabel('Reconstruction accuracy')
ylabel('Number of Events')
title('Localized Events Usually Require Many Control Signals')

f5 = figure(fig_opt{:});
scatter(accuracy_m, num_events_m, '*', 'LineWidth',2)
hold on
scatter(accuracy_d, num_events_d, 'LineWidth',2)%,...
%     'markeredgecolor', 	[0.3010, 0.7450, 0.9330])
scatter(accuracy_l, num_events_l, 'LineWidth',2)%, ...
%     'markeredgecolor', [0.4940, 0.1840, 0.5560]	)
ylim([0,5])
xlim([0,1])
title('There Exist Mortar-Like Events within Localized')

xlabel('Reconstruction accuracy')
ylabel('Number of Events')
legend('Mortar', 'Distributed', 'Localized')

saveas(f4, pp.paper_foldername + "localized_scatter" + ".png");
saveas(f5, pp.paper_foldername + "localized_distributed_mortar_scatter" + ".png");
%
%% Plot 6: Inset

% Load two datasets
fname = "localized_inset1_good";
load(pp.intermediate_foldername + fname + ".mat")
title_str = "Mortar-like event in Distributed";
plot_reconstruction_helper

saveas(f1, pp.paper_foldername + fname + ".png");

%
fname = "localized_inset2_good";
load(pp.intermediate_foldername + fname + ".mat")

title_str = "Mortar-like event in Distributed";
plot_reconstruction_helper

saveas(f1, pp.paper_foldername + fname + ".png");
%% Plot 7: Multi-event insets

fname = "localized_inset1_multi";
load(pp.intermediate_foldername + fname + ".mat")
title_str = "Multiple Control Signals in Distributed";
plot_reconstruction_helper

saveas(f1, pp.paper_foldername + fname + ".png");

%
fname = "localized_inset2_multi";
load(pp.intermediate_foldername + fname + ".mat")
title_str = "Multiple Control Signals in Distributed";
plot_reconstruction_helper

saveas(f1, pp.paper_foldername + fname + ".png");


%% Figure 3 ALT: violin + examples

f1a = figure(fig_opt{:});
fname = "localized_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;
subplot(121)
violinplot([acc_m], {"Localized"})
ylabel("Reconstruction Accuracy")
ylim([-0.1, 0.6])
% title("Variance Explained")

example_ind = [10, 1, 2, 3];
num_rows = length(example_ind);
num_delays = 10; % TODO: not hardcode
cmap = colormap(parula(6));

for i = 1:num_rows
    i_subplot = 2*i;
    subplot(num_rows, 2, i_subplot);
    
    i_data = example_ind(i);
    
    X = all_X(i_data,:);
    X_recon = all_X_reconstructed(i_data,:);
    U = all_paths{i_data}.U;
    
    plot(ts, X, line_opt{:}, 'color', cmap(i,:))
    hold on
    plot(ts2, X_recon, line_opt{:})
    
    plot(ts3, 0.025*U - max(abs(X)), 'k', line_opt{:})
    xlim([1,200])
    ylim([-max(abs(X)), max(abs(X))])
    xticks([])
    yticks([])
    
    this_accuracy = acc_m(i_data);
    y = this_accuracy;
    
    if i==1
        title("Example time series")
    end
    if i<num_rows
        xticks([])
    else
        xlabel("Time (\mu s)")
        legend({"Time series", "Reconstruction", "Control signal"})
    end
    
    subplot(121)
    plot([0.5, 1.5], [y, y], 'color', cmap(i, :), 'linewidth', 3)
end

%
out_fname = pp.paper_foldername + "fig2/" + "alternate_violin_distributed";
saveas(f1a, out_fname + ".png");

%% Figure 3 ALT ALT: violin of all 3 datasets + examples

f1a = figure(fig_opt{:});
fname = "mortar_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;
fname = "distributed_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_d = accuracy;
fname = "localized_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_l = accuracy;

to_plot = struct(...
    'Mortar', acc_m,...
    'Distributed', acc_d,...
    'Localized', acc_l...
);

subplot(121)
violinplot(to_plot);
ylabel("Variance Explained")
ylim([-0.1, 0.6])
xlim([0.5, 3.5])

% Examples will be from each, in columns
cmap = colormap(parula(6));
num_cols = 6;

%
fname = "mortar_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;
acc_targets = [0.55, 0.4, 0.25, 0.1];
example_ind = get_indices_at_y_level(acc_targets, acc_m);
inset_subplots = [4, 10, 16, 22];
x_for_violin_plot = 1;
plot_time_series_column(all_X, all_X_reconstructed, all_paths,...
    example_ind, ts, ts2, ts3, acc_m, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Mortar")

%
fname = "distributed_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_targets = [0.4, 0.3, 0.2, 0.1];
example_ind = get_indices_at_y_level(acc_targets, acc_d);
inset_subplots = [5, 11, 17, 23];
x_for_violin_plot = 2;
plot_time_series_column(all_X, all_X_reconstructed, all_paths,...
    example_ind, ts, ts2, ts3, acc_d, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Distributed")
%
fname = "localized_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_l = accuracy;
acc_targets = [0.55, 0.4, 0.25, 0.1];
example_ind = get_indices_at_y_level(acc_targets, acc_l);
inset_subplots = [6, 12, 18, 24];
x_for_violin_plot = 3;

plot_time_series_column(all_X, all_X_reconstructed, all_paths,...
    example_ind, ts, ts2, ts3, acc_l, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    true, "Localized")

set_times_new_roman
%
out_fname = pp.paper_foldername + "fig2/" + "alternate_violin_distributed";
saveas(f1a, out_fname + ".png");
saveas(f1a, out_fname + ".pdf");


%% Figure 4: Explained variance + absolute variance

f4 = figure(fig_opt{:});
% Load accuracy
fname = "mortar_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;
% num_m = num_events;
% num_m = num_m(1:250);
fname = "mortar_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
abs_m = vecnorm(all_X');

fname = "distributed_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_d = accuracy;
% acc_d = acc_d(1:250);
% num_d = num_events;
% num_d = num_d(1:250);
fname = "distributed_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
abs_d = vecnorm(all_X');
% abs_d = abs_d(1:250);

fname = "localized_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_l = accuracy;
% num_l = num_events;
fname = "localized_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
abs_l = vecnorm(all_X');

% Load FFT
% freq_threshold = 0.1;
freq_threshold = -inf;

fname = "\spectrogram_mortar_fnames";
load(pp.intermediate_foldername + fname + ".mat")
freq_m = f(all_locs);
to_keep_m = freq_m>freq_threshold;
freq_m = freq_m(to_keep_m);

%
fname = "\spectrogram_distributed_fnames";
load(pp.intermediate_foldername + fname + ".mat")
freq_d = f(all_locs);
to_keep_d = freq_d>freq_threshold;
freq_d = freq_d(to_keep_d);

%
fname = "\spectrogram_localized_fnames";
load(pp.intermediate_foldername + fname + ".mat")
freq_l = f(all_locs);
to_keep_l = freq_l>freq_threshold;
freq_l = freq_l(to_keep_l);
% Scratch Plot: accuracy vs absolute amplitude
thresh = 0.4;

plot(num_m(acc_m>thresh), acc_m(acc_m>thresh), 'o')
hold on
plot(num_d(acc_d>thresh), acc_d(acc_d>thresh), 'o')
plot(num_l(acc_l>thresh), acc_l(acc_l>thresh), '*')
% plot(acc_m, abs_m, 'o')
% hold on
% plot(acc_d, abs_d, 'o')
% plot(acc_l, abs_l, '*')
xlabel("Num")
ylabel("acc")

legend({"M", "D", "L"})
%% Actual Plot: frequency vs. accuracy
% f1 = figure();
thresh = 0.0;

% plot(freq_m, acc_m(to_keep_m), 'o')
% hold on
% plot(freq_d, acc_d(to_keep_d), 'o')
% plot(freq_l, acc_l(to_keep_l), '*')
% xlabel("Frequency")
% ylabel("Reconstruction accuracy")
% legend({"Mortar", "Distributed", "Localized"})
% title("Accuracy is not correlated to peak frequency")
% 
% out_fname = pp.paper_foldername + "meeting_10_12/" + "freq_vs_acc";
% saveas(f1, out_fname + ".png");

% plot(freq_m(acc_m>thresh), acc_m(acc_m>thresh), 'o')
% hold on
% plot(freq_d(acc_d>thresh), acc_d(acc_d>thresh), 'o')
% plot(freq_l(acc_l>thresh), acc_l(acc_l>thresh), '*')
%
f4 = figure(fig_opt{:});

subplot(1,2,1)
x_sz = 50;
freq_hist_bins = linspace(9e4, 2.75e5, x_sz);
x_hist_bins = linspace(-0.1, 0.6, x_sz);
[x_freq, x_acc] = meshgrid(freq_hist_bins, x_hist_bins);
x1_vec = x_freq(:);
x2_vec = x_acc(:);
x_support = [x1_vec, x2_vec];

dens = ksdensity([freq_m, acc_m(to_keep_m)], x_support);
contour(x_freq, x_acc, reshape(dens, [x_sz,x_sz]), 'r', 'linewidth', 1.5)
hold on
dens = ksdensity([freq_d, acc_d(to_keep_d)], x_support);
contour(x_freq, x_acc, reshape(dens, [x_sz,x_sz]), 'b', 'linewidth', 1.5)
dens = ksdensity([freq_l, acc_l(to_keep_l)], x_support);
contour(x_freq, x_acc, reshape(dens, [x_sz,x_sz]), 'g', 'linewidth', 1.5)

% ax = gca;
% labels = string(ax.XAxis.TickLabels); % extract
% labels(2:2:end) = ' '; % remove every other one
% ax.XAxis.TickLabels = labels; % set

xlabel("Hz")
ylabel("Variance Explained")

%
fname = "mortar_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_targets = [0.5, 0.4, 0.2, 0.1];
example_ind = get_indices_at_y_level(acc_targets, acc_m);
inset_subplots = [4, 10, 16, 22];
x_for_violin_plot = freq_m(example_ind);
plot_time_series_column(all_X, [], [],...
    example_ind, ts, ts2, ts3, acc_m, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Mortar", false, 's')
%

%
fname = "distributed_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_targets = [0.5, 0.3, 0.1, 0.0];
example_ind = get_indices_at_y_level(acc_targets, acc_d);
inset_subplots = [5, 11, 17, 23];
x_for_violin_plot = freq_d(example_ind);
plot_time_series_column(all_X, [], [],...
    example_ind, ts, ts2, ts3, acc_d, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Distributed", false, 'o')

%
fname = "localized_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_targets = [0.5, 0.3, 0.1, 0.0];
example_ind = get_indices_at_y_level(acc_targets, acc_l);
inset_subplots = [6, 12, 18, 24];
x_for_violin_plot = freq_l(example_ind);
plot_time_series_column(all_X, [], [],...
    example_ind, ts, ts2, ts3, acc_l, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Localized", false, 'v')

leg = legend({"Mortar", "Distributed", "Localized"});
leg.String(4:end) = [];


set_times_new_roman

out_fname = pp.paper_foldername + "fig4/" + "freq_vs_acc_contour";
saveas(f4, out_fname + ".png");
saveas(f4, out_fname + ".pdf");
