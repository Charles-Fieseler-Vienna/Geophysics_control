%% Shared setup
%%
pp = ControlProject();

offsets = normrnd(0, 0.05, 250, 1);
% Time delay embedding shortens the time series, so this should align the
% control signals and the original (un-embedded) data
ts = 0.5*(1:401);
ts2 = 0.5*(1:393);
ts3 = 0.5*(1:392);

% Options for display and saving
tikz_opt = {'width', '0.4\textwidth', 'height', '0.2\textheight'};
fig_opt = {'DefaultAxesFontSize', 20, 'WindowState', 'Maximized'};
line_opt = {'LineWidth',1};

actually_save = false;

error("The script is initialized. Next, please run individual sections")

%% Figure 6: violin of all 3 datasets + examples

figure6 = figure(fig_opt{:});
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
% Note; the accuracies should be reloaded to keep the indexing aligned
%
fname = "mortar_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;
acc_targets = [0.55, 0.4, 0.25, 0.1];
example_ind_mortar = get_indices_at_y_level(acc_targets, acc_m);
inset_subplots = [4, 10, 16, 22];
x_for_violin_plot = 1;
plot_time_series_column(all_X, all_X_reconstructed, all_paths,...
    example_ind_mortar, ts, ts2, ts3, acc_m, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Mortar")

%
fname = "distributed_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_targets = [0.4, 0.3, 0.2, 0.1];
example_ind_distributed = get_indices_at_y_level(acc_targets, acc_d);
inset_subplots = [5, 11, 17, 23];
x_for_violin_plot = 2;
plot_time_series_column(all_X, all_X_reconstructed, all_paths,...
    example_ind_distributed, ts, ts2, ts3, acc_d, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Distributed")
%
fname = "localized_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_l = accuracy;
acc_targets = [0.55, 0.4, 0.25, 0.1];
example_ind_localized = get_indices_at_y_level(acc_targets, acc_l);
inset_subplots = [6, 12, 18, 24];
x_for_violin_plot = 3;

plot_time_series_column(all_X, all_X_reconstructed, all_paths,...
    example_ind_localized, ts, ts2, ts3, acc_l, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    true, "Localized")

set_times_new_roman
%

if actually_save
    out_fname = pp.paper_foldername + "fig2/" + "alternate_violin_distributed";
    saveas(figure6, out_fname + ".png");
    saveas(figure6, out_fname + ".pdf");
    
    % Also save some metadata
    out_fname = pp.paper_foldername + "../intermediate_raw/" + "fig6_indices.mat";
    save(out_fname,...
        'example_ind_localized',...
        'example_ind_distributed',...
        'example_ind_mortar')
end

%% Figure 7: cannot disentangle events using FFT
figure7 = figure(fig_opt{:});

% Optional threshold to remove outliers
% freq_threshold = 0.1;
% Convert to kHz

fname = "\spectrogram_mortar_fnames";
load(pp.intermediate_foldername + fname + ".mat")
freq_m = f(all_locs) / 1000;
% freq_m = freq_m(freq_m>freq_threshold);

%
fname = "\spectrogram_distributed_fnames";
load(pp.intermediate_foldername + fname + ".mat")
freq_d = f(all_locs) / 1000;
% freq_d = freq_d(freq_d>freq_threshold);

%
fname = "\spectrogram_localized_fnames";
load(pp.intermediate_foldername + fname + ".mat")
freq_l = f(all_locs) / 1000;
% freq_l = freq_l(freq_l>freq_threshold);

to_plot = struct(...
    'Mortar', freq_m,...
    'Distributed', freq_d,...
    'Localized', freq_l...
);

subplot(121)
violinplot(to_plot);
hold on
ylabel("kHz")
% title("Dominant frequency")
ylim([9e1, 2.75e2])
xlim([0.5, 3.5])

num_cols = 6;
cmap = colormap(parula(6));

%
fname = "mortar_fnames_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
freq_targets = [2.4e2, 2.3e2, 2.2e2, 2e2];
example_ind_mortar = get_indices_at_y_level(freq_targets, freq_m);
inset_subplots = [4, 10, 16, 22];
x_for_violin_plot = 1;
plot_time_series_column(all_X, [], [],...
    example_ind_mortar, ts, ts2, ts3, freq_m, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Mortar")

%
fname = "distributed_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
freq_targets = [2.4e2, 2.2e2, 1.4e2, 1e2];
example_ind_distributed = get_indices_at_y_level(freq_targets, freq_d);
inset_subplots = [5, 11, 17, 23];
x_for_violin_plot = 2;
plot_time_series_column(all_X, [], [],...
    example_ind_distributed, ts, ts2, ts3, freq_d, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Distributed")

%
fname = "localized_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
freq_targets = [2.4e2, 2.2e2, 1.4e2, 1e2];
example_ind_localized = get_indices_at_y_level(freq_targets, freq_l);
inset_subplots = [6, 12, 18, 24];
x_for_violin_plot = 3;
plot_time_series_column(all_X, [], [],...
    example_ind_localized, ts, ts2, ts3, freq_l, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Localized")

%
set_times_new_roman

if actually_save
    out_basename = pp.paper_foldername + "fig3\triple_frequencies";
    saveas(figure7, out_basename + ".png");
    saveas(figure7, out_basename + ".pdf");
    
    % Requires external package
    % cleanfigure();
    % matlab2tikz(char(out_basename + ".tex"), tikz_opt{:});
    % set(gcf, 'color', 'none');
    % export_fig(out_basename + ".pdf");
    % export_fig(out_basename + ".png");
    
    % Also save some metadata
    out_fname = pp.paper_foldername + "../intermediate_raw/" + "fig7_indices.mat";
    save(out_fname,...
        'example_ind_localized',...
        'example_ind_distributed',...
        'example_ind_mortar')
end

%% Figure 8: Contour plot of Explained variance + FFT frequency

figure8 = figure(fig_opt{:});
%
% Load accuracy
%
fname = "mortar_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;
fname = "mortar_fnames_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
abs_m = vecnorm(all_X');

fname = "distributed_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_d = accuracy;
fname = "distributed_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
abs_d = vecnorm(all_X');

fname = "localized_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_l = accuracy;
fname = "localized_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
abs_l = vecnorm(all_X');

%
% Load FFT, and optionally apply a threshold (to keep the axes sensible)
%
% freq_threshold = 0.1;
freq_threshold = -inf;

fname = "\spectrogram_mortar_fnames";
load(pp.intermediate_foldername + fname + ".mat")
freq_m = f(all_locs);
to_keep_m = freq_m>freq_threshold;
freq_m = freq_m(to_keep_m) / 1000;

%
fname = "\spectrogram_distributed_fnames";
load(pp.intermediate_foldername + fname + ".mat")
freq_d = f(all_locs);
to_keep_d = freq_d>freq_threshold;
freq_d = freq_d(to_keep_d) / 1000;

%
fname = "\spectrogram_localized_fnames";
load(pp.intermediate_foldername + fname + ".mat")
freq_l = f(all_locs);
to_keep_l = freq_l>freq_threshold;
freq_l = freq_l(to_keep_l) / 1000;

%
% Actually set up the plot
%

% Set up the histograms for the contour (density) plot
subplot(1,2,1)
x_sz = 50;
freq_hist_bins = linspace(9e1, 2.75e2, x_sz);
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

%
% Labels etc. for plot
%
xlabel("kHz")
ylabel("Variance Explained")

%
fname = "mortar_fnames_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_targets = [0.5, 0.4, 0.2, 0.1];
example_ind_mortar = get_indices_at_y_level(acc_targets, acc_m);
inset_subplots = [4, 10, 16, 22];
x_for_violin_plot = freq_m(example_ind_mortar);
plot_time_series_column(all_X, [], [],...
    example_ind_mortar, ts, ts2, ts3, acc_m, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Mortar", false, 's')
%
fname = "distributed_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_targets = [0.5, 0.3, 0.1, 0.0];
example_ind_distributed = get_indices_at_y_level(acc_targets, acc_d);
inset_subplots = [5, 11, 17, 23];
x_for_violin_plot = freq_d(example_ind_distributed);
plot_time_series_column(all_X, [], [],...
    example_ind_distributed, ts, ts2, ts3, acc_d, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Distributed", false, 'o')
%
fname = "localized_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_targets = [0.5, 0.3, 0.1, 0.0];
example_ind_localized = get_indices_at_y_level(acc_targets, acc_l);
inset_subplots = [6, 12, 18, 24];
x_for_violin_plot = freq_l(example_ind_localized);
plot_time_series_column(all_X, [], [],...
    example_ind_localized, ts, ts2, ts3, acc_l, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Localized", false, 'v')

leg = legend({"Mortar", "Distributed", "Localized"}, ...
    'Location', 'northwest');
leg.String(4:end) = [];

set_times_new_roman
if actually_save
    out_fname = pp.paper_foldername + "fig4/" + "freq_vs_acc_contour";
    saveas(figure8, out_fname + ".png");
    saveas(figure8, out_fname + ".pdf");
    
    % Also save some metadata
    out_fname = pp.paper_foldername + "../intermediate_raw/" + "fig8_indices.mat";
    save(out_fname,...
        'example_ind_localized',...
        'example_ind_distributed',...
        'example_ind_mortar')
end

%% Figure 9?: violin of 3 new mortar datasets + original + examples

include_original = false;

len = 500;
ts = 0.5*(1:len);
ts2 = 0.5*(1:len-8);
ts3 = 0.5*(1:len-9);

figure6 = figure(fig_opt{:});
fname = "mortar_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_original = accuracy;
fname = "mortar_new1_fnames_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;
fname = "mortar_new4_fnames_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_d = accuracy;
fname = "mortar_new3_fnames_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_l = accuracy;

% There are outliers in some
to_keep_d = abs(acc_d) < 1;
acc_d = acc_d(to_keep_d);

if include_original
    to_plot = struct(...
        'Original', acc_original,...
        'Mortar1', acc_m,...
        'Mortar2', acc_d,...
        'Mortar3', acc_l...
    );
else
    to_plot = struct(...
        'Mortar1', acc_m,...
        'Mortar2', acc_d,...
        'Mortar3', acc_l...
    );
end

subplot(121)
violinplot(to_plot);
ylabel("Variance Explained")
ylim([-0.1, 0.6])
if include_original
    xlim([0.5, 4.5])
else
    xlim([0.5, 3.5])
end

% Examples will be from each, in columns
cmap = colormap(parula(6));
if include_original
    num_cols = 8;
else
    num_cols = 6;
end

%
% Note; the accuracies should be reloaded to keep the indexing aligned
%
if include_original
    len = 401;
    ts_raw = 0.5*(1:len);
    ts2_raw = 0.5*(1:len-8);
    ts3_raw = 0.5*(1:len-9);
    fname = "mortar_fnames_all_raw_data";
    load(pp.intermediate_foldername + fname + ".mat")
    acc_m = accuracy;
    acc_targets = [0.55, 0.4, 0.25, 0.1];
    example_ind_original = get_indices_at_y_level(acc_targets, acc_m);
    inset_subplots = [5, 13, 21, 29];
    x_for_violin_plot = 1;
    plot_time_series_column(all_X, all_X_reconstructed, all_paths,...
        example_ind_original, ts_raw, ts2_raw, ts3_raw, acc_m, num_cols,...
        inset_subplots, line_opt, cmap, x_for_violin_plot,...
        false, "Mortar Original")
end
%
fname = "mortar_new1_fnames_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;
acc_targets = [0.3, 0.25, 0.2, 0.1];
example_ind_mortar1 = get_indices_at_y_level(acc_targets, acc_m);
if include_original
    inset_subplots = [6, 14, 22, 30];
    x_for_violin_plot = 2;
else
    inset_subplots = [4, 10, 16, 22];
    x_for_violin_plot = 1;
end
plot_time_series_column(all_X, all_X_reconstructed, all_paths,...
    example_ind_mortar1, ts, ts2, ts3, acc_m, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Mortar1")

%
fname = "mortar_new2_fnames_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_targets = [0.4, 0.3, 0.2, 0.1];
% For some reason the violin plot won't work with all
all_X = all_X(to_keep_d,:);
all_X_reconstructed = all_X_reconstructed(to_keep_d,:);
all_paths = all_paths(to_keep_d);
example_ind_mortar2 = get_indices_at_y_level(acc_targets, acc_d);
if include_original
    inset_subplots = [7, 15, 23, 31];
    x_for_violin_plot = 3;
else
    inset_subplots = [5, 11, 17, 23];
    x_for_violin_plot = 2;
end
plot_time_series_column(all_X, all_X_reconstructed, all_paths,...
    example_ind_mortar2, ts, ts2, ts3, acc_d, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Mortar2")
%
fname = "mortar_new3_fnames_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_l = accuracy;
acc_targets = [0.55, 0.4, 0.25, 0.1];
example_ind_mortar3 = get_indices_at_y_level(acc_targets, acc_l);
if include_original
    inset_subplots = [8, 16, 24, 32];
    x_for_violin_plot = 4;
else
    inset_subplots = [6, 12, 18, 24];
    x_for_violin_plot = 3;
end
plot_time_series_column(all_X, all_X_reconstructed, all_paths,...
    example_ind_mortar3, ts, ts2, ts3, acc_l, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    true, "Mortar3")

set_times_new_roman
%

if actually_save
    out_fname = pp.paper_foldername + "fig9/" + "alternate_violin_distributed";
    saveas(figure6, out_fname + ".png");
    saveas(figure6, out_fname + ".pdf");
    
    % Also save some metadata
    out_fname = pp.paper_foldername + "../intermediate_raw/" + "fig9_indices.mat";
    save(out_fname,...
        'example_ind_mortar3',...
        'example_ind_mortar2',...
        'example_ind_mortar1',...
        'example_ind_original')
end
