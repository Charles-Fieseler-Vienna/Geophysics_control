%% Shared setup
%%
pp = PurdueProject();

offsets = normrnd(0, 0.05, 250, 1);
% Time delay embedding shortens the time series, so this should align the
% control signals and the original (un-embedded) data
len = 499;
ts = 0.5*(1:len);
ts2 = 0.5*(1:len-7);
ts3 = 0.5*(1:len-8);

cmap = colormap(parula(6));
num_cols = 6;
% subset_ind = 1:5:2500;

% Options for display and saving
tikz_opt = {'width', '0.4\textwidth', 'height', '0.2\textheight'};
fig_opt = {'DefaultAxesFontSize', 20, 'WindowState', 'Maximized'};
line_opt = {'LineWidth',1};

actually_save = false;

error("The script is initialized. Next, please run individual sections")

%% Figure ?: accuracy 2 datasets + examples
len = 499;
ts = 0.5*(1:len);
ts2 = 0.5*(1:len-7);
ts3 = 0.5*(1:len-8);

figure6 = figure(fig_opt{:});
fname = "unknown_fnames1_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;
fname = "unknown_fnames2_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_d = accuracy;

acc_no_outliers = acc_d;
acc_no_outliers = acc_no_outliers(acc_d>0);
to_plot = struct(...
    'Unknown_1', acc_m,...
    'Unknown_2', acc_no_outliers...
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
fname = "unknown_fnames1_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_m = accuracy;
% acc_targets = [0.0, 0.05, 0.1, 0.2];
acc_targets = [0.2, 0.1, 0.05, 0.0];
example_ind_mortar = get_indices_at_y_level(acc_targets, acc_m);
inset_subplots = [4, 10, 16, 22];

all_X = all_X(:, 1:length(ts));

x_for_violin_plot = 1;
plot_time_series_column(all_X, all_X_reconstructed, all_paths,...
    example_ind_mortar, ts, ts2, ts3, acc_m, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Unknown_1")

%
fname = "unknown_fnames2_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
acc_targets = [0.2, 0.1, 0.05, 0.0];

% acc_d = acc_d(1:end-1);

example_ind_distributed = get_indices_at_y_level(acc_targets, acc_d);
inset_subplots = [5, 11, 17, 23];
x_for_violin_plot = 2;

all_X = all_X(:, 1:length(ts));

plot_time_series_column(all_X, all_X_reconstructed, all_paths,...
    example_ind_distributed, ts, ts2, ts3, acc_d, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Unknown_2")

set_times_new_roman
%

if actually_save
    out_fname = pp.paper_foldername + "fig_unknown/" + "fig_unknown_acc";
    saveas(figure6, out_fname + ".png");
    saveas(figure6, out_fname + ".pdf");
    
    % Also save some metadata
    out_fname = pp.paper_foldername + "../intermediate_raw/" + "fig_unknown.mat";
    save(out_fname,...
        'example_ind_distributed',...
        'example_ind_mortar')
end

%% Figure ?: FFT of unknown
figure7 = figure(fig_opt{:});

len = 500;
ts = 0.5*(1:len);
ts2 = 0.5*(1:len-7);
ts3 = 0.5*(1:len-8);

% Optional threshold to remove outliers
% freq_threshold = 0.1;

fname = "\spectrogram_unknown_fnames1";
load(pp.intermediate_foldername + fname + ".mat")
freq_1 = f(all_locs);
% freq_m = freq_m(freq_m>freq_threshold);

%
fname = "\spectrogram_unknown_fnames2";
load(pp.intermediate_foldername + fname + ".mat")
freq_2 = f(all_locs);
% freq_d = freq_d(freq_d>freq_threshold);

%

to_plot = struct(...
    'Unknown1', freq_1,...
    'Unknown2', freq_2...
);

subplot(121)
violinplot(to_plot);
hold on
ylabel("Hz")
% title("Dominant frequency")
% ylim([9e4, 2.75e5])
xlim([0.5, 3.5])

num_cols = 6;
cmap = colormap(parula(6));

%
fname = "unknown_fnames1_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
freq_targets = [2.4e5, 2.3e5, 2.2e5, 2e5];
example_ind_mortar = get_indices_at_y_level(freq_targets, freq_1);
inset_subplots = [4, 10, 16, 22];
x_for_violin_plot = 1;

% all_X = all_X(:, subset_ind);

plot_time_series_column(all_X, [], [],...
    example_ind_mortar, ts, ts2, ts3, freq_1, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Unknown1")

%
fname = "unknown_fnames2_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
freq_targets = [2.4e5, 2.2e5, 1.4e5, 1e5];
example_ind_distributed = get_indices_at_y_level(freq_targets, freq_2);
inset_subplots = [5, 11, 17, 23];
x_for_violin_plot = 2;

plot_time_series_column(all_X, [], [],...
    example_ind_distributed, ts, ts2, ts3, freq_2, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Unknown2")

%
set_times_new_roman

if actually_save
    out_basename = pp.paper_foldername + "fig_unknown\fig_unknown_fft";
    saveas(figure7, out_basename + ".png");
    saveas(figure7, out_basename + ".pdf");
    
    % Also save some metadata
    out_fname = pp.paper_foldername + "../intermediate_raw/" + "fig_unknown_fft.mat";
    save(out_fname,...
        'example_ind_distributed',...
        'example_ind_mortar')
end

%% Figure ?: Contour plot


figure8 = figure(fig_opt{:});
%
% Load accuracy
%
fname = "unknown_fnames1_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_1 = accuracy;
fname = "unknown_fnames1_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
abs_1 = vecnorm(all_X');

fname = "unknown_fnames2_acc";
load(pp.intermediate_foldername + fname + ".mat")
acc_2 = accuracy;
fname = "unknown_fnames2_all_raw_data";
load(pp.intermediate_foldername + fname + ".mat")
abs_2 = vecnorm(all_X');

%
% Load FFT, and optionally apply a threshold (to keep the axes sensible)
%
% freq_threshold = 0.1;
freq_threshold = -inf;

fname = "\spectrogram_unknown_fnames1";
load(pp.intermediate_foldername + fname + ".mat")
freq_1 = f(all_locs);
to_keep_m = freq_1>freq_threshold;
freq_1 = freq_1(to_keep_m);

%
fname = "\spectrogram_unknown_fnames2";
load(pp.intermediate_foldername + fname + ".mat")
freq_2 = f(all_locs);
to_keep_d = freq_2>freq_threshold;
freq_2 = freq_2(to_keep_d);

%
% Actually set up the plot
%

% Set up the histograms for the contour (density) plot
subplot(1,2,1)
x_sz = 50;
% freq_hist_bins = linspace(9e4, 2.75e5, x_sz);
freq_hist_bins = linspace(0, 4e5, x_sz);
x_hist_bins = linspace(-0.1, 0.6, x_sz);
[x_freq, x_acc] = meshgrid(freq_hist_bins, x_hist_bins);
x1_vec = x_freq(:);
x2_vec = x_acc(:);
x_support = [x1_vec, x2_vec];

dens = ksdensity([freq_1, acc_1(to_keep_m)], x_support);
contour(x_freq, x_acc, reshape(dens, [x_sz,x_sz]), 'r', 'linewidth', 1.5)
hold on
dens = ksdensity([freq_2, acc_2(to_keep_d)], x_support);
contour(x_freq, x_acc, reshape(dens, [x_sz,x_sz]), 'b', 'linewidth', 1.5)

% ax = gca;
% labels = string(ax.XAxis.TickLabels); % extract
% labels(2:2:end) = ' '; % remove every other one
% ax.XAxis.TickLabels = labels; % set

%
% Labels etc. for plot
%
xlabel("Hz")
ylabel("Variance Explained")

%
fname = "\spectrogram_unknown_fnames1";
load(pp.intermediate_foldername + fname + ".mat")
acc_targets = [0.1, 0.05, 0.01, 0.0];
example_ind_mortar = get_indices_at_y_level(acc_targets, acc_1);
inset_subplots = [4, 10, 16, 22];
x_for_violin_plot = freq_1(example_ind_mortar);
plot_time_series_column(all_X(:, 1:end-1), [], [],...
    example_ind_mortar, ts, ts2, ts3, acc_1, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Unknown1", false, 's')
%
fname = "\spectrogram_unknown_fnames2";
load(pp.intermediate_foldername + fname + ".mat")
acc_targets = [0.1, 0.05, 0.01, 0.0];
example_ind_distributed = get_indices_at_y_level(acc_targets, acc_2);
inset_subplots = [5, 11, 17, 23];
x_for_violin_plot = freq_2(example_ind_distributed);
plot_time_series_column(all_X(:, 1:end-1), [], [],...
    example_ind_distributed, ts, ts2, ts3, acc_2, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    false, "Distributed", false, 'o')

leg = legend({"Unknown 1", "Unknown 1"});
leg.String(4:end) = [];

set_times_new_roman
if actually_save
    out_fname = pp.paper_foldername + "fig_unknown/" + "freq_vs_acc_contour";
    saveas(figure8, out_fname + ".png");
    saveas(figure8, out_fname + ".pdf");
    
    % Also save some metadata
    out_fname = pp.paper_foldername + "../intermediate_raw/" + "fig_unknown_contour_indices.mat";
    save(out_fname,...
        'example_ind_localized',...
        'example_ind_distributed',...
        'example_ind_mortar')
end