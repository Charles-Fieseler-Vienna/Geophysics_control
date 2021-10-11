%% Shared setup
%%
pp = PurdueProject();

fig_opt = {'DefaultAxesFontSize', 24, 'WindowState', 'Maximized'};

offsets = normrnd(0, 0.05, 250, 1);
ts = 0.5*(1:401);
ts2 = 0.5*(1:393);
ts3 = 0.5*(1:392);

error("Please run individual sections")

%% Intro
%%

%% Plot 1: 3 time series
fname = "mortar_ex1";
load(pp.intermediate_foldername + fname + ".mat")

f1 = figure(fig_opt{:});
plot(ts, X, 'Linewidth', 2)
title("Example Mortar (No Clay) Event")
xlim([1,200])
ylabel('Amplitude (Volts)')
xlabel('Time (micro seconds)')

saveas(f1, pp.presentation_foldername + "mortar_intro" + ".png");

%
fname = "distributed_ex1";
load(pp.intermediate_foldername + fname + ".mat")

f1 = figure(fig_opt{:});
plot(ts, X, 'Linewidth', 2)
title("Example Distributed Clay Event")
xlim([1,200])
ylabel('Amplitude (Volts)')
xlabel('Time (micro seconds)')

saveas(f1, pp.presentation_foldername + "distributed_intro" + ".png");

%
fname = "localized_inset1_multi";
load(pp.intermediate_foldername + fname + ".mat")

f1 = figure(fig_opt{:});
plot(ts, X, 'Linewidth', 2)
title("Example Localized Clay Event")
xlim([1,200])
ylabel('Amplitude (Volts)')
xlabel('Time (micro seconds)')

% saveas(f1, pp.presentation_foldername + "localized_intro" + ".png");

%% Mortar
%%

%% Plot 1: Example good dataset
f1 = figure(fig_opt{:});

fname = "mortar_ex1";
load(pp.intermediate_foldername + fname + ".mat")

subplot(211)
plot(ts, X, 'Linewidth', 2)
hold on
plot(ts2, X_recon, 'LineWidth',2)
title("Example Mortar Event (Successful Reconstruction)")
legend({'Data', 'Reconstruction'}, 'Location', 'Southeast')
xlim([1,200])
ylabel('Amplitude (Volts)')

subplot(212)
plot(ts3, U, 'LineWidth',2)
legend('Learned Control Signal')
xlim([1,200])
xlabel('Time (micro seconds)')
ylabel('Amplitude (A.U.)')

saveas(f1, pp.presentation_foldername + fname + ".png");
%
%% Plot 2: Example good dataset
f2 = figure(fig_opt{:});

fname = "mortar_ex2";
load(pp.intermediate_foldername + fname + ".mat")

subplot(211)
plot(ts, X, 'Linewidth', 2)
hold on
plot(ts2, X_recon, 'LineWidth',2)
title("Example Mortar Event (Successful Reconstruction)")
legend({'Data', 'Reconstruction'}, 'Location', 'Southeast')
xlim([1,200])
ylabel('Amplitude (Volts)')

subplot(212)
plot(ts3, U, 'LineWidth',2)
legend('Learned Control Signal')
xlim([1,200])
xlabel('Time (micro seconds)')
ylabel('Amplitude (A.U.)')

saveas(f2, pp.presentation_foldername + fname + ".png");
%% Plot 3: Reconstruction accuracy
f3 = figure(fig_opt{:});

fname = "mortar_acc";
load(pp.intermediate_foldername + fname + ".mat")

violinplot(accuracy)
title("Mortar Events Are Consistently Reconstructed")
ylabel("Reconstruction Accuracy")
xticklabels([])

saveas(f3, pp.presentation_foldername + fname + ".png");
%% Plot 4: Scatterplot
f4 = figure(fig_opt{:});

fname = "mortar_scatter";
load(pp.intermediate_foldername + fname + ".mat")

scatter(accuracy, num_events, '*', 'LineWidth',2)
xlabel('Reconstruction Accuracy')
ylabel('Number of Control Signals')
title("Mortar Reconstructions Require Few Control Signals")
ylim([0,5])
xlim([0,1])

saveas(f4, pp.presentation_foldername + fname + ".png");


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

saveas(f1, pp.presentation_foldername + fname + ".png");
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

saveas(f2, pp.presentation_foldername + fname + ".png");
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

saveas(f3, pp.presentation_foldername + fname + ".png");
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

saveas(f3, pp.presentation_foldername + "distributed_acc" + ".png");
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

saveas(f4, pp.presentation_foldername + "distributed_scatter" + ".png");
saveas(f5, pp.presentation_foldername + "distributed_mortar_scatter" + ".png");
%
%% Plot 6: Inset good examples

% Load two datasets
fname = "distributed_inset1_max";
load(pp.intermediate_foldername + fname + ".mat")

title_str = "";
plot_reconstruction_helper

saveas(f1, pp.presentation_foldername + fname + ".png");


% Load two datasets
% fname = "distributed_inset2_max_single_event";
% load(pp.intermediate_foldername + fname + ".mat")
% 
% title_str = "";
% plot_reconstruction_helper
% 
% saveas(f1, pp.presentation_foldername + fname + ".png");

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

saveas(f3, pp.presentation_foldername + "localized_accuracy" + ".png");
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

saveas(f4, pp.presentation_foldername + "localized_scatter" + ".png");
saveas(f5, pp.presentation_foldername + "localized_distributed_mortar_scatter" + ".png");
%
%% Plot 6: Inset

% Load two datasets
fname = "localized_inset1_good";
load(pp.intermediate_foldername + fname + ".mat")
title_str = "Mortar-like event in Distributed";
plot_reconstruction_helper

saveas(f1, pp.presentation_foldername + fname + ".png");

%
fname = "localized_inset2_good";
load(pp.intermediate_foldername + fname + ".mat")

title_str = "Mortar-like event in Distributed";
plot_reconstruction_helper

saveas(f1, pp.presentation_foldername + fname + ".png");
%% Plot 7: Multi-event insets

fname = "localized_inset1_multi";
load(pp.intermediate_foldername + fname + ".mat")
title_str = "Multiple Control Signals in Distributed";
plot_reconstruction_helper

saveas(f1, pp.presentation_foldername + fname + ".png");

%
fname = "localized_inset2_multi";
load(pp.intermediate_foldername + fname + ".mat")
title_str = "Multiple Control Signals in Distributed";
plot_reconstruction_helper

saveas(f1, pp.presentation_foldername + fname + ".png");
