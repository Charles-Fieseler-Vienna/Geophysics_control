
f1 = figure(fig_opt{:});

subplot(211)
plot(ts, X, 'Linewidth', 2)
hold on
plot(ts2, X_recon, 'LineWidth',2)
title(title_str)
legend({'Data', 'Reconstruction'}, 'Location', 'Southeast')
xlim([1,200])
ylabel('Amplitude (Volts)')

subplot(212)
plot(ts3, U, 'LineWidth',2)
legend('Learned Control Signal')
xlim([1,200])
xlabel('Time (micro seconds)')
ylabel('Amplitude (A.U.)')