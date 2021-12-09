function plot_time_series_column(all_X, all_X_reconstructed, all_paths,...
    example_ind, ts, ts2, ts3, acc, num_cols,...
    inset_subplots, line_opt, cmap, x_for_violin_plot,...
    to_plot_legend, top_title_str, violin_not_scatter, marker)
% Helper for plotting a column of example time series with controls
if ~exist('violin_not_scatter', 'var')
    violin_not_scatter = true;
end

if isempty(all_X_reconstructed)
    use_reconstructed = false;
else
    use_reconstructed = true;
end
if violin_not_scatter
    violin_x = [x_for_violin_plot-0.5, x_for_violin_plot+0.5];
else
    scatter_x = x_for_violin_plot;
end
num_rows = length(example_ind);

for i = 1:num_rows
    i_subplot = inset_subplots(i);
%     subplot(num_rows, num_cols, i_subplot);
    subplot_tight(num_rows, num_cols, i_subplot, [0.01]);
    
    i_data = example_ind(i);
    
    X = all_X(i_data,:);
    if use_reconstructed
        X_recon = all_X_reconstructed(i_data,:);
        U = all_paths{i_data}.U;
    end
    
    plot(ts, X, line_opt{:}, 'color', cmap(i,:))
    hold on
    if use_reconstructed
        plot(ts2, X_recon, line_opt{:})
        plot(ts3, 0.5*U*max(X) - max(abs(X)), 'k', line_opt{:})
    end
    xlim([1,200])
    ylim([-max(abs(X)), max(abs(X))])
    yticks([])
    
    this_accuracy = acc(i_data);
    y = this_accuracy;
    
    if i==1
        title(top_title_str)
    end
    if i<num_rows
        xticks([])
    else
        xlabel("Time (\mu s)")
        if to_plot_legend
            leg = legend({"Time series", "Reconstruction", "Control signal"});
        end
    end
    
    subplot(121)
    if violin_not_scatter
        plot(violin_x, [y, y], 'color', cmap(i, :), 'linewidth', 3)
    else
        scatter(scatter_x(i), y, [], cmap(i, :), 'linewidth', 10, ...
            'marker', marker)
    end
end


end

