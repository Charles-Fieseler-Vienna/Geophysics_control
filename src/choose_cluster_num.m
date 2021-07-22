function [best_gmm,all_AIC,all_gmm] = choose_cluster_num(X, k_max, opt)
% Chooses the "correct" number of gmm clusters
%   Maximizes AIC
%   Starts at k=2

all_gmm = cell(k_max-1,1);
all_AIC = zeros(k_max-1,1);

for i = 2:k_max
    all_gmm{i-1} = fitgmdist(X, i, opt{:});
    all_AIC(i-1) = all_gmm{i-1}.BIC;
end

[val, i] = min(all_AIC);
fprintf('Best cluster determined to be %d\n', i+1)
best_gmm = all_gmm{i};

end

