function [perc, X_reconstruction] = var_explained_by_dmdc(X, U, A, B)
% Percentage of variance of a signal explained by a DMDc model
%   If A and B are not passed, then they are calculated using '/'

% Loosely get the "actual" signal
% f = @(X) movmean(movmean(X,11),3); % TODO
f = @(X) X;
X_clean = f(X);


if ~exist('A', 'var')
    [~, ~, X_reconstruction] = simple_dmdc(X_clean, U);
end

% Test: correlation
perc = corrcoef(X_reconstruction(1,:), X_clean(1,:));
% perc = cov(X_reconstruction(1,:), X_clean(1,:));
perc1 = perc(1,2);

% Test: just return L2 error
total_norm = norm(X_clean);
% fprintf("The signal was determined to have %.2f percent signal\n",...
%     total_norm / norm(X));
perc2 = 1 - (norm(X_reconstruction - X) / total_norm);

perc = mean([perc1, perc2]);
end

