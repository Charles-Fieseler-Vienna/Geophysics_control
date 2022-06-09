function [perc, X_reconstruction] = var_explained_by_dmdc(X, U, A,...
    norm_of_only_first_dim, average_with_correlation_coef)
% Currently: returns the average of the variance explained and the
% correlation coefficient
%
% Percentage of variance of a signal explained by a DMDc model
%   If A and B are not passed, then they are calculated using '/'

% Loosely get the "actual" signal
% f = @(X) movmean(movmean(X,11),3); % TODO
f = @(X) X;
X_clean = f(X);


if ~exist('A', 'var') || isempty(A)
    [~, ~, X_reconstruction] = simple_dmdc(X_clean, U);
end
if ~exist('average_with_correlation_coef', 'var')
    average_with_correlation_coef = true;
end

if average_with_correlation_coef
    % Test: correlation
    perc = corrcoef(X_reconstruction(1,:), X_clean(1,:));
    % perc = cov(X_reconstruction(1,:), X_clean(1,:));
    perc1 = perc(1,2);
end

% Test: just return L2 error
if norm_of_only_first_dim
    total_norm = norm(X_clean(1,:));
    perc2 = 1 - (norm(X_reconstruction(1,:) - X(1,:)) / total_norm);
else
    total_norm = norm(X_clean);
    perc2 = 1 - (norm(X_reconstruction - X) / total_norm);
end

if average_with_correlation_coef
    perc = mean([perc1, perc2]);
else
    perc = perc2;
end

end

