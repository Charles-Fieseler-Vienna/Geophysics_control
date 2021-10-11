function [perc, X_reconstruction] = var_explained_by_dmdc(X, U, A,...
    average_with_correlation_coef)
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
total_norm = norm(X_clean);
% fprintf("The signal was determined to have %.2f percent signal\n",...
%     total_norm / norm(X));
perc2 = 1 - (norm(X_reconstruction - X) / total_norm);

if average_with_correlation_coef
    perc = mean([perc1, perc2]);
else
    perc = perc2;
end

end

