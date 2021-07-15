%% Data
% Build example dynamics matrix
seed = 17;
n = 2;                  % Matrix size; should be even
m = 1000;               % Number of data points
eigenvalue_min = 0.95;  % Minimum eigenvalue; 1.0 = stable
[X_dmd, A] = test_dmd_dat(n, m, 0, eigenvalue_min, seed);

% Build several random controllers
ctr_timing = [100:105, 200:220, 500:505];
U = zeros(1,m-1);
U(ctr_timing) = 1.0;

% Generate controlled data
x0 = X_dmd(:,1);
B = ones(n,1);
X_true = real(calc_reconstruction_dmd(x0, [], A, B, U));
% Add noise
noise = 0.05;
X_raw = X_true + noise*randn(size(X_true));

% Plot
figure;

subplot(2,1,1)
plot(1:m,X_true(1,:), 'LineWidth',2)
hold on
plot(1:m,X_raw(1,:))
title("Data")
legend(["Truth", "Noisy Observations"])

subplot(2,1,2)
plot(1:m-1,U)
title("True Control signal")



%% Embedding
% We only see one dimension
X = X_raw(1,:);

% Augment with (excessive) delays
num_delays = 10;
X_embed = time_delay_embed(X, num_delays);

% Reduce dimensionality with SVD; try to learn optimal dimension
num_dimensions = optimal_truncation(X_embed);
[~, D, V] = svd(X_embed); 
X = D(1:num_dimensions,1:num_dimensions)*V(:,1:num_dimensions)';

% Plot new dimensions
figure;
plot(X')
title(sprintf("%d Dimensions revealed using Time Delay Embedding", num_dimensions))
legend()

figure;
plot(X_true')
title('The original dimensions')

%% Plot frequencies of data

figure;

subplot(2,1,1)
fs = 100;               % sampling frequency
n = length(X_true(1,:));
X_hat = fft(X_true(1,:));
Y = fftshift(X_hat);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y).^2/n;     % zero-centered power
plot(fshift,powershift)
title('Power spectrum of the data')
xlabel('Frequency (check units)')

subplot(2,1,2)
n = length(U(1,:));
U_hat = fft(U(1,:));
Y = fftshift(U_hat);
powershift = abs(Y).^2/n;     % zero-centered power
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
plot(fshift,powershift)
title('Power spectrum of the control signal')
xlabel('Frequency (check units)')

%% Power spectra of just a sine wave + controller
fs = 100;               % sampling frequency
t = 0:(1/fs):(10-1/fs); % time vector
S = cos(2*pi*t(2:end));

figure;

% JUST DATA
subplot(2,3,1)
plot(S)
title('Data')

subplot(2,3,4)
n = length(S);
X_hat = fft(S);
Y = fftshift(X_hat);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y).^2/n;     % zero-centered power
plot(fshift,powershift)
title('Power spectrum of the data')
xlabel('Frequency (check units)')


% JUST CONTROLLER
subplot(2,3,2)
plot(U)
title('Controller')

subplot(2,3,5)
X_hat = fft(U);
Y = fftshift(X_hat);
powershift = abs(Y).^2/n;     % zero-centered power
plot(fshift,powershift)
title('Power spectrum of the controller')
xlabel('Frequency (check units)')

% SUM
subplot(2,3,3)
plot(S + U)
title('Data + controller')

subplot(2,3,6)
X_hat = fft(S + U);
Y = fftshift(X_hat);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y).^2/n;     % zero-centered power
plot(fshift,powershift)
title('Power spectrum of the data')
xlabel('Frequency (check units)')

%% Test: what about notch-filter embedded signal?
% X = S + U; % true cosine + controller
X = X_true(1,:); % true cosine + controller

% First, get the frequency of the signal
n = length(X);
X_hat = fft(X);
Y = fftshift(X_hat);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y).^2/n;     % zero-centered power

% Get frequency to remove
figure
findpeaks(powershift, 'NPeaks', 1, 'MinPeakHeight',10);
[~, peak_loc] = findpeaks(powershift, ...
    'NPeaks', 1, 'MinPeakHeight',max(powershift)/2);
peak_freq = abs(fshift(peak_loc));
title(sprintf('Peak found at %.2f', peak_freq))

band = peak_freq*[0.9,1.1];%peak_freq + [-2, 2];
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',band(1),'HalfPowerFrequency2',band(2), ...
               'DesignMethod','butter','SampleRate',fs);

% Actually filter
X_filt = filtfilt(d, X);

% Plot
figure
subplot(2,1,1)
plot(X')
title('Original data')

subplot(2,1,2)
plot(X_filt')
title('Filtered data')
%% Test: TONS of embeddings

% Augment with (excessive) delays
num_delays_vec = [2, 4, 10];
for i = 1:length(num_delays_vec)
    num_delays = num_delays_vec(i);
    
    X_embed = time_delay_embed(X_raw(1,:), num_delays);

    % Test: What is the embedded residual?
    X1 = X_embed(:, 1:end-1);
    X2 = X_embed(:, 2:end);
    A2 = X2/X1;

    res2 = X2 - A2*X1;

    figure;
    subplot(2,1,1)
    plot(res2')
    title(sprintf('%d-delay embedded residual', num_delays))
    subplot(2,1,2)
    plot(U(num_delays:end))
    title('True controller (offset to align)')
end

%% Plot integral of the embedded signal as well

% Augment with (excessive) delays
num_delays_vec = [2, 4, 10];
for i = 1:length(num_delays_vec)
    num_delays = num_delays_vec(i);
    
    X_embed = time_delay_embed(X_raw(1,:), num_delays);

    % Test: What is the embedded residual?
    X1 = X_embed(:, 1:end-1);
    X2 = X_embed(:, 2:end);
    A2 = X2/X1;

    res2 = X2 - A2*X1;

    figure;
    subplot(3,1,1)
    plot(res2')
    title(sprintf('%d-delay embedded residual', num_delays))
    
    subplot(3,1,2)
    plot(U(num_delays:end))
    title('True controller (offset to align)')
    
    subplot(3,1,3)
    plot(cumsum(res2'))
    title('Integrated learned controller')
end
