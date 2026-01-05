clear; clc; close all;

%% Constants
c = 3e8;            % Speed of light
fs = 2e9;           % Sampling frequency
fc = 76.5e9;        % Carrier frequency
BW = 1e9;           % Bandwidth
Tc = 3e-6;          % Chirp duration
N = 4096;           % Number of chirps
A_tx = 1;           % Amplitude

% FMCW Parameters
SNR_dB= 10;
s = BW/Tc;
M = round(Tc*fs);             % No. of samples per chirp
R_max = c*fs/(4*s);           % Max Range
v_max = c/(4*fc*Tc);          % Max Velocity
delta_v = c/(2*fc*N*Tc);      % Velocity Resolution
delta_R = c/(2*BW);           % Range Resolution

fprintf('Range Res: %.2fm\n',delta_R);
fprintf('Max Range: %.2fm\n',R_max);
fprintf('Velocity Res: %.2fm/s\n', delta_v);
fprintf('Max Velocity: %.2fm/s\n',v_max);

% we should have a time vector to represent the samples on the time axis
% each chirp is of length=M, to represent it, we need to multiply it by Ts so that we have 1Ts, 2Ts,...
% this represents the fast time axis, it is vertical so we take the transpose of the samples matrix
t= (0:M-1)'/fs;

% Transmitted signal
% we use baseband of -0.5 to 0.5 GHz, as if we used 76 GHz to 77 GHz as it is, it would require an impossible fs to achieve Nyquist rate
% fmin=-0.5GHz
Tx= A_tx*exp(2j*pi*(-BW*t/2 + s*t.^2/2));

% noise amp to be added to received signal
signal_power = mean(abs(Tx).^2);
noise_power = signal_power / 10^(SNR_dB/10);
noise_amp = sqrt(noise_power / 2); %In a balanced system, half of the power goes into the Real component and the other half goes into the Imaginary component

% Amplitude vs time of transmitted signal
figure(1);
plot(t,real(Tx));
title(sprintf('1. Transmitted Signal (SNR = %d dB)', SNR_dB));
xlabel('Time'); ylabel('Amplitude');

% Frequency vs time of transmitted signal
freq=diff(unwrap(angle(Tx))) * fs/(2*pi); % change in phase=frequency
% we use unwrap as angle ranges from -pi to pi only but we want to see the whole change
% we mult. it by fs to get change per second instead of change per sample
% diff makes the size of freq matrix less than the size of t matrix by 1
% we assume that the last element of freq matrix is like its end
freq=[freq;freq(end)];
freq_rep= repmat(freq,5); % display 5 sawtooth chirps
t_freq=(0:length(freq_rep)-1)/ fs; % x-axis must be the same size as y-axis

figure(2);
plot(t_freq,freq_rep);
title(sprintf('1. Transmitted Signal (Frequency)(SNR = %d dB)', SNR_dB));
xlabel('Time'); ylabel('Frequency');

% Targets
targets.R = [50; 180];
targets.v = [+100; -80];

% Received signal generation
% The received signal is only the trans. signal but shifted and delayed
% we need to find Rx for each chirp

mix= zeros(M,N); % allocate memory for mixing output
% this represesnts the slow time (rows) vs fast time (columns)
% each chirp has M samples, so column must have M rows to fit all the data 
% we repeat this for N chirps, so we need N columns

for n= 1:N % we get Rx for each chirp
    % allocate memory for rx
    Rx= zeros(M,1); %reset rx for every new chirp
    % the sz of Rx is exactly the same as Tx (M)
for i= 1:length(targets.R) % loop on all targets 

% 1st: get RTT (delay)
% we need to make account for moving targets too
    range = targets.R(i) + targets.v(i)*(n-1)*Tc; %used to calculate the distance of the target for the nth chirp in a sequence
    RTT = 2*range/c;

% 2nd: get doppler shift
fd= 2*targets.v(i)*fc/c;
dopp_shft= exp(2j*pi*fd*(n-1)*Tc);

% 3rd: we only want the overlap btn the Tx and Rx only
for k=1:M
if t(k)>=RTT
Rx(k)= Rx(k) + A_tx*exp(2j*pi*(-BW*(t(k)-RTT)/2 + s*(t(k)-RTT).^2/2))*dopp_shft;
else
        Rx(k) = 0; 
    end
end
end
% 4th: add noise (SNR=10db)
Rx = Rx + noise_amp * (randn(M,1) + 1j*randn(M,1));

% 5th: mixing
% Take the results of the nth chirp and put them into the nth column
mix(:, n) = Tx.* conj(Rx); 
end

% we need a single chirp of mixer output and Rx for display
mix_chirp= mix(:, 1); 
Rx_chirp= Rx(:, 1);

% Amplitude vs time of received Signal 
figure(3);
plot(t, real(Rx_chirp));
title(sprintf('3. Received Signal (SNR = %d dB)', SNR_dB));
xlabel('Time'); ylabel('Amplitude');

% Amplitude vs time of mixer Output 
figure(4); 
plot(t, real(mix_chirp)); 
title(sprintf('4. Mixed Signal (SNR = %d dB)', SNR_dB));
xlabel('Time'); ylabel('Amplitude');

% Range-FFT Plot
% to decrease the error, we apply padding which is adding zeros to the end of the fft matrix, to get more resolution points
M_pad = M * 4; % as if we zoomed the fft to 4x to have better resolution and smaller bins
mix_fft= fft(mix(:,1), M_pad); % apply fft on the 1st chirp with adding zeros
k= 0:(M_pad/2)-1; % length of k is the length of time signal (no.of samples)
% we divide by 2 to take the +ve part only
% we need to define the frequency bin/ resolution
freq_bin= fs/M_pad;
freq_axis= k*freq_bin; % frequency axis: frequency that corresponds to each k
range_axis= (freq_axis * c)/(2*s); % convert frequency axis to range axis

% find the peaks
[peak, loc_index] = findpeaks(abs(mix_fft(1:M_pad/2)),'MinPeakHeight', 0.3 * max(abs(mix_fft(1:M_pad/2))));
peak_ranges = range_axis(loc_index); % we want to see the loc. index corresponds to how many meters
% we lookup for the (loc. index)^th item in the range_axis

figure(5);
plot(range_axis, abs(mix_fft(1:M_pad/2))); 
 % we take the abs bc we care about the mag only above, we took one chirp (column) and plotted the its full range profile
hold on;
plot(peak_ranges, peak, 'rv', 'MarkerFaceColor','r'); % puts a red reversed triangle over the peaks
title(sprintf('5. Range FFT (SNR = %d dB)', SNR_dB));
xlabel('Range (m)'); ylabel('Amplitude');

% print detected targets
fprintf('---------------------------------------------------\n');
fprintf('Detected Targets Range | Original Targets Range\n');
for idx = 1:length(loc_index)
fprintf('%-7.2fm               | %-7.2fm\n', peak_ranges(idx),targets.R(idx));
end

%% Range-Doppler Map & Velocity

% to get the delta phase shift, we perform 2D fft across the slow & fast time axes
fft_2D= fft2(mix,M_pad,N); %computes fft over the rows and then over the columns

% the problem here is that fft2 sets the zero at the first bin so we use
% fftshift to set the zero at the middle again
fft_2D_shifted= abs(fftshift(fft_2D,2));
fft_2D_shifted= fft_2D_shifted(1:M_pad/2, :);

% we need to generate the velocity axis to plot on it (we already generated the range axis)
% we can use the velocity resolution we calculated at the beginning
vel_axis = linspace(v_max, -v_max, N); %create axis from +VMax to -VMax as the fftshift shifts the polarity


% Doppler Spectrum
% we need to take a slice of the 2D map of each detected target
figure(6);
xlabel('Velocity (m/s)');
ylabel('Amplitude (m)');
title(sprintf('6. Doppler Spectrum (SNR = %d dB)', SNR_dB));
hold on;
% print detected targets
fprintf('---------------------------------------------------\n');
fprintf('Detected Targets Velocity | Original Targets Velocity\n');
est_velocities = zeros(length(loc_index), 1);
for i=1:length(loc_index)
    dopp_spect = fft_2D_shifted(loc_index(i), :); %Go to the row where my target is located (Range) and give me the entire horizontal line of data (Velocity) across all chirps
    plot(vel_axis, dopp_spect);
    % Find peaks
    [v_peak, v_index] = max(dopp_spect);
    est_v= vel_axis(v_index);
    est_velocities(i) = est_v;
    if est_v > 0, status = 'Approaching';
    elseif est_v < 0, status = 'Receding';
    else, status = 'Stationary'; end
    fprintf('%+8.2f m/s (%-10s)| %+8.2f m/s\n', est_v, status, targets.v(i));
     plot(est_v, v_peak, 'rv', 'MarkerFaceColor','r'); % puts a red reversed triangle over the peaks
end

% Range-Doppler Map
figure(7);
imagesc(vel_axis, range_axis, 20*log10(fft_2D_shifted));
% we used db scale instead of linear to give it more power so that the targets appear clearer in the plot
axis xy;     % Put range 0 at the bottom 
colormap hot;           
xlabel('Velocity (m/s)');
ylabel('Range (m)');
title(sprintf('7. Range-Doppler Map(SNR = %d dB)', SNR_dB));
hold on; 
plot(est_velocities, peak_ranges, 'wx', 'LineWidth', 2, 'MarkerSize', 10);