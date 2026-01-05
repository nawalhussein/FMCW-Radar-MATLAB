clear; clc; close all;

%% Constants
c = 3e8;            % Speed of light
fs = 2e9;           % Sampling frequency  fs >= 2fbmax
fc = 76.5e9;        % Carrier frequency
BW = 1e9;           % Bandwidth
Tc = 3e-6;          % Chirp duration
N = 2048;           % Number of chirps
A_tx = 1;           % Amplitude

% FMCW Parameters
SNR_dB= 10;
s = BW/Tc;
M = round(Tc*fs);             % No. of samples per chirp
R_max = c*fs*Tc/(4*BW);       % Max Range 
v_max = c/(4*fc*Tc);          % Max Velocity 
delta_v = c/(2*fc*(N*Tc));    % Velocity Resolution
delta_R = c/(2*BW);           % Range Resolution

fprintf('Range Res: %.2fm\n',delta_R);
fprintf('Max Range: %.2fm\n',R_max);
fprintf('Velocity Res: %.2fm/s\n', delta_v);
fprintf('Max Velocity: %.2fm/s\n',v_max);

% we should have a time vector to represent the samples on the time axis
% each chirp is of length=M, to represent it, we need to multiply it by Ts so that we have 1Ts, 2Ts,...
% this represents the fast time axis, it is vertical so we take the transpose of the samples matrix
t= (0:M-1)'/fs;

%Transmitted signal
% we use baseband of -0.5 to 0.5 GHz, as if we used 76 GHz to 77 GHz as it is, it would require an impossible fs to achieve Nyquist rate
% fmin=-0.5GHz
Tx= A_tx*exp(2j*pi*(-BW*t/2 + s*t.^2/2));

% noise amp to be added to received signal
signal_power = mean(abs(Tx).^2);
noise_power = signal_power / 10^(SNR_dB/10);
noise_amp = sqrt(noise_power / 2);

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
grid on; 
xlabel('Time'); ylabel('Frequency');

% Targets
targets.R = [50; 122; 200];
targets.v = [+70; -80; +25];

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
%we need to make account for moving targets too
    range = targets.R(i) + targets.v(i)*(n-1)*Tc; 
    RTT = 2*range/c;

% 2nd: get doppler shift
fd= 2*targets.v(i)*fc/c;
dopp_shft= exp(2j*pi*fd*(n-1)*Tc);

% 3rd: we only want the overlap btn the Tx and Rx only
for k=1:length(t)
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
mix(:, n) = Tx.* conj(Rx);  % mixed signal most i 
end

% we need a single chirp of mixer output and Rx for display
mix_chirp= mix(:, 1); 
Rx_chirp= Rx(:, 1);

% Amplitude vs time of received Signal 
figure(3);
plot(t, real(Rx_chirp));
title(sprintf('3. Received Signal (SNR = %d dB)', SNR_dB));
% Amplitude vs time of mixer Output 
figure(4); 
plot(t, real(mix_chirp)); 
title(sprintf('4. Mixed Signal (SNR = %d dB)', SNR_dB));
grid on; 
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

% CFAR Setup 
% The same idea of moving average filter but in a way that is more adaptive
% it makes a new threshold for each point specific only for it 
% RP reference point for averaging , GP guard point ,
% we make something like a band (RP RP .. GP GP .. CUT .. GP GP .. RP RP)
% we ignore GP to avoid leak of target with surroundings (target is not a thin discrete signal it is a thin sine)
% RP the averaged points  
Px_used = 20*log10(abs(mix_fft(1:M_pad/2))); 
Pfa = 1e-8;   % For every 10^6 detected target one of them would be a noise
N_reference_value = 20;  % Total reference cells (left , right)
N_guard_value = 4;     % Total guard cells (left and right)
% Calculate Bias (Alpha)
bias =  Pfa^(-1/N_reference_value) - 1;
% Create Sliding Window Mask to use it to convolve
mask = ones(1, N_reference_value + N_guard_value + 1);
mask(N_reference_value/2 + 1 : N_reference_value/2 + N_guard_value + 1) = 0; % Zero out Guard + CUT
% --- Calculate Threshold ---
% Convolve to get noise floor
noise_estimate = conv(Px_used, mask, 'same');
manual_threshold = 20*log10(noise_estimate * bias);   % dynamic threshold 

% logic to check if point greater than threshold it is a true target, other
% it is not a target it is noise
true_peak = Px_used > manual_threshold;  
% 2. Find peaks ONLY in the true targets that are detected above

[~, loc_index] = findpeaks(abs(true_peak));  % we want to see the loc. index corresponds to how many meters
 peak_values = Px_used(loc_index);

% we lookup for the (loc. index)^th item in the range_axis something like mapping it to meters
peak_ranges = range_axis(loc_index);    

figure(5);

% 1. Plot the Original Signal 
G4 = plot(range_axis, Px_used, 'b', 'LineWidth', 1, 'DisplayName', 'Signal Power'); 
hold on;

% 2. Plot the CFAR Threshold (Red Dashed Line)
G5 = plot(range_axis, manual_threshold, 'r--', 'LineWidth', 2, 'DisplayName', 'CFAR Threshold');

% 3. Plot the Detected Markers
% We plot them at the locations found, using the peak values
G6 = plot(peak_ranges, peak_values, 'gv', 'MarkerFaceColor', 'g', 'DisplayName', 'Detected Targets');

title(sprintf('5. CFAR Detection Results (SNR = %d dB)', SNR_dB));
xlabel('Range (m)'); 
ylabel('Power (Linear)');
grid on;

% --- Printing Results ---
fprintf('---------------------------------------------------\n');
fprintf('Detected Targets Range | Original Targets Range\n');
    for idx = 1:length(loc_index)
        fprintf('%-7.2fm               | %-7.2fm\n', peak_ranges(idx), targets.R(idx));
    end

    legend([G4, G5, G6], 'Raw Signal', 'CFAR Threshold', 'Target Peak');

%% Range-Doppler Map & Velocity

% to get the delta phase shift, we perform 2D fft across the slow & fast time axes
fft_2D= fft2(mix,M_pad,N); %computes fft over the rows and then over the columns

% the problem here is that the zero is at the first by default as a result
% of fft so we use fft shift to put the zero at the middle
% so we need to shift the fft
fft_2D_shifted= abs(fftshift(fft_2D,2));
fft_2D_shifted= fft_2D_shifted(1:M_pad/2, :); % here we put zero at the middle and removed the negative more as it is ignored for range

% we need to generate the velocity axis to plot on it (we already generated the range axis)
% we can use the velocity resolution we calculated at the beginning
vel_axis = linspace(v_max, -v_max, N); %create axis from +VMax to -VMax as the fftshift shifts the polarity

% Range-Doppler Map
figure(6);
imagesc(vel_axis, range_axis, 20*log10(fft_2D_shifted)); % here we used this function to plot a map that uses the full range of colors in the colormap.
% we used db scale instead of linear to give it more power so that the targets appear clearer in the plot
axis xy;                 % Put range 0 at the bottom 
colormap hot;           
xlabel('Velocity (m/s)');
ylabel('Range (m)');
title(sprintf('6. Range-Doppler Map(SNR = %d dB)', SNR_dB));

% Doppler Spectrum
% we need to take a slice of the 2D map of each detected target
figure(7);
xlabel('Velocity (m/s)');
ylabel('Amplitude (m)');
title(sprintf('7. Doppler Spectrum (SNR = %d dB)', SNR_dB));

hold on;
% print detected targets
fprintf('---------------------------------------------------\n');
fprintf('Detected Targets Velocity | Original Targets Velocity\n');

%cfar for velocity :
for i=1:length(loc_index)
dopp_spect = fft_2D_shifted(loc_index(i), :);
Px_vel = (dopp_spect).^2;
Pfa = 1e-8;  
N_reference_vel = 25; 
N_guard_vel = 16 ;
bias_vel = Pfa^(-1/N_reference_vel) - 1; 
% Create a mask for the sliding window
mask_vel = ones(1, N_reference_vel + N_guard_vel + 1);
mask_vel(N_reference_vel/2 + 1 : N_reference_vel/2 + N_guard_vel + 1) = 0; % Clear guard + CUT
% Calculate noise power using sliding window (convolution)
noise_estimate_vel = conv(Px_vel, mask_vel, 'same');
manual_threshold_vel = noise_estimate_vel * bias_vel;
% Target detection
detections_vel = Px_vel > manual_threshold_vel;

% find the peaks
% 1. Apply the mask to the real signal (Zero out noise)
masked_signal_for_peak = Px_vel .* detections_vel;  
[~, vel_indices] = max(masked_signal_for_peak);

peak_ranges = vel_axis(vel_indices);    
% we lookup for the (loc. index)^th item in the vel_axis
     G1 = plot(vel_axis,20*log10(Px_vel));   % signal plot
     G2 = plot(vel_axis, 20*log10(manual_threshold_vel), 'r--', 'LineWidth', 2, 'DisplayName', 'CFAR Threshold');

    if peak_ranges > 0, status = 'Approaching';
    elseif peak_ranges < 0, status = 'Receding';
    else, status = 'Stationary'; end
    fprintf('%+8.2f m/s (%-10s)| %+8.2f m/s\n', peak_ranges, status, targets.v(i));

        real_peak_heights = Px_vel(vel_indices); % to get real height of the pulse 
    G3 = plot(peak_ranges, 20*log10(real_peak_heights), 'bo', 'MarkerFaceColor','b'); % puts circle over the peaks
end % the for loop of getting velocity slice
legend([G1, G2, G3], 'Raw Signal', 'CFAR Threshold', 'Target Peak');
hold off;