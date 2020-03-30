function error_rate = cal_OFDM(interval_ratio, N_s, SNR)
% Definition and Declaration 
% Assume there are N_s subcarrier wave
% N_s different symbols represented by integer (1, 2)
S_k = randi(2, N_s, 1)';

% impulse response of the channel
h_channel = [0.2, 0, 0, 0.7, 0, 0.1];
H_k = fft(h_channel, N_s);
L = 3;

% signal-to-noise ratio(dB)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transmitter
S_n = N_s*ifft(S_k);
S_n = [S_n(N_s-L+2:end), S_n];  % Cyclic prefix

% Transmission through the channel
C_r = zeros(1, interval_ratio*(N_s+L-1));
for i = 1:(N_s+L-1)
    C_r(interval_ratio*i-interval_ratio+1:interval_ratio*i) = S_n(i)*ones(1, interval_ratio); % Padding S_n to the corresponding size 
end

S_r = conv(C_r, h_channel);           % Multipath effect
S_r_withNoise = awgn(S_r, SNR);       % Gaussian Noise

% Receiver
Receiver_Sample = zeros(1, N_s+L-1);
for i = 1:(N_s+L-1) % Sample discarding L-1 elements caused by the convolution
    Receiver_Sample(i) = S_r_withNoise(interval_ratio*i-floor(interval_ratio/2));
end
Demodulation = fft(Receiver_Sample(L:end))./N_s;

S_n_bar = round(abs(Demodulation./H_k));


% Calculate bit error rate
error = sum(S_n_bar ~= S_k);
error_rate = error/N_s;

end

