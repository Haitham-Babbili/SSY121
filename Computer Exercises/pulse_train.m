clc; close all; clear all

% Input parameters
Tsymb=0.02;       % symbol time
fsymb = 1/Tsymb;  % Symbol rate [symb/s]
Tsamp = 0.002;    % Sampling time
fs = 1/Tsamp;     % sampling frequency [Hz]
fsfd = fs/fsymb;  % Number of samples per symbol (choose fs such that fsfd is an integer for simplicity) [samples/symb]
% Transmitted symbols
x = [(1 + 1i), (1 + 1i), (-1 +1i), (-1 - 1i), (1 - 1i)]/sqrt(2);% x: 1 - QPSK/4-QAM

% Pulse train using "for loop" or "Convolution"
span = 6; %how many symbol times to we want of pulse (note that RC and RRC is not zero outside of symbol time!)
[pulse, t] = rcpuls(0.35,1/fsymb,fs,span);
figure; stem(t,pulse); legend('Pulse signal') % plot pulse signal
figure; stem(eps:Tsymb:5*Tsymb,real(x)); % plot first 5 symbols
hold on;
s_loop = zeros(1,(length(x)+2*span)*Tsymb*fs); 
for i = 1:5
    t_current = t+ (i-1)*Tsymb;
    plot(t_current, real(x(i)*pulse)); % Each symbol times pulse
    idx = (i-1)*Tsymb*fs+1:ceil((i-1+2*span)*Tsymb*fs-1); % The current time index. Time interval is one symbol time between each symbol
    s_loop(idx) = s_loop(idx) + x(i)*pulse;
end
p1 = plot([t,t(end):1/fs:t(end)+5*Tsymb], real(s_loop), 'r', 'LineWidth', 3); % Plot s using for loop
legend(p1, 's using for loop');

% Using convolution method
x_upsample = upsample(x, fsfd);      % Space the symbols fsfd apart, to enable pulse shaping using conv.
figure; stem(real(x(1:5))); hold on
figure; stem(real(x_upsample))

s = conv(pulse,x_upsample); % Convolute 'pulse' with 'x_upsample'
figure;
p2 = plot(Tsamp*(-floor(length(pulse)/2):length(s)-ceil(length(pulse)/2)), real(s), 'b', 'LineWidth', 3);
legend(p2,'s using convolution')

