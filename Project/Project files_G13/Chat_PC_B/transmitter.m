% COMPLETE THE TRANSMITTER!

% pack = message to be transmitted (consists of 432 bits from the GUI, always!)
% fc = carrier frequency
function transmitter(pack,fc)
%% Transmitter

% Input parameters
data = pack;                                                                % The transmitted data
fs = 22050;                                                                 % Sampling frequency [Hz]
Tsamp = 1/fs;                                                               % Sampling time
const = [(1 + 1i), (1 - 1i), (-1 -1i), (-1 + 1i)]/sqrt(2);                  % Constellation 1 - QPSK/4-QAM
M = length(const);                                                          % Number of symbols in the constellation
bitperSymb = log2(M);                                                       % Number of bits per symbol
alpha = 0.4;                                                                % Roll off factor
fsfd = 120;                                                                 % Number of samples per symbol [samples/symb] (it better to have it as interger)
Rs = fs/fsfd;                                                               % Calculate Symbole rate [symbole/Hz]
Ts = 1/Rs;                                                                  % Symbole Time/Duration [Symbole/sec] 
Rb = Rs*bitperSymb;                                                         % bit rate [bit/sec]
Bw = (1+alpha)/(2*Ts);                                                      % Bandwidth
msg = buffer(data, bitperSymb)';                                            % Bits to Messages, Group bits into bits per symbol
msg_idx = bi2de(msg, 'left-msb')'+1;                                        % Message to symbol index
x = const(msg_idx);                                                         % Look up symbols on constellation
preamb = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];     % 26 bits from Barker code
x = [preamb, x];                                                            % tx_symbol, TX vector

% Upsampling, more samble,i.e, space, for the symbols by fsfd apart
x_upsample = upsample(x, fsfd);                                             % to enable pulse shaping using conv.

% Convert the symbole to base-band signal (train of pulses)
span = 6;                                                                   % Pulse width
pulse = rtrcpuls(alpha,1/Rs,fs,span);                                       % Create root raised cosine pulse shape
s_train = conv(pulse,x_upsample);                                           % Convolve symbole with pulse ==> Baseband signal.

% Passband Signal 
s_train_TX = sqrt(2)*s_train.*exp(-1i*2*pi*fc*(0:length(s_train)-1)*Tsamp); % Modulation, Put the baseband signal on carrier signal 
s_train_TX_real = real(s_train_TX);                                         % Only Real part of the signal will be send
s_train_TX_real = s_train_TX_real./max(abs(s_train_TX_real));               % Normalize signal to 1 to prevent clipping of waveforms

tx_signal = s_train_TX_real;
player = audioplayer(tx_signal, fs);       %create an audioplayer object to play the noise at a given sampling frequency
playblocking(player); % Play the noise 
