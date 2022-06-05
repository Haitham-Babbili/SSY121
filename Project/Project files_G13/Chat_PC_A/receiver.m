% RECEIVER 
%
% This is the receiver structure that you will have to complete.
% The function: receiver(fc) is a setup function for the receiver. Here,
% the audiorecorder object is initialized (see help audiorecorder or
% MATLAB's homepage for more information about the object).
% 
% The callback function audioTimerFcn() is a callback function that is
% triggered on a specified time interval (here it is determined in the
% setup function, by the variable time_value)
% 
% Your task is to extend this code to work in the project!
%%

function [audio_recorder] = receiver(fc)
fs = 22050;                                                                 %sampling frequency
audio_recorder = audiorecorder(fs,24,1);                                    %create the recorder

audio_recorder.UserData.fs = fs;                                            %allocate for constellation
audio_recorder.UserData.fc  = fc;                                           %allocate for eye diagram
audio_recorder.UserData.lastStop = 0;                                       %added to remove extra signal after previous success

%attach callback function
time_value = 0.5;                                                           % how often the function should be called in seconds
set(audio_recorder,'TimerPeriod',time_value,'TimerFcn',@audioTimerFcn);     % attach a function that should be called every second, the function that is called is specified below.

%ADD USER DATA FOR CALLBACK FUNCTION (DO NOT CHANGE THE NAMES OF THESE VARIABLES!)
audio_recorder.UserData.receive_complete = 0; % this is a flag that the while loop in the GUI will check
audio_recorder.UserData.pack  = []; %allocate for data package
audio_recorder.UserData.pwr_spect = []; %allocate for PSD
audio_recorder.UserData.const = []; %allocate for constellation
audio_recorder.UserData.eyed  = []; %allocate for eye diagram


record(audio_recorder); %start recording
end


% CALLBACK FUNCTION
% This function will be called every [time_value] seconds, where time_value
% is specified above. Note that, as stated in the project MEMO, all the
% fields: pwr_spect, eyed, const and pack need to be assigned if you want
% to get outputs in the GUI.
function audioTimerFcn(recObj, event, handles)
%% Receiver

disp('Callback triggered')
disp(' ')

fc = recObj.UserData.fc;                                                    % Carrier frequency
fs = recObj.UserData.fs;                                                    % Sampling frequency [Hz]
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
preamb = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];     % 26 bits from Barker code
span = 6;                                                                   % Pulse width
pulse = rtrcpuls(alpha,1/Rs,fs,span);                                       % RRC Pulse

disp('Parameters are set')
disp(' ')

y = getaudiodata(recObj)';                                                  % Get the data from the recorder
y = y(recObj.UserData.lastStop+1:end);                                      % Remove the previous data chunk
y = y./max(abs(y));                                                         % Normalize the RX data

disp('Signal is captured')
disp(' ')

rx = sqrt(2)*y.*exp(1j*2*pi*fc*(0:length(y)-1)*Tsamp);                      % Demodulation, Put the baseband signal on carrier signal 

disp('Signal is demodulated')
disp(' ')

%% Time sync

preamb_upsample = upsample(preamb, fsfd);                                   % upsamling the known preamble
preamb_train = conv(pulse, preamb_upsample);                                % shape the known preamble like the rx preamble
corr = conv(rx, fliplr(preamb_train));                                      % correlate the processed preamble and demodulated signal to find the similarity

corrThreshold = 10;                                                         % noise threshold
[maxCorr, Tmax] = max(corr);                                                % find location of max correlation
Tx_hat = Tmax - length(preamb_train);                                       % find delay (it will take length preamble)

% pass the signal only if the correlation is not with noise and the delay
% is positive and the whole signal is captured
if((corrThreshold > abs(maxCorr)) || (Tmax < length(preamb_train)) || ((length(y) - Tmax) < fsfd * 216))
    disp('No signal found')
    disp(' ')
    return;
end

disp(['Delay is ', num2str(Tx_hat), ' samples'])
disp(' ')

recObj.UserData.lastStop = length(y);                                       % keep track for the last stop to remove the last proessed data

%% Match filter

MF = fliplr(conj(pulse));                                                   % Create matched filter impulse response

MF_output = filter(MF, 1, rx);                                              % Run received signal through matched filter
MF_output = MF_output(length(MF):end);                                      % remove transient 

disp('The signal is filtered')
disp(' ')

%% Phase sync

rx_preamb = MF_output(Tx_hat+1:fsfd:Tx_hat+length(preamb_upsample));        % extract the recieved preamble symbols
rx_preamb = rx_preamb./mean(abs(rx_preamb));                                % normalize the preamle symbols

calcPhase = (angle(rx_preamb) - angle(preamb));                             % calculated the difference in angle between the known preamble and the received preamble
calcPhase = mod(calcPhase, 2*pi); 
if(var(calcPhase)>5)
    error = 5;
    calcPhase(calcPhase < (2*pi+error) & calcPhase > (2*pi-error)) = calcPhase(calcPhase < (2*pi+error) & calcPhase > (2*pi-error)) - 2*pi;
end                                          % express the angle using the same reference (all angles starts from zero with direction CCW)
averageCalcPhase = mean(calcPhase);                                         % average the phase shift angles to get one value

MF_output = MF_output(Tx_hat+length(preamb_upsample)+1 : Tx_hat+length(preamb_upsample)+216*fsfd); % get exactly the recieved signal
MF_output = MF_output * exp(-1j * averageCalcPhase);                        % correct the phase

disp(['The signal is Phase Shifted with ', num2str(averageCalcPhase * 180/pi), ' degrees'])
disp(' ')

%% Signal to symbols

rx_vec = MF_output(1:fsfd:end);                                             % sample the matched filter signal at symbole frequency
rx_vec = rx_vec./mean(abs(rx_vec));                                         % normalize the received vector

disp('Symbols acquired')
disp(' ')

%% Symbol to bits

% Minimum Eucledian distance detector
% Relate the detection to Detection region
metric = abs(repmat(rx_vec.',1,4) - repmat(const, length(rx_vec), 1)).^2;   % Compute the distance to each possible symbol
[tmp, msg_hat] = min(metric, [], 2);                                        % find the closest for each received symbol
msg_hat = msg_hat'-1;                                                       % Get the index of the symbol in the constellation
msg_hat = de2bi(msg_hat, 2, 'left-msb')';                                   % Convert symbols to bits
data_hat = msg_hat(:)';                                                     % Write as a vector

disp('Data acquired')
disp(' ')

bits = data_hat(1:432);
x = rx_vec(1:216);
pulse_train = MF_output;

disp('Finished processing');
disp(' ')
%------------------------------------------------------------------------------
% HOW TO SAVE DATA FOR THE GUI
%   NOTE THAT THE EXAMPLE HERE IS ONLY USED TO SHOW HOW TO OUTPUT DATA
%------------------------------------------------------------------------------

% Step 1: save the estimated bits
recObj.UserData.pack = bits;

% Step 2: save the sampled symbols
recObj.UserData.const = x;

% Step 3: provide the matched filter output for the eye diagram
recObj.UserData.eyed.r = pulse_train;
recObj.UserData.eyed.fsfd = fsfd;

% Step 4: Compute the PSD and save it. 
% !!!! NOTE !!!! the PSD should be computed on the BASE BAND signal BEFORE matched filtering
[pxx, f] = pwelch(pulse_train,1024,768,1024, fs); % note that pwr_spect.f will be normalized frequencies
f = fftshift(f); %shift to be centered around fs
f(1:length(f)/2) = f(1:length(f)/2) - fs; % center to be around zero
p = fftshift(10*log10(pxx/max(pxx))); % shift, normalize and convert PSD to dB
recObj.UserData.pwr_spect.f = f;
recObj.UserData.pwr_spect.p = p;

% In order to make the GUI look at the data, we need to set the
% receive_complete flag equal to 1:
recObj.UserData.receive_complete = 1; 
    
end
