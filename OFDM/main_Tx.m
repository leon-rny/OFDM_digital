clear; close all;
mod_modus = 3;  % 0: BPSK, 1: QPSK, 2: 16-QAM, 3: 64-QAM
N = 12;          % Anzahl der OFDM-Symbole also mit Servie+DAta

%% Preamble (SYNC) Generator
% Short Preamble
S = sqrt(13/6) * [0, 0, 1+1j, 0, 0, 0, -1-1j, 0, 0, 0, 1+1j, 0, 0, 0, -1-1j, 0, 0, 0, -1-1j, 0, 0, 0, 1+1j, 0, 0, 0, 0, 0, 0, 0, -1-1j, 0, 0, 0, -1-1j, 0, 0, 0, 1+1j, 0, 0, 0, 1+1j, 0, 0, 0, 1+1j, 0, 0, 0, 1+1j, 0,0];
S_filled = [0 S(28:end) zeros(1,11) S(1:26)];
S_ifft = ifft(S_filled, 64);
r_short = [S_ifft S_ifft S_ifft(1:32)];

% Long Preamble 
L = [1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 0, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1];
L_filled = [zeros(1,6) L(1:end) zeros(1,5)];
L_shifted = fftshift(L_filled);
L_ifft = ifft(L_shifted, 64);
r_long = [L_ifft(33:64) L_ifft L_ifft];

%% Datascrambler
seed = [1 0 1 0 1 0 1];
output = scrambler(seed, mod_modus, N);
scrambledData = output.scrambledData;
numBits = output.numBits;

%% Signal Field Generator           
rate = [1 1 0 1];                                       % Datenrate: 6 Mbits/s aus Tabelle 80 S. 14
length = [0 0 0 0 0 0 1 0 0 0 0 0];                     % Datenlänge: 128 bits
tail = zeros(1,6);                                      % Auffüllen mit 0 als Platzhalter
signal = [rate 0 length 0 tail zeros(1,24)];            % Zusammenführen zu einem Frame mit 0 als Guards
signal_mapped = mapping(0, signal);                     % Signal: Mapping 
signal_pilot = pilot(signal_mapped);                    % Signal: Pilot Insertion
%signal_reshaped = [signal_pilot(33:64) signal_pilot(1:32)];  % Signal: Reshaped für IFFT
signal_shifted = fftshift(signal_pilot);                
signal_ifft = ifft(signal_shifted);                     % Signal: IFFT
GI_signal = signal_ifft(end-15:end);                    % Signal: GI

Init = zeros(1, 7);                                     % Sync für Descrambler
reserved = zeros(1, 9);                                 % Platzhalter for future use 
service = [Init reserved]; 
service_data = [service scrambledData(1:numBits/N-16)]; % Service Daten
service_data_mapped = mapping(mod_modus, service_data); % Servicedaten: Mappen
service_data_pilot = pilot(service_data_mapped);        % Servicedaten: Pilot Insertion
%service_reshaped = [service_data_pilot(33:64) service_data_pilot(1:32)];  % Servicedaten: Reshaped für IFFT
service_data_shifted = fftshift(service_data_pilot);    
service_data_ifft = ifft(service_data_shifted);         % Servicedaten: IFFT
GI_service = service_data_ifft(end-15:end);             % Servicedaten: GI

%% Mapper
scrambled_remaining = scrambledData(numBits/N-15:end);
data_ifft = [];
quatsch = [];
for i = 1:N-1
    tmp = mapping(mod_modus, scrambled_remaining((numBits/N*(i-1))+1:numBits/N*i)); % Mapping
    tmp_pilot = pilot(tmp);                                                         % Pilot Insertion
    quatsch = [quatsch tmp_pilot];

    %tmp_pilot_reshaped = [tmp_pilot(33:64) tmp_pilot(1:32)];
    tmp_shifted = fftshift(tmp_pilot);
    tmp_ifft = ifft(tmp_shifted, 64);                                               % IFFT
    GI = tmp_ifft(end-15:end);
    tmp_with_GI = [GI tmp_ifft];                                                    % GI
    data_ifft = [data_ifft tmp_with_GI];
end

quatsch = reshape(quatsch, [64, N-1]); 

%% Zusammenfügen
r_data = [r_short r_long GI_signal signal_ifft GI_service service_data_ifft data_ifft];
%y_descrambled = descrambler(scrambledData, mod_modus, N);

df = exp(-1i * 2 * pi * 300e3 * (1:size(r_data,2)) / 20e6);

r_data = r_data .* df;


%% Plot
%% Plot
%% Plot

%Zeitginal long Präampel
f1=figure(1);

plot(real(r_short));
hold on;
plot(imag(r_short));
ylabel('Amplitude');
xlabel('Sample');
legend('Realteil', 'Imaginärteil');
grid on;

%Zeitginal long Präampel
f2=figure(2);

plot(real(r_long));
hold on;
plot(imag(r_long));
xline(80, '--k');
ylabel('Amplitude');
xlabel('Sample');
legend('Realteil', 'Imaginärteil');
grid on;

% Zeitsignal komplett
f3=figure(3);

plot(real(r_data));
hold on;
plot(imag(r_data));
ylim([-0.3 0.3]);
ylabel('Amplitude');
xlabel('Sample');
legend('Realteil', 'Imaginärteil');
xline(160, '--k');
xline(320, '--k');
length_line = size(r_data,2);
for i=320:80:length_line-80
   xline(i, '--k');
end
grid on;
legend('Realteil', 'Imaginärteil');

% exportgraphics(f1, 'shortpreambel_tx.pdf', 'ContentType','vector');
% exportgraphics(f2, 'longpreambel_tx.pdf', 'ContentType','vector');
% exportgraphics(f3, 'komplett_tx.pdf', 'ContentType','vector');

%% Sendevorgang
% p_Tx = 60;
% f_Tx = 2500;
% 
% OpenFPGA;
% Open_RFFE_UART;
% 
% RFFE_SetLO_Frequency(f_Tx);
% SetTxPower(p_Tx);
% r_data = 512 * r_data;
% 
% UploadAWG_Signal(Data_TX);
% StartAWG_single;
