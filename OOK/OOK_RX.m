clear; close all;
load('rx_difference.mat');
%load('rx_signal_one_burst_right.mat');
%load('OOK_Tx_Signal.mat');
dataIN10Folge = repmat([1, 0], 1, 64);
%load('rx_differnce_artificial.mat');
% y=xdbb_if.';
%% Parameter RX
lna_mode = 1;                                   % ca. 15dB
g = 20;                                          % Gain von RX
freq = 2500;                                    % Einstellung der Trägerfrequenz                        
Threshold = 121;
RecordLength = 3820;
RecordLengthReal = 16384;

f_s=20e6;                                       % Schwelle der Empfangssignals
M = 20;
f_IF = 500e3;                                   % Zwischenfrequenz von 500 kHz
n = 0:RecordLength-1;                               % Index-Vektor 

% %% Empfangen Settings
% 
% OpenFPGA
% Open_RFFE_UART
% 
% %% Funktionsaufrufe
% while 1
%     if RFFE_ReadLO_LockDetect()                 % ganze Zeit detektieren
%         break
%     end
% end
% 
% SetRxGain(lna_mode, g);                                  % Gain setzen
% RFFE_SetLO_Frequency(freq);                              % Empfangsfrequenz setzen
% t_achse = linspace(0, RecordLengthReal./f_s, RecordLengthReal);  % Zeitvektor
% 
% y = RecordSignal(RecordLengthReal,Threshold);                % Signal empfangen

%% abschnibeln
y = y(1:3820);

%% Heruntermischen -> Betragsbildung -> DC Offset Entfernung

exp_signal = (1/500)*exp(-1j*2*pi*500e3*n/f_s);       % Erzeugen der komplexen Exponentialfunktion für Runtersetzen

y1 = y.' .* exp_signal;                             % heruntermischen 
y_abs = abs(y1);                                    % absolutwert nehmen
y_DC_cleaned = (y_abs - mean(y_abs));                 % DC Offset entfernen

%% Erstellen der Preambelsequenz
N = 63;                             % Länge der m-Sequenz
registerLength = 6;                 % Länge des Registers
seed = [1 1 0 1 1 0];               % Seed

register = seed;
taps = [1, registerLength];         % Die Taps sind an den Positionen 1 (x^6) und 6 (x)
p = zeros(N, 1);                    % Initialisieren der Sequenz

for i = 1:N

    p(i) = register(end);                           % Die Sequenz wird durch das letzte Element des Registers generiert
    feedback = mod(sum(register(taps)), 2);         % Berechnen des nÃ¤chsten Bit (Feedback)
    register = [feedback, register(1:end-1)];       % Aktualisieren des Registers
end

trans_p = p.';
pmod = 2*p - 1;                                     % Umwandlung der Bitsequenz: 0 wird zu -1, 1 wird zu +1

deltaImpuls = [1; zeros(M-1, 1)];                       % Erstellen eines Delta-Impulses | [1 0 0 0 0 ...]
hr = ones(1, M);                                        % Impulsantwort des FIR-Filters [1 1 1 1 1 ...]

pmod_uberabgetastet = kron(pmod, deltaImpuls);          % Anwenden der Überabtastung auf die preambelsequenz
pmod_gefiltert = filter(hr, 1, pmod_uberabgetastet);    % Anwenden des FIR-Filters auf die überabgetastete Sequenz

%% Korrelieren, um Daten Beginn zu finden
pmod_gefiltert_trans = pmod_gefiltert.';                    % transponiert

y_conv = conv(y_DC_cleaned, flip(pmod_gefiltert_trans));    % Faltung von empfangensignal  (Betrag & DC entfernt) & Preambel
[y_conv_max, y_conv_max_index] = max(y_conv);               % Maximalen Wert und dessen Index finden, also ist da vor unsere Daten beginnen

%% ohne SAF
% %% bestimmen der richtigen Abtastzeitpunkte
t_neu = (0:RecordLength-1) / f_s;
abtastwerte = [];
y_sample_new = [];

lower = y_conv_max_index + 10;              % erster Abtastwert
upper = lower + 20 * 127;                   % letzert Abtastwert

for j = lower:20:upper                      % richtige Abtastwerte  begin+[10 30 50 70...]
    abtastwerte(end + 1) = j;
end

for i = 1:length(abtastwerte)               % die gesampleten Werte [-0.9  0.85 ...]
    index = abtastwerte(i);
    y_sample_new(i) = y_DC_cleaned(index);
end

%% Entscheider
y_sample_entscheider(y_sample_new>0)=1;
y_sample_entscheider(y_sample_new<0)=0;

%% Descrambler
lfsrState1 = [1 0 1 0 1 0 1];
numBits = 128;
descrambledData = zeros(1, numBits);

for i = 1:numBits
    feedbackBit = xor(lfsrState1(7), lfsrState1(4));
    descrambledData(i) = xor(y_sample_entscheider(i), feedbackBit);
    lfsrState1 = [feedbackBit lfsrState1(1:6)];
end

%% SAF
hr = ones(1, 20);                           
y_SAF = filter(hr, 1, y_DC_cleaned);

%% bestimmen der richtigen Abtastzeitpunkte
t_neu = (0:RecordLength-1) / f_s;
abtastwerte_SAF = [];
y_sample_new_SAF = [];

lower = y_conv_max_index+20;              % erster Abtastwert
upper = lower + 20 * 127;                   % letzert Abtastwert

for j = lower:20:upper                      % richtige Abtastwerte  begin+[10 30 50 70...]
    abtastwerte_SAF(end + 1) = j;
end

for i = 1:length(abtastwerte_SAF)               % die gesampleten Werte [-0.9  0.85 ...]
    index1 = abtastwerte_SAF(i);
    y_sample_new_SAF(i) = y_SAF(index1);
end

%% Entscheider
y_sample_entscheider_SAF(y_sample_new_SAF>0)=1;
y_sample_entscheider_SAF(y_sample_new_SAF<0)=0;

%% Descrambler
lfsrState1 = [1 0 1 0 1 0 1];
numBits = 128;
descrambledData_SAF = zeros(1, numBits);

for i = 1:numBits
    feedbackBit = xor(lfsrState1(7), lfsrState1(4));
    descrambledData_SAF(i) = xor(y_sample_entscheider_SAF(i), feedbackBit);
    lfsrState1 = [feedbackBit lfsrState1(1:6)];
end

%% BER
fehler = xor(dataIN10Folge, descrambledData);   % ohne SAF
anzahl_fehler = sum(fehler);
BER = anzahl_fehler / 128;

fehler_SAF = xor(dataIN10Folge, descrambledData_SAF);   % mit SAF
anzahl_fehler_SAF = sum(fehler_SAF);
BER_SAF = anzahl_fehler_SAF / 128;
