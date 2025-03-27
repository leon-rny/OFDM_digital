%% GUI aufrufen
clear; close all;
OpenFPGA
Open_RFFE_UART

%% Parameter
p_Tx = 20;
f_Tx = 2500;
RecordLength = 3820;
Symbol_Length = 1;
T_D = 1e-6 * Symbol_Length;
f_s = 20e6;
M = T_D/f_s;

%% allgemeine Sendereinstellungen

OpenFPGA
Open_RFFE_UART

SetTxPower(p_Tx);
RFFE_SetLO_Frequency(f_Tx);

%% Datenquelle und Scrambler
lfsrState1 = [1 0 1 0 1 0 1];
numBits = 128;
% dataIn = randi([0 1], 1, numBits);
dataIn= repmat([1, 0], 1, 64);
scrambledData = zeros(1, numBits);

for i = 1:numBits
    feedbackBit = xor(lfsrState1(7), lfsrState1(4));
    scrambledData(i) = xor(dataIn(i), feedbackBit);
    lfsrState1 = [feedbackBit lfsrState1(1:6)];
end

%% Pr�ambelsequenz
% Parameter definieren
N = 63; % L�nge der m-Sequenz
registerLength = 6; % L�nge des Registers

% Initialisieren des LFSR mit einem zuf�lligen Seed, der nicht Null ist
% Der Seed sollte eine zuf�llige bin�re Sequenz der L�nge 6 sein
seed = [1 1 0 1 1 0];
% Initialisieren des Registers mit dem Seed
register = seed;

% Definieren der Tap-Positionen basierend auf dem Polynom x^6 + x + 1
% MATLAB-Indizes beginnen bei 1, daher entsprechen die Positionen 6 und 1
taps = [1, registerLength]; % Die Taps sind an den Positionen 1 (x^6) und 6 (x)

% Generieren der Pr�ambelsequenz
p = zeros(N, 1); % Initialisieren der Sequenz

for i = 1:N
    % Die Sequenz wird durch das letzte Element des Registers generiert
    p(i) = register(end);
    
    % Berechnen des n�chsten Bit (Feedback)
    feedback = mod(sum(register(taps)), 2);
    
    % Aktualisieren des Registers
    register = [feedback, register(1:end-1)];
end

trans_p = p.';
dataALL = [trans_p, scrambledData];

%% 2-ASK Modulation
pmod = 2*p - 1; % Umwandlung der Bitsequenz: 0 wird zu -1, 1 wird zu +1
scrambledData_mod = 2*scrambledData -1;
% �berpr�fen der periodischen Autokorrelationseigenschaften
PAKF = ifft(fft(pmod).*conj(fft(pmod)));
trans_pmod = pmod.';
dataALL = [trans_pmod, scrambledData_mod];
dataALL = dataALL.';

%% �bertabtastung
% Definieren des �berabtastungsfaktors
M = 20;

% Erstellen eines Delta-Impulses f�r die �berabtastung
deltaImpuls = [1; zeros(M-1, 1)]; % Ein Vektor mit einer '1' gefolgt von M-1 Nullen

% Anwenden der �berabtastung auf die 2-ASK modulierte Sequenz pmod
pmod_uberabgetastet = kron(dataALL, deltaImpuls);

%% FIR-Sendepulsfilter
hr = ones(1, M); % Impulsantwort des FIR-Filters / Rechteckpuls

% Anwenden des FIR-Filters auf die �berabgetastete Sequenz
pmod_gefiltert = filter(hr, 1, pmod_uberabgetastet);

%% DC-Offset
pmod_dc_offset = pmod_gefiltert + 1;

%% Hochsetzen auf niedrige Zwischenfrequenz
% Angenommene Parameter
f_IF = 500e3; % Zwischenfrequenz von 500 kHz
f_s = 20e6; % Angenommene Abtastrate von 1 MHz
N = length(pmod_dc_offset); % L�nge des Signals mit DC-Offset
n = 0:N-1; % Index-Vektor

% Erzeugen der komplexen Exponentialfunktion f�r die Hochsetzung
exp_signal = 250*exp(1j*2*pi*f_IF*n/f_s);

% Hochsetzen des Signals auf die Zwischenfrequenz
xdbb_if = pmod_dc_offset.' .* exp_signal;
% xdbb_if = xdbb_if + 250*randn(1, length(xdbb_if));

f_achse = f_s/RecordLength*(0:RecordLength-1);  % Frequenzachse
fx = f_achse - max(f_achse)/2;

pmod_fft = fft(pmod_dc_offset);                                                           % Skalierung hinzuf�gen /RecordLength
pmod_fft_shifted = abs(fftshift(pmod_fft));                                                % Verschieben des Spektrums

xdbb_if_fft = fft(xdbb_if);                                                           % Skalierung hinzuf�gen /RecordLength
xdbb_if_fft_shifted = abs(fftshift(xdbb_if_fft));                                                % Verschieben des Spektrums
