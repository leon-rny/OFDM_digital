clear; close all;
%load('Tx_data_mod0_theorie.mat')
load('Tx_data_offset_mod1.mat');
N = 12;             % Anzahl der OFDM-Symbole also mit Service+Data
mod_modus = 0;
f_s = 20e6;
N_P = 16;
y = [zeros(1,500) r_data];

%% Empfangsvorgang
%{
OpenFPGA;
Open_RFFE_UART;
SetRxGain(lna_mode, g);
RFFE_SetLO_Frequency(f);
LOf = RFFE_ReadLO_Frequency();
Ld = RFFE_ReadLO_LockDetect();
%r_data = RecordSignal(16384,200);
%}

%% Frame Start Detector
for n = 1:length(y)-32
    tmp = 0;
    for k = 0:N_P-1
        tmp = tmp + (y(n+k)) * conj(y(n+k+N_P));        
    end
    y_scx(n) = tmp;
end
index_start = find(abs(y_scx) > max(abs(y_scx))/2);
index_start = index_start(1);

y_scx_plot = y_scx;
y_scx = y_scx(index_start+7:index_start+129);
%% Frequenzoffset Detector
% y_max = max(abs(y_scx));
angles1 = angle(y_scx);
angle1 = mean(angles1);

df = -angle1*f_s/(2*pi*N_P);


%% Frequenzoffset Kompensation
for n = 1:length(y)
    y(n) = y(n) * exp(1j * (angle1+0.01)*n/N_P);
end

% for n = 1:length(y)
%     y(n) = y(n) * exp(1j *0.01*n/N_P);
% end








%%
% %% Frequenzoffset Detector
% % y_max = max(abs(y_scx));
% y_aa = angle(y_scx);
% y_aa1 = mean(y_aa);
% %y_avg = mean(y_scx);
% % y_avg = mean(y_scx(abs(y_scx) > (y_max / 2.2))); % Mittelwertbildung
% df = -angle(y_avg)*f_s/(2*pi*N_P);
% 
% 
% %% Frequenzoffset Kompensation
% % for n = 1:length(y)
% %     y(n) = y(n) * exp(1j * (angle(y_avg)+0.01)*n/N_P);
% % end
% 
% for n = 1:length(y)
%     y(n) = y(n) * exp(1j *0.01*n/N_P);
% end













%% Preamble und Guard Band Entfernung
y_data = y(index_start+327-4:end-4);                    % Entfernung Short und Long Preamble
l = length(y_data);
%r = rem(l,80);                                     % was nicht passt, wird passend gemacht
y_only_OFDM = reshape(y_data, [80, N+1]);           % Angeordnete OFDM Symbole 
y_only_OFDM(1:16,:) = [] ;                          % Guard ausschneiden

%% Fouriertrafo
y_data_fft = [];
for i = 1:N+1
    tmp_fft = fft(y_only_OFDM(:,i), 64);
    tmp_fft_shifted = fftshift(tmp_fft);
    y_data_fft = [y_data_fft tmp_fft_shifted];
end

%% Nullen Entfernung
y_reshaped = [y_data_fft(7:32,:); y_data_fft(34:59,:)];

%% Channel Estimator
y_long_1_raw = y(index_start+199-4:index_start+262-4);
y_long_1_fft = fft(y_long_1_raw, 64);
y_long_1_shift = fftshift(y_long_1_fft);
y_long_1 = [y_long_1_shift(7:32) y_long_1_shift(34:59)];

y_long_2_raw = y(index_start+263-4:index_start+326-4);
y_long_2_fft = fft(y_long_2_raw, 64);
y_long_2_shift = fftshift(y_long_2_fft);
y_long_2 = [y_long_2_shift(7:32) y_long_2_shift(34:59)];

L = [1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1];

for i = 1:52
    H(i) = 1/2 * ((y_long_1(i) + y_long_2(i)) / L(i));
end

figure(6);
plot(unwrap(angle(H)));
hold on;
plot(abs(H));
grid on;

%% Channel Equalizer

for n = 1:N+1
    for i = 1:52
        D(i,n) = y_reshaped(i,n)/H(i);
    end 
end

%% Phase Derotation
y_pilot = D([6, 20, 33, 47], :);
D([6, 20, 33, 47], :) = [];
y_pilot(4, :) = y_pilot(4, :) * -1;
H_phase = H([6, 20, 33, 47]);

% CPE berechnen 13x1
for i = 1:4
    sum_H = (abs(H_phase(i)).^2);
end
sum_H = sum_H*4;

for n = 1:N+1
    sum_end = 0;
    for i = 1:4
        sumHH = abs(H_phase(i))^2*y_pilot(i,n);
        sum_end = sum_end + sumHH;
    end
    CPE(n) = 1/sum_H .* sum_end; 
end

% CPE / Datensignal 
for n = 1:N+1
    for i = 1:48
        y_corrected(i, n) = D(i, n) / CPE(n);
    end
end


y_corrected_plot = y_corrected;
y_corrected_plot(:,1) = [];

%% Demapping
y_binary = [];
for i = 2:N+1
    tmp = demapping(mod_modus, y_corrected(:,i));
    y_binary = [y_binary; tmp];
end
y_binary = y_binary .';

%% Descramler
y_wo_service = y_binary(17:end);
descrambledData = descrambler(y_wo_service, mod_modus, N);

%% Plotting


% % Schmidl-Cox-Korrelators
% f1 = figure(1);
% 
% % plot(real(y_scx_plot));
% % hold on;
% % plot(imag(y_scx_plot));
% plot(abs(y_scx_plot));
% ylabel('Amplitude');
% xlabel('Sample');
% %legend('Imaginärteil','Realteil');
% grid on;


% % Raw IQ-Ebene
% f2 = figure(2);
% for i = 1:N
%     plot(real(y_reshaped(:,i)), imag(y_reshaped(:,i)), 'o', 'Color', 'r');
%     hold on;
%     grid on;
%     axis equal;
% end
% xlabel('In-Phase (I)')
% ylabel('Quadratur (Q)')
% 

% nach Frequenz-Offset-Entfernen IQ-Ebene
f3 = figure(3);
for i = 1:N
    plot(real(y_reshaped(:,i)), imag(y_reshaped(:,i)), 'o', 'Color', 'r');
    hold on;
    grid on;
    axis equal;
end
xlabel('In-Phase (I)')
ylabel('Quadratur (Q)')


% nach Kanalschätzung 
f4 = figure(4);

for i = 1:N
    plot(real(D(:,i)), imag(D(:,i)), 'o', 'Color', 'r');
    hold on; 
    grid on;
    axis equal;
end

for i = 1:N
    h1 = plot(real(CPE(:,i)), imag(CPE(:,i)), 'o', 'Color', 'b', 'LineWidth', 10);
end


xlabel('In-Phase (I)');
ylabel('Quadratur (Q)');
legend(h1, 'CPE'); 


% nach allem IQ-Ebene
f5 = figure(5);
for i = 1:N
    plot(real(y_corrected_plot(:,i)), imag(y_corrected_plot(:,i)), 'o', 'Color', 'r');
    hold on;
    grid on;
    axis equal;
end
xlabel('In-Phase (I)')
ylabel('Quadratur (Q)')


%exportgraphics(f1,'Plots/Rx_Theorie/3_2_Schmidl_Cox_Korrelator_thoerie.pdf','ContentType','vector');

% exportgraphics(f2,'Plots/Rx_Theorie/3_2_Theorie_BPSK_raw.pdf','ContentType','vector');
% exportgraphics(f3,'Plots/Rx_Theorie/3_2_Theorie_BPSK_nach_Frequenz_Offset_Entfernen.pdf','ContentType','vector');
% exportgraphics(f4,'Plots/Rx_Theorie/3_2_Theorie_BPSK_nach_Kanalschätzung.pdf','ContentType','vector');
% exportgraphics(f5,'Plots/Rx_Theorie/3_2_Theorie_BPSK_nach_Derotation.pdf','ContentType','vector');


% exportgraphics(f2,'Plots/Rx_Theorie/3_2_Theorie_QPSK_raw.pdf','ContentType','vector');
% exportgraphics(f3,'Plots/Rx_Theorie/3_2_Theorie_QPSK_nach_Frequenz_Offset_Entfernen.pdf','ContentType','vector');
% exportgraphics(f4,'Plots/Rx_Theorie/3_2_Theorie_QPSK_nach_Kanalschätzung.pdf','ContentType','vector');
% exportgraphics(f5,'Plots/Rx_Theorie/3_2_Theorie_QPSK_nach_Derotation.pdf','ContentType','vector');


% exportgraphics(f2,'Plots/Rx_Theorie/3_2_Theorie_16QAM_raw.pdf','ContentType','vector');
% exportgraphics(f3,'Plots/Rx_Theorie/3_2_Theorie_16QAM_nach_Frequenz_Offset_Entfernen.pdf','ContentType','vector');
% exportgraphics(f4,'Plots/Rx_Theorie/3_2_Theorie_16QAM_nach_Kanalschätzung.pdf','ContentType','vector');
% exportgraphics(f5,'Plots/Rx_Theorie/3_2_Theorie_16QAM_nach_Derotation.pdf','ContentType','vector');

% exportgraphics(f2,'Plots/Rx_Theorie/3_2_Theorie_64QAM_raw.pdf','ContentType','vector');
% exportgraphics(f3,'Plots/Rx_Theorie/3_2_Theorie_64QAM_nach_Frequenz_Offset_Entfernen.pdf','ContentType','vector');
% exportgraphics(f4,'Plots/Rx_Theorie/3_2_Theorie_64QAM_nach_Kanalschätzung.pdf','ContentType','vector');
% exportgraphics(f5,'Plots/Rx_Theorie/3_2_Theorie_64QAM_nach_Derotation.pdf','ContentType','vector');





