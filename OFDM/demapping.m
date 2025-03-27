function data_demapped = demapping(mod_modus, data)

    if mod_modus == 0                           % BPSK
        result = zeros(48, 1);
        result(data < 0) = 0;
        result(data > 0) = 1;
        data_demapped = result.';

    elseif mod_modus == 1                       % QPSK
        data = sqrt(2) .* data;
        % Initialisierung der Ausgabematrix
        data_demapped = zeros(length(data), 2);
        
        % Definition der Zuordnung entsprechend der 64-QAM Tabelle
        I_mapping = [-1, 1];
        Q_mapping = I_mapping; % Identisch für I und Q
        bit_mapping = {'0', '1'};
        
        for i = 1:length(data)
            % Extrahieren des Real- und Imaginärteils
            I = real(data(i));
            Q = imag(data(i));
            
            % Zuordnung zu Binärcodes
            [~, I_index] = min(abs(I_mapping - I));
            [~, Q_index] = min(abs(Q_mapping - Q));
            
            % Umwandlung in Binärstring und dann in Binärvektor
            binaryStr = strcat(bit_mapping{I_index}, bit_mapping{Q_index});
            binaryVec = arrayfun(@(c) str2double(c), binaryStr);
            
            tempMatrix(i, :) = binaryVec;

        end
        data_demapped = reshape(tempMatrix', [], 1);

    elseif mod_modus == 2                       % 16-QAM
        data = sqrt(10) .* data;
        % Initialisierung der Ausgabematrix
        data_demapped = zeros(length(data), 4);
        
        % Definition der Zuordnung entsprechend der 64-QAM Tabelle
        I_mapping = [-3, -1, 1, 3];
        Q_mapping = I_mapping; % Identisch für I und Q
        bit_mapping = {'00', '01', '11', '10'};
        
        for i = 1:length(data)
            % Extrahieren des Real- und Imaginärteils
            I = real(data(i));
            Q = imag(data(i));
            
            % Zuordnung zu Binärcodes
            [~, I_index] = min(abs(I_mapping - I));
            [~, Q_index] = min(abs(Q_mapping - Q));
            
            % Umwandlung in Binärstring und dann in Binärvektor
            binaryStr = strcat(bit_mapping{I_index}, bit_mapping{Q_index});
            binaryVec = arrayfun(@(c) str2double(c), binaryStr);
            
            tempMatrix(i, :) = binaryVec;

        end
        data_demapped = reshape(tempMatrix', [], 1);

    elseif mod_modus == 3                   % 64-QAM
        data = sqrt(42) .* data;
        % Initialisierung der Ausgabematrix
        data_demapped = zeros(length(data), 6);
        
        % Definition der Zuordnung entsprechend der 64-QAM Tabelle
        I_mapping = [-7, -5, -3, -1, 1, 3, 5, 7];
        Q_mapping = I_mapping; % Identisch für I und Q
        bit_mapping = {'000', '001', '011', '010', '110', '111', '101', '100'};
        
        for i = 1:length(data)
            % Extrahieren des Real- und Imaginärteils
            I = real(data(i));
            Q = imag(data(i));
            
            % Zuordnung zu Binärcodes
            [~, I_index] = min(abs(I_mapping - I));
            [~, Q_index] = min(abs(Q_mapping - Q));
            
            % Umwandlung in Binärstring und dann in Binärvektor
            binaryStr = strcat(bit_mapping{I_index}, bit_mapping{Q_index});
            binaryVec = arrayfun(@(c) str2double(c), binaryStr);
            
            tempMatrix(i, :) = binaryVec;

        end
        data_demapped = reshape(tempMatrix', [], 1);
    end
 end
