function data_mapped = mapping(mod_modus, data_mapped)

    if mod_modus == 0                           % BPSK
        result = zeros(48, 1);
        result(:, 1) = data_mapped(1:1:end);    % I-Teil
        result(data_mapped < 0.5) = -1;
        data_mapped = result.';

    elseif mod_modus == 1                       % QPSK
        result = zeros(48, 2);
        result(:, 1) = data_mapped(1:2:end);    % I-Teil
        result(:, 2) = data_mapped(2:2:end);    % Q-Teil
        result(result < 0.5) = -1;
        result_complex = result(:, 1) + 1i * result(:, 2);
        data_mapped = 1/sqrt(2) .* result_complex;
        data_mapped = data_mapped.';
    
        elseif mod_modus == 2                       % 16-QAM
            result = zeros(1, size(data_mapped, 2)./2); 
        
            for i = 1:2:size(data_mapped, 2) 
                pair = data_mapped(i:i+1);  
                if isequal(pair, [0, 0])
                    result((i+1)/2) = -3;  
                elseif isequal(pair, [0, 1])
                    result((i+1)/2) = -1;  
                elseif isequal(pair, [1, 1])
                    result((i+1)/2) = 1;   
                elseif isequal(pair, [1, 0])
                    result((i+1)/2) = 3;   
                end
            end
            result_2 = zeros(48,2);
        
            result_2(:, 1) = result(1:2:end);    % I-Teil
            result_2(:, 2) = result(2:2:end);    % Q-Teil
            result_complex = result_2(:, 1) + 1i * result_2(:, 2);
            data_mapped = 1/sqrt(10).*result_complex;
            data_mapped = data_mapped.';
        
    elseif mod_modus == 3                       %64-QAM
        result = zeros(1, size(data_mapped, 2)./8);  
        for i = 1:3:size(data_mapped, 2)
            pair = data_mapped(i:i+2);  
            
            index = (i-1)/3 + 1;
        
            if isequal(pair, [0, 0, 0])
                result(index) = -7;  
            elseif isequal(pair, [0, 0, 1])
                result(index) = -5;  
            elseif isequal(pair, [0, 1, 1])
                result(index) = -3;   
            elseif isequal(pair, [0, 1, 0])
                result(index) = -1;   
            elseif isequal(pair, [1, 1, 0])
                result(index) = 1;  
            elseif isequal(pair, [1, 1, 1])
                result(index) = 3;   
            elseif isequal(pair, [1, 0, 1])
                result(index) = 5; 
            elseif isequal(pair, [1, 0, 0])
                result(index) = 7; 
            end
        end
    
        result_2 = zeros(48,2);
    
        result_2(:, 1) = result(1:2:end);    % I-Teil
        result_2(:, 2) = result(2:2:end);    % Q-Teil
        result_complex = result_2(:,1) + 1i * result_2(:,2);
        data_mapped = 1/sqrt(42).*result_complex;
        data_mapped = data_mapped.';
    end 
end