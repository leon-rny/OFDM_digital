function r_pilot = pilot(data)
    insertPositions = [12, 26, 40];
    originalVecPos = 1;

    r_pilot = zeros(1, 64);
   
    for i = 7:59
        if ismember(i, insertPositions)
            r_pilot(:, i) = 1;
        elseif ismember(i, 54)
            r_pilot(:, i) = -1;
        elseif ismember(i, 33)
            r_pilot(:, i) = 0;
        else  
            r_pilot(:, i) = data(:, originalVecPos);
            originalVecPos = originalVecPos + 1;
        end
    end
end 