function output = scrambler(seed, mod_modus, N)

    if mod_modus == 0
        numBits = 48*N-16;
    elseif mod_modus == 1
        numBits = 96*N-16;
    elseif mod_modus == 2
        numBits = 192*N-16;
    elseif mod_modus == 3
        numBits = 288*N-16;
    end

    %dataIn = randi([0 1], 1, numBits);
    dataIn= repmat([1, 0], 1, numBits/2);

    scrambledData = zeros(1, numBits);
    for i = 1:numBits
        feedbackBit = xor(seed(7), seed(4));
        scrambledData(i) = xor(dataIn(i), feedbackBit);
        seed = [feedbackBit seed(1:6)];
    end

    output.scrambledData = scrambledData;
    output.numBits = numBits+16;

end