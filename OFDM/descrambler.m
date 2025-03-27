function descrambledData = descrambler(data, mod_modus, N)
    %N = N-1;
    if mod_modus == 0
        numBits = 48*N-16;
    elseif mod_modus == 1
        numBits = 96*N-16;
    elseif mod_modus == 2
        numBits = 192*N-16;
    elseif mod_modus == 3
        numBits = 288*N-16;
    end
   
    seed= [1 0 1 0 1 0 1];
    descrambledData = zeros(1, numBits);
    for i = 1:numBits
        feedbackBit = xor(seed(7), seed(4));
        descrambledData(i) = xor(data(i), feedbackBit);
        seed = [feedbackBit seed(1:6)];
    end

end