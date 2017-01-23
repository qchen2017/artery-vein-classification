function displayGraphLengths(T, true_length, straight_length)

    I_trueLength = zeros(T.w, T.l);
    I_straightLength = zeros(T.w, T.l);
    for i = 1 : length(T.link)

        I_trueLength(T.link(i).point) = true_length(i);
        I_straightLength(T.link(i).point) = straight_length(i);
        
    end
    figure, imagesc(I_trueLength);
    figure, imagesc(I_straightLength);

end