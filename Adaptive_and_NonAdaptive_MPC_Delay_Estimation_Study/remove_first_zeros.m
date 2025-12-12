function P = remove_first_zeros(P)  
    index = -1;
    for i=1:length(P)
        if abs(P(i)) > 0.0000000001
            index = i;
            break
        end
    end
    if index ~=-1
        P = P(index:end);
    end
    if P == zeros(1, length(P))
        P = 0;
    end
end