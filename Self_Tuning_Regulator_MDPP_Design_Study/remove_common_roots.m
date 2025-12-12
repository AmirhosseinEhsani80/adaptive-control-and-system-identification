function [R, S] = remove_common_roots(BR,BS,precision)
    
    roots_BR = roots(BR);
    roots_BS = roots(BS);
    
    index_pairs = [];
    for i=1:length(roots_BS)
        same_roots_R = discretize(real(roots_BR), [real(roots_BS(i))-precision,real(roots_BS(i))+precision])==1 & discretize(imag(roots_BR), [imag(roots_BS(i))-precision,imag(roots_BS(i))+precision])==1;
        index_R = find(same_roots_R==1);
        if ~isempty(index_R)
            index_S = i;
            index_pairs = [index_pairs, [index_R; index_S]]; % first row is for R and second row is for S
        end
    end
    
    if ~isempty(index_pairs)
        common_R = poly(roots_BR(index_pairs(1, :)));
        common_S = poly(roots_BS(index_pairs(2, :)));
        R = deconv(BR, common_R);
        S = deconv(BS, common_S);
    else
        R = BR;
        S = BS;
    end

end