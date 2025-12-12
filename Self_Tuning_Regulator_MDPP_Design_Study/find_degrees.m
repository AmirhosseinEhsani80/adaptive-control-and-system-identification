function [deg_R, deg_Ao] = find_degrees(Am, A, Ao, Bminus, Bplus, BmPrime)
    syms q; 
    assert(Bplus(1)==1, "B plus is not monic");
    Ao_poly = poly2sym(Ao,q);
    A_poly = poly2sym(A,q);
    Am_poly = poly2sym(Am,q);
    Bminus_poly = poly2sym(Bminus,q);
    Bplus_poly = poly2sym(Bplus,q);
    BmPrime_poly = poly2sym(BmPrime,q);
    Bm_poly = BmPrime_poly * Bminus_poly;
    
    deg_Ao = polynomialDegree(A_poly,q)-polynomialDegree(Bplus_poly,q)-1;
    
%     assert(polynomialDegree(Ao_poly,q)==deg_A0, "degree of Ao is not correct");
    
    Ac_poly = Ao_poly * Am_poly * Bplus_poly;
    deg_R = polynomialDegree(Ac_poly,q) - polynomialDegree(A_poly,q);
    deg_S = deg_R;
    deg_Rprime = deg_R - polynomialDegree(Bplus_poly);

    
    

end