function [R_solved,S_solved,T] =Diophantine(Am, A, Ao, Bminus, Bplus, BmPrime)
    syms q; 
    assert(Bplus(1)==1, "B plus is not monic");
    Ao_poly = poly2sym(Ao,q);
    A_poly = poly2sym(A,q);
    Am_poly = poly2sym(Am,q);
    Bminus_poly = poly2sym(Bminus,q);
    Bplus_poly = poly2sym(Bplus,q);
    BmPrime_poly = poly2sym(BmPrime,q);
    Bm_poly = BmPrime_poly * Bminus_poly;
    
    deg_A0 = polynomialDegree(A_poly,q)-polynomialDegree(Bplus_poly,q)-1;
    
    assert(polynomialDegree(Ao_poly,q)==deg_A0, "degree of Ao is not correct");
    
    Ac_poly = Ao_poly * Am_poly * Bplus_poly;
    deg_R = polynomialDegree(Ac_poly,q) - polynomialDegree(A_poly,q);
    deg_S = deg_R;
    deg_Rprime = deg_R - polynomialDegree(Bplus_poly);
    
    Rprime = sym("r", [1, deg_Rprime]);
    Rprime = [1,Rprime];
    
    S = sym("s", [1, deg_S+1]);
    
    Rprime_poly = poly2sym(Rprime,q);
    S_poly = poly2sym(S,q);
    
    LHS = A_poly * Rprime_poly + Bminus_poly * S_poly;
    RHS = Ao_poly * Am_poly;
    assert(polynomialDegree(LHS,q)==polynomialDegree(RHS,q), "LHS and RHS are not in same degree");
    
    LHS_coef = coeffs(LHS, q, 'all');
    RHS_coef = coeffs(RHS, q, 'all');
    
    
    if polynomialDegree(Bminus_poly * S_poly,q) < polynomialDegree(A_poly * Rprime_poly,q)
        results = double(vpa(struct2array(solve(LHS_coef(2:end)-RHS_coef(2:end)))));
    else
        results = double(vpa(struct2array(solve(LHS_coef-RHS_coef))));
    end
    
    Rprime_solved = results(1:length(Rprime)-1);
    S_solved = results(length(Rprime):length(Rprime) + length(S)-1);
    
    Rprime_solved = [1, Rprime_solved];
    Rprime_solved_poly = poly2sym(Rprime_solved,q);
    
    S_solved_poly = poly2sym(S_solved,q);
    R_solved_poly = Rprime_solved_poly * Bplus_poly;
    
    R_solved = double(vpa(coeffs(R_solved_poly, q, 'all')));
    
    T_poly = Ao_poly * BmPrime_poly;
    T = double(vpa(coeffs(T_poly, q, 'all')));

end