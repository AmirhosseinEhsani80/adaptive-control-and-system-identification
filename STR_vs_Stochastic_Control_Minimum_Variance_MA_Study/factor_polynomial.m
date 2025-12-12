function [poly_plus, poly_minus] = factor_polynomial(polynomial)
    poly_roots = roots(polynomial);
    plus = poly_roots(abs(poly_roots)<=1); %inside circle
    minus =  poly_roots(abs(poly_roots)>1);
    

    poly_plus = poly(plus);
    poly_minus = poly(minus);

    poly_minus = poly_minus * polynomial(1)  * poly_plus(1);
    poly_plus = poly_plus / poly_plus(1);
 
end

    