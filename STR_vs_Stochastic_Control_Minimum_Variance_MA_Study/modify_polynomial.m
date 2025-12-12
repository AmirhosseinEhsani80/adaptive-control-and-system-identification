function modifed_poly = modify_polynomial(polynomial)
    poly_roots = roots(polynomial);
    poly_plus = poly_roots(abs(poly_roots)<=1); %inside circle
    poly_minus =  poly_roots(abs(poly_roots)>1);
    poly_minus_modif = fliplr(poly(poly_minus));
    poly_plus_modif = poly(poly_plus);
    modifed_poly = conv(poly_minus_modif, poly_plus_modif);
    modifed_poly = modifed_poly/modifed_poly(1);
end