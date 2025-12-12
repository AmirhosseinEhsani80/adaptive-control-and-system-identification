function u_f = filter_signal(F, u, start)
    % with considering the initial conditions
    u_f = F * [u(start:-1:max(start-(length(F)-1),1)), zeros(1, 1 -(start-(length(F)-1)))].';
end