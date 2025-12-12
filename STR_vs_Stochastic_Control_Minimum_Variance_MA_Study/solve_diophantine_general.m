function [alpha, beta] = solve_diophantine_general(A, B, D, delay)
    n=length(A)-1;
    m=length(B)-1+delay; 
    k=length(D)-1; 
    B_extended=[zeros([1,n-m]),B,zeros([1,delay])];
    D_extended=[zeros([1,(2*n-1)-k]),D]; 
    E =0*ones([2*n,2*n]);
    for i=1:2*n
       if i<=n
        E(i,1:n)=[A(end-i+1:1:end),zeros([1,n-i])];
        E(i,n+1:2*n)=[B_extended(end-i+1:1:end),zeros([1,n-i])];
       else
        E(i,1:n)=[zeros([1,-(n-i+1)]),A(1:end-i+(n))]; 
        E(i,n+1:2*n)=[zeros([1,-(n-i+1)]),B_extended(1:end-i+(n))];
       end
    end
    D_to_solve=transpose(fliplr(D_extended)); 
    M=inv(E)*D_to_solve;
    alpha=fliplr(M(1:n).');
    beta=fliplr(M(n+1:2*n).');
end