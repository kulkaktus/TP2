function [R,S]=poleplace(B,A,Hr,Hs,P)
    q = tf('q',0.04);
    d = 2;
    n_a = length(A)-1;
    n_b = length(B);

    n_r = n_a-1;
    n_s = n_b + d -1;
    R = zeros(length(n_r),1);
    S = zeros(length(n_s+1),1);
    % M1, M2 are left and right halves of sylvester matrix
    M1 = tril(toeplitz([A,zeros(1,n_b+d-1)]));
    M1 = M1(:, 1:n_b+d);
    M2 = tril(toeplitz([zeros(1,n_a+d-n_b+1), B, zeros(1,n_b-1)]));
    M2 = M2(:, 1:n_a);
    p = transpose([P,zeros(1, n_a+n_b+d-length(P))]);
    x = M\p;
    
    S = x(1:n_s+1);
    R = x(n_s+2:end);

end