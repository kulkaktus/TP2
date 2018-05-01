function [R,S]=poleplace(B,A,Hr,Hs,P)
    Ts = 0.04;
    d = 2;
    
    n_a = length(A)-1;
    n_b = length(B);
    n_hr = length(Hr)-1;
    n_hs = length(Hs)-1;
    n_a1 = n_a+n_hs;
    n_b1 = n_b+n_hr;
    
    n_s1 = n_b+n_hr+d-1;
    
    %Convert to TFs
    A = tf(A, 1, Ts, 'Variable', 'q^-1');
    B = tf(B, 1, Ts, 'Variable', 'q^-1');
    Hr = tf(Hr, 1, Ts, 'Variable', 'q^-1');
    Hs = tf(Hs, 1, Ts, 'Variable', 'q^-1');
    
    A1 = A*Hs;
    B1 = B*Hr;
    
    A1_coefs = cell2mat(A1.numerator);
    B1_coefs = cell2mat(B1.numerator);
    
    % M1, M2 are left and right halves of sylvester matrix
    M1 = tril(toeplitz([A1_coefs,zeros(1,n_b1+d-1)]));
    M1 = M1(:, 1:n_b1+d);
    M2 = tril(toeplitz([zeros(1,n_a1+d-n_b1+1), B1_coefs, zeros(1,n_b1-1)]));
    M2 = M2(:, 1:n_a1);
    M = [M1, M2];
    p = ([P,zeros(1, n_a1+n_b1+d-length(P))])';
    x = M\p;
    
    S1 = x(1:n_s1+1)';
    R1 = x(n_s1+2:end)';
    
    S1 = tf(S1, 1, Ts, 'Variable', 'q^-1');
    R1 = tf(R1, 1, Ts, 'Variable', 'q^-1');
    
    S = S1*Hs;
    R = R1*Hr;
    
    S = cell2mat(S.numerator);
    R = cell2mat(R.numerator);
end