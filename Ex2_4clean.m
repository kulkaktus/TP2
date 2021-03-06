%% 2.4 Q-parametrization
load('TorMod.mat');
Ts = 0.04;
A = G3.f;
B = G3.b;
p1 = -1.7689; 
p2 = 0.7939;
P = [1 p1 p2];

Hs = [1 -1];
Hr = [1 1]; 

% Set auxiliary poles
a = 0.16;
poles_aux = [a,a,a,a,a,a,a,a,a,a];
coefs = poly(poles_aux);
P_new = conv(P, coefs);

[R, S] = poleplace(B, A, Hr, Hs, P_new) ;
P_end = conv(A,S) + conv(B,R);

B = B(2:end); %Separate delay and B
d = [0 1];

M_m = 0.4;
U_max = 56.2; %35 dB  = 56.2

% Create R and S functions
Q_length = 8;
R0 = [R, zeros(1, Q_length)];
S0 = [S, zeros(1, Q_length)];
R_new = @(Q) R0 + conv(A, conv(Hr, conv(Hs, Q)));
S_new = @(Q) S0 - conv(d, conv(B, conv(Hs, conv(Hr, Q))));

%Set inequality (c) and equality (ceq) constraints
c = @(Q) [norm(M_m*S_new(Q), Inf) - 1;
     norm(tf(conv(A, R_new(Q)), P_end, ...
            Ts, 'variable', 'z^-1'), Inf) - U_max];
ceq = @(Q) [];
Nonlincon = @(Q)deal(c(Q), ceq(Q));

% Define modulus margin and optimize
Mod_marg = @(Q) norm(tf(S_new(Q), 1 ,Ts, 'variable', 'z^-1'), Inf)^(-1);
Q_opt = fmincon(Mod_marg, zeros(1, Q_length),[],[],[],[], ...
                [-Inf,-Inf],[Inf,Inf], Nonlincon);

% Evaluate R,S,T using optimal Q
R_final = R_new(Q_opt);
S_final = S_new(Q_opt);
T = sum(R_final);

%% Check constraints
CL = tf(conv(T,G3.b), P_end, Ts,'variable','z^-1'); %Closed loop system
U = tf(conv(A,R_final), P_end, Ts, 'variable', 'z^-1'); %Input sensitivity function
input = tf(conv(T,A), conv(A,S_final) + conv(conv(d, B), R_final),...
          Ts,'variable','z^-1'); 
      
figure(1);  
step(input); 
figure(2);
step(CL); 
figure(3);
bodemag(U);
MM = inv(norm(S_final,inf)) %Modulus margin