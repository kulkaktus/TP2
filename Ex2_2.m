%% 2.2
load('TorMod.mat');
A = G3.f;
B = G3.b;
p1 = -1.7689; 
p2 = 0.7939;
P = [1 p1 p2];


Hs = [1 -1];
Hr = [1]; 
[R_old, S_old] = poleplace(B, A, Hr, Hs, P);

P_end_old = conv(A,S_old) + conv(B,R_old);
T_old = sum(R_old)


%% Ex2_3
Ts = 0.04;
CL = tf(conv(T_old,B), P_end_old, Ts,'variable','z^-1')
figure(2)
step(CL)
stepinfo(CL)

K = tf(R,S, Ts,'variable','z^-1')
OL = tf(G3*K);
%figure(2)
%nyquist(OL) % 3 encirclements around -1 corresponds to 3 unstable poles 



figure(3)
input_old = tf(conv(T_old,A), conv(A,S_old) + conv(B, R_old), Ts,'variable','z^-1');
step(input)
figure(4)

U=tf(conv(A,R), P_end, Ts, 'variable', 'z^-1');
bodemag(U)


bode(U)
figure(5)
bode(OL)

S3 = feedback(1,K*G3);
M_m3 = inv(norm(S3,inf)) % 0.1566
%% 2.4
load('TorMod.mat');
Ts = 0.04;
A = G3.f;
B = G3.b;
p1 = -1.7689; 
p2 = 0.7939;
P = [1 p1 p2];

Hs = [1 -1];
Hr = [1 1]; 
a = 0.2;
nq = 9;


Margins = zeros(1, 3);
hail_marys = zeros(1, 4); %These satisfy all constraints and specs
index = 1;  
for a = 0.16
    for nq = 8
        B = G3.B;
        poles_aux = [a,a,a,a,a,a,a,a,a,a];
        coefs = poly(poles_aux);
        P_new = conv(P, coefs);

        [R, S] = poleplace(B, A, Hr, Hs, P_new) ;
        P_end = conv(A,S) + conv(B,R);

        B = B(2:end); %Separate delay and B
        q_delay = [0 1];

        M_m = 0.4;
        U_max = 56.2; %35 dB  = 56.2

        R = [R, zeros(1, nq)];
        S = [S, zeros(1, nq)];

        R_new = @(Q) R + conv(A, conv(Hr, conv(Hs, Q)));
        S_new = @(Q) (S - conv(q_delay, conv(B, conv(Hs, conv(Hr, Q)))));
        
        K_new = @(Q) tf(R_new(Q),S_new(Q), Ts,'variable','z^-1');
        output_sensitivity = @(Q) feedback(1,K_new(Q)*G3);
        Mod_marg = @(Q) norm(output_sensitivity(Q), Inf)^(-1);

        %Set inequality (c) and equality (ceq) constraints
        c = @(Q) [norm(M_m*output_sensitivity(Q), Inf) - 1;
             norm(tf(conv(A, R_new(Q)), P_end, Ts, 'variable', 'z^-1'), Inf) - U_max];
        ceq = @(Q) [];
        Nonlincon = @(Q)deal(c(Q), ceq(Q));


        Q_opt = fmincon(Mod_marg, zeros(1, nq),[],[],[],[],[-Inf,-Inf],[Inf,Inf], Nonlincon);

        R_final = R_new(Q_opt);
        S_final = S_new(Q_opt);
        T = sum(R_final);

        K_final = tf(R_final,S_final, Ts,'variable','z^-1');
        Sens_out_new = feedback(1,K_final*G3);
        MM_new = Mod_marg(Q_opt)
     
        Margins(index, :) = [MM_new, a, nq];
        index = index + 1
        CL = tf(conv(T,G3.b), P_end, Ts,'variable','z^-1'); 
        cake = stepinfo(CL)
    end
end
hail_marys
%% Check constraints
CL = tf(conv(T,G3.b), P_end, Ts,'variable','z^-1'); 
U=tf(conv(A,R_final), P_end, Ts, 'variable', 'z^-1');

input= tf(conv(T,A), conv(A,S_final) + conv(conv(q_delay, B), R_final), Ts,'variable','z^-1'); figure(1);  
step(input); 
stepinfo(CL) 
%figure(3);
%bodemag(U);


%% Plots
hold on
figure(1)
subplot(2,2,1)
CL_old = tf(conv(T_old,G3.b), P_end_old, Ts,'variable','z^-1'); 
CL = tf(conv(T,G3.b), P_end, Ts,'variable','z^-1'); 
title("Tracking step response, original and improved controller")
xlabel( "Time [s]")
legend("Original", "Improved")
step(CL_old, CL)


U=tf(conv(A,R), P_end, Ts, 'variable', 'z^-1');
figure(3)
step(U)% Does not fullfill the criterium


M_m = 0.4;
U_max = 100; %35 dB  = 56.2
q_delay = [0 1];

R = [R, 0, 0];
S = [S, 0, 0, 0];
new_R = @(Q) R + conv(A,conv(Hr, conv(Hs, Q)));
new_S = @(Q) (S - conv(q_delay, conv(B, conv(Hs, conv(Hr, Q)))));

c = @(Q) [norm(M_m*(S + conv(q_delay, conv(B, conv(Hs, conv(Hr, Q))))), Inf) - 1;
     norm(deconv(conv(A, new_R(Q)), P_end), Inf) - U_max];
ceq = @(Q) [];
zero = @(Q) 0;
Mod_marg = @(Q) -(norm(S + conv(q_delay, conv(B, conv(Hs, conv(Hr, Q)))), Inf))^(-1);
Nonlincon = @(Q)deal(c(Q), ceq(Q));
Q_opt = fmincon(zero, [1,1], [], [], [],[], [-Inf,-Inf], [Inf,Inf], Nonlincon);

final_R = new_R(Q_opt);
final_S = new_S(Q_opt);


subplot(2,2,2)
B = G3.B;
B = B(2:end);


title("Control signal, original and improved controller")
xlabel("Time [s]")
ylabel("u")
legend("Original", "Improved")
step(input_old, input)


subplot(2,2,3)
Sens_out = feedback(1,K*G3);
K_final = tf(R_final,S_final, Ts,'variable','z^-1');
Sens_out_new = feedback(1,K_final*G3);
MM_old = norm(Sens_out, Inf)^-1
MM_new = (norm(Sens_out_new, Inf))^(-1)

title("Output sensitivity function, original and improved controller")
legend("Original","Improved")
bodemag(Sens_out, Sens_out_new)


subplot(2,2,4)
U=tf(conv(A,R), P_end_old, Ts, 'variable', 'z^-1');
U_new=tf(conv(A,R_final), P_end, Ts, 'variable', 'z^-1');
title("Input sensitivity function, original and improved controller")
legend("Original","Improved")
bodemag(U, U_new)
