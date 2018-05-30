%% 2.2
load('TorMod.mat');
A = G3.f;
B = G3.b;
p1 = -1.7689;
p2 = 0.7939;
P = [1 p1 p2];


Hs = [1 -1];
Hr = [1]; 
[R, S] = poleplace(B, A, Hr, Hs, P);

P_end = conv(A,S) + conv(B,R);
T = sum(R)


%% Ex2_3
Ts = 0.04;
CL = tf(conv(T,B), P_end, Ts,'variable','z^-1')
figure(1)
step(CL)
stepinfo(CL)

K = tf(R,S, Ts,'variable','z^-1')
OL = tf(G3*K);
figure(2)
nyquist(OL) % 3 encirclements around -1 corresponds to 3 unstable poles 



figure(3)
input = tf(conv(T,A), conv(A,S) + conv(B, R), Ts,'variable','z^-1');
step(input)
figure(4)
U=tf(conv(A,R), P_end, Ts, 'variable', 'z^-1');
bodemag(U)

S3 = feedback(1,K*G3);
M_m3 =inv(norm(S3,inf)) % 0.1566
%% 2.4
load('TorMod.mat');
A = G3.f;
B = G3.b;
p1 = -1.7689;
p2 = 0.7939;
P = [1 p1 p2];

q = tf('q', 0.04);

Hs = [1 -1];
Hr = [1 1]; 
a = 0.3;
peaks = zeros(20, 3);
current_best = 100;
i = 1
for nq = 9
    for a=0.4
        B = G3.b;
        poles_aux = [a,a,a,a,a,a,a,a,a,a];
        coefs = poly(poles_aux);
        P_new = conv(P, coefs);

        [R, S] = poleplace(B, A, Hr, Hs, P_new) ;

        P_end = conv(A,S) + conv(B,R)

        T = sum(R);


        Ts = 0.04;
        CL = tf(conv(T,B), P_end, Ts,'variable','z^-1'); 
        figure(1)
        step(CL)
        U=tf(conv(A,R), P_end, Ts, 'variable', 'z^-1');
        figure(3)
        bodemag(U)% Does not fullfill the criterium
        figure(5)
        B = B(2:end);
        M_m = 0.4;
        U_max = 56.2; %35 dB  = 56.2
        q_delay = [0 1];

        R = [R, zeros(1, nq)];
        S = [S, zeros(1, nq)];
        new_R = @(Q) R + conv(A, conv(Hr, conv(Hs, Q)));
        new_S = @(Q) (S - conv(q_delay, conv(B, conv(Hs, conv(Hr, Q)))));

        c = @(Q) [norm(M_m*(S - conv(q_delay, conv(B, conv(Hs, conv(Hr, Q))))), Inf) - 1;
             norm(tf(conv(A, new_R(Q)), P_end, Ts, 'variable', 'z^-1'), Inf) - U_max];
        ceq = @(Q) [];
        zero = @(Q) 0;
        Mod_marg = @(Q) norm(tf(S + conv(q_delay, conv(B, conv(Hs, conv(Hr, Q)))), 1 ,Ts, 'variable', 'z^-1'), Inf)^(-1);
        Nonlincon = @(Q)deal(c(Q), ceq(Q));
        Q_opt = fmincon(Mod_marg, zeros(1, nq), [], [], [],[], [-Inf,-Inf], [Inf,Inf], Nonlincon);

        final_R = new_R(Q_opt);
        final_S = new_S(Q_opt);

        T = sum(final_R);


        Ts = 0.04;

         CL = tf(conv(T,G3.b), P_end, Ts,'variable','z^-1'); 
         U=tf(conv(A,final_R), P_end, Ts, 'variable', 'z^-1');
         bode(U)
         input= tf(conv(T,A), conv(A,final_S) + conv(conv(q_delay, B), final_R), Ts,'variable','z^-1');
         step(input)

    %     figure(11)

        aaa = stepinfo(CL);
       
       peaks(i,:) = [nq, a, aaa.RiseTime];
        i = i + 1;
    end
end
input = tf(conv(T,A), conv(A,S) + conv(conv(q_delay, B), R), Ts,'variable','z^-1');
current_best
peaks