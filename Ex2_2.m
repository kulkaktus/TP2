%% 2.2
load('TorMod.mat');
A = G3.f;
B = G3.b;
p1 = -1.7689;
p2 = 0.7939;
P = [1 p1 p2];


Hs = [1 -1];
Hr = [1]; 
[R, S] = poleplace(B, A, Hr, Hs, P) 

P_end = conv(A,S) + conv(B,R)
T = sum(R)


%% Ex2_3
Ts = 0.04;
CL = tf(conv(T,B), P_end, Ts,'variable','z^-1')
figure(1)
step(CL)
stepinfo(CL)

K = tf(R,S, Ts,'variable','z^-1')
OL = G3*K;
figure(2)
nyquist(OL) % 3 encriclments around -1 corresponds to 3 unstable poles 

U=tf(conv(A,R), P, Ts, 'variable', 'z^-1');
figure(3)
step(U)% Does not fullfill the criterium
figure(4)
bode(U)

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

[R, S] = poleplace(B, A, Hr, Hs, P) 

P_end = conv(A,S) + conv(B,R)

T = sum(R)

Ts = 0.04;
CL = tf(conv(T,B), P_end, Ts,'variable','z^-1')
figure(1)
step(CL)
stepinfo(CL)
