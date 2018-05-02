load('TorMod.mat');
A = G3.f;
B = G3.b;
d = [0 0 1];
p1 = 0.7939;
p2 = -1.7689;
P = [1 p1 p2];

Hs = [1 -1];
Hr = [1]; 

[R, S] = poleplace(B, A, Hr, Hs, P) 

P_end = conv(A,S) + conv(B,R)
T = sum(R)


%% Ex2_3
CL = tf(conv(T,B), P_end, Ts,'variable','z^-1')
step(CL)

K = deconv(R,S); %Comment on the result
U=tf(conv(A,R), P, Ts, 'variable', 'z^-1');