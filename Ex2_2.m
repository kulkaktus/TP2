load('TorMod.mat');
A = G3.f;
B = G3.b;

p1 = -0.5;
p2 = 0.5;
P = [1 p1];

Hs = [1 -1];
Hr = [1];

[R, S] = poleplace(B, A, Hr, Hs, P)