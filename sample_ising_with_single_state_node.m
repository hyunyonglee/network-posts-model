clear; clc;

%%
%-------------------%
% 1. define network %
%-------------------%

% adjacency matrix
N = 20;
M = diag( ones(N-1,1),1);
M(1,N) = 1;
M = M + M';

Vs = zeros(N,1);
Vs(1:2:N) = 1;

Ms = M;
Ms(1:N,N+1) = Vs(:);
Ms(N+1,1:N) = Vs(:);


figure;
plot(graph(M))
title('without the special node')

figure;
plot(graph(Ms))
title('with the special node')

%% 
%-----------------%
% 2. define model %
%-----------------%

% interaction matrix
h = 0.2;
H = eye(2);
Hs = [0; h];

% temperature
T = 1;

%%
%-----------------------------------------------------------------%
% 3. calculate the free energy for a given network & model H at T %
% note: currently, maximum connectivity should not be over 5      %
%-----------------------------------------------------------------%

F = free_energy_with_single_state_node( T, H, M, Hs, Vs )


%%
% (Benchmark) N-Ising chain with pbc using transfer matrix

R = [1 0; 0 exp(-h/T)];
B = exp(-H/T);
F_exact = -T * log( sum( diag( (B*R*B)^(N/2) )  ) )
