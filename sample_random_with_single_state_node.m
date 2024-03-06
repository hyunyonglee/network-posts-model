clear; clc;

%%
%------------------------------------%
% 1. define model & hyper parameters %
%------------------------------------%

% random adjacency matrix M
% N: dimensoin
% p: link probability (0 < p < 1)
% note: p should be small enough, 
%       otherwise the connectivity will be too large 
%       and thus calculation will be too heavy or impossible

N = 20; 
p = 0.05; 
M = rand_adj_matrix( N, p );
Vs = randi([0,1],[N,1]);

Ms = M;
Ms(1:N,N+1) = Vs(:);
Ms(N+1,1:N) = Vs(:);

figure;
plot(graph(M))
title('without the special node')

figure;
plot(graph(Ms))
title('with the special node')

% interaction matrix
H = randn(5);
H = H + H';

Hs = randn(5,1);

% temperature
T = 0.1;


%%
%-----------------------------------------------------------------%
% 2. calculate the free energy for a given network & model H at T %
% note: currently, maximum connectivity should not be over 5      %
%-----------------------------------------------------------------%
clc;
F = free_energy_with_single_state_node( T, H, M, Hs, Vs )


