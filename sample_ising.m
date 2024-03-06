clear; clc;

%%
%------------------------------------%
% 1. define model & hyper parameters %
%------------------------------------%

% adjacency matrix
N = 20;
M = diag( ones(N-1,1),1);
M(1,N) = 1;
M = M + M';

% interaction matrix
H = eye(2);

% temperature
T = 1;

%%
%-----------------------------------------------------------------%
% 2. calculate the free energy for a given network & model H at T %
% note: currently, maximum connectivity should not be over 5      %
%-----------------------------------------------------------------%

F = free_energy( T, H, M )


%%
% (Benchmark) N-Ising chain with pbc using transfer matrix

M = exp(-H/T);
F_exact = -T * log( sum( diag(M^N) ) )
