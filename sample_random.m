clear; clc;

%%
%------------------------------------%
% 1. define model & hyper parameters %
%------------------------------------%

% random adjacency matrix M
% N: number of nodes
% p: link probability (0 < p < 1)
% note: p should be small enough, 
%       otherwise the connectivity will be too large 
%       and thus calculation will be too heavy or impossible

N = 10; 
p = 0.05; 
M = rand_adj_matrix( N, p );

% interaction matrix
H = randn(5);
H = H + H';

% temperature
T = 0.1;


%%
%-----------------------------------------------------------------%
% 2. calculate the free energy for a given network & model H at T %
% note: currently, maximum connectivity should not be over 5      %
%-----------------------------------------------------------------%

F = free_energy( T, H, M )


