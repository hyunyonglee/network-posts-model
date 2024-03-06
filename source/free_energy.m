function [F, cost] = free_energy( temp, H, M )
% Copyright 2020, Hyun-Yong Lee (hyunyong@korea.ac.kr)
% Function: Calculating free energy
% (Input) temp: temperature
% (Input) H: interaction matrix
% (Input) M: (symmetric) adjacency matrix
% (Ouput) F: free energy
% (Ouput) cost: cost for the tensor network contraction

    % Calculate all necessary variables for the tensor network 
    % for a given adjacency matrix M
    [sequence, cost, ~, ~, ~, bond_leg, all_leg] = adj_matrix_to_contraction_order(M);

    % Boltzman factor matrix
    expH = exp(-H/temp);
    
    % normalization factor
    z0 = max(expH(:));
    
    % (normalized) boltzman factor matrix
    expH = expH/z0;
    
    % Trivial link tensors
    % Here, I assume the largest connectivity in the network is 5
    % If you have larger connectivity, just add T6, T7 and so on.
    % But, the complexity will grows exponentially with the connectivity
    q = size(H,1);
    T1 = ones(q,1);
    T2 = eye(q);
    T3 = zeros(q,q,q);
    T4 = zeros(q,q,q,q);
    T5 = zeros(q,q,q,q,q); % Maximum connectivity
    
    for i=1:q
        T3(i,i,i) = 1;
        T4(i,i,i,i) = 1;
        T5(i,i,i,i,i) = 1;
    end
    
    % prepare tensor contraction
    Ts = cell(1,length(all_leg));
    N1 = length(all_leg);
    N2 = length(bond_leg);
    
    for i=1:(N1-N2)
        
        if( length(all_leg{i})>5 )
            error('Error: Maximum connectivity of the network is over 5')
        end

        if( length(all_leg{i}) == 1); T = T1;
        elseif( length(all_leg{i}) == 2); T = T2;
        elseif( length(all_leg{i}) == 3); T = T3;
        elseif( length(all_leg{i}) == 4); T = T4;
        elseif( length(all_leg{i}) == 5); T = T5;
        end
        Ts{i} = T;
        
    end
    
    for i=1:N2 
        Ts{ N1-N2+i } = expH; 
    end
    
    % (normalized) partition function
    Z = ncon(Ts, all_leg, sequence);
    
    % free energy
    F = - temp * ( N2 * log(z0) + log(Z) );


end

