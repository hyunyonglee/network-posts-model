function [sequence, cost, subnet_sequence, subnet_tensors, site_leg, bond_leg, all_leg] = adj_matrix_to_contraction_order(M)
% Copyright 2020, Hyun-Yong Lee (hyunyong@korea.ac.kr)
% Function: Calculating all tensor contraction (ncon) requirements for a given adjacency matrix
% (Input) M: Adjacency Matrix
% (Outputs) ...: All necessary variables for ncon function

    
    [N,~] = size(M);
    
    site_leg = cell(1,N);
    bond_leg = cell(1,sum(M(:))/2);
    I = 1;
    
    for i=1:N
        
        leg = [];
        for j=1:N
            
            if( M(i,j) > 1.0e-6 )
                
                if( i < j )
                    C = 10000000 * i + 1000 * j + 1;
                    bond_leg{I} = [C, C+1];
                    I = I+1;
                else
                    C = 10000000 * j + 1000 * i + 2;
                end
                leg = [leg, C];
                
            end
            
        end
        site_leg{i} = leg;
        
    end
    
    
    all_leg = [site_leg(:)', bond_leg(:)'];
    [sequence, cost, subnet_sequence, subnet_tensors] = netcon(all_leg,0);

end

