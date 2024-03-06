function M = rand_adj_matrix( N, p )
% Copyright 2020, Hyun-Yong Lee (hyunyong@korea.ac.kr)
% Function: Generating Random Adjacency Matrix
% (Input) N: Size of matrix
% (Input) p: Link probability
% (Output) M: Randon adjacency matrix (symmetric)
    
    M = zeros(N);
    
    for i=1:(N-1)
        
        for j=(i+1):N
            x = rand;
            if( x < p )
                M(i,j) = 1;
            end
        end
        
        if( sum(M(i,:))==0 )
            M(i,randi([i+1,N])) = 1;
        end
        
    end
    
    M = M+M';

end

