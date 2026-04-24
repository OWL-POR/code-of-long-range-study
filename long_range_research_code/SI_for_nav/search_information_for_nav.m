function SI = search_information_for_nav(W, path_all, has_memory)
% SEARCH_INFORMATION                    Search information (Modified based on search_information.m (in Brain connectivity toolbox))
%
%   SI = search_information(W, L,has_memory)
%
%   Computes the amount of information (measured in bits) that a random
%   walker needs to follow the shortest path between a given pair of nodes.
%
%   Inputs:
%
%       W
%           Weighted/unweighted directed/undirected
%           connection weight matrix.

    %   path_all 
%
%
%      	has_memory,
%           This flag defines whether or not the random walker "remembers"
%           its previous step, which has the effect of reducing the amount
%           of information needed to find the next state. If this flag is
%           not set, the walker has no memory by default.
%
%
%   Outputs:
%
%       SI,
%           pair-wise search information of navigation strategy (matrix). Note that SI(i,j) may be
%           different from SI(j,i), hense, SI is not a symmetric matrix
%           even when adj is symmetric.


if ~exist('has_memory','var') 
    has_memory = false;
end

N = size(W,1); 

T = diag(sum(W,2))\W; 

SI = zeros(N,N); 
SI(eye(N)>0) = nan; 

for i = 1:N

    for j = 1:N

        if i ~= j

            path = path_all{i, j}; 
            
            lp = length(path); 
    
            if ~isempty(path) % Only calculate search information of the paths with successful navigation.
                pr_step_ff = nan(1,lp-1);
                if has_memory
                    pr_step_ff(1) = T(path(1),path(2));
                    for z=2:lp-1
                        pr_step_ff(z) = T(path(z),path(z+1))/(1 - T(path(z-1),path(z)));
                    end
                else
                    for z=1:length(path)-1
                        pr_step_ff(z) = T(path(z),path(z+1));
                    end
                end
                prob_sp_ff = prod(pr_step_ff);
                SI(i,j) = -log2(prob_sp_ff);
            else
                SI(i,j) = inf;
            end

        end

    end
end
