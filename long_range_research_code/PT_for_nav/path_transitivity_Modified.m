function [T] = path_transitivity_Modified(W)
% PATH_TRANSITIVITY             Transitivity based on navigation paths (Modified based on path_transitivity.m (in Brain connectivity toolbox))
%
%   T = path_transitivity(path)
%
%   This function computes the density of local detours (triangles) that
%   are available along the navigation paths between all pairs of nodes.
%
%   Inputs:
%
%       path,
%           path node sequence of navigation
%
%   Output:
%
%       T,
%           matrix of pairwise path transitivity.
%

n=length(W);
m=zeros(n,n);
T=zeros(n,n);

for i=1:n-1
    for j=i+1:n
        x=0;
        y=0;
        z=0;
        for k=1:n
            if W(i,k)~=0 && W(j,k)~=0 && k~=i && k~=j
                x=x+W(i,k)+W(j,k);
            end
            if k~=j
                y=y+W(i,k);
            end
            if k~=i
                z=z+W(j,k);
            end
        end
        m(i,j)=x/(y+z);
    end
end
m=m+m';


L_matrix = W;

n = size(W, 1);

L_matrix(L_matrix ~= 0) = 1./L_matrix(L_matrix~=0); 

ED = table2array(readtable("stand_ED_BN246.csv"));

[~, ~, ~, ~, paths] = navigation_wu_my(L_matrix, ED);

T=inf(n,n);


% --- path transitivity ---%%

for i=1:n
    for j=1:n
        x=0;
        path = paths{i,j};
        if ~isempty(path)

            K=length(path);
            
            for t=1:K-1
                for l=t+1:K
                    x=x+m(path(t),path(l));
                end
            end
            T(i,j)=2*x/(K*(K-1)); 

        else

            continue

        end
    end
end
