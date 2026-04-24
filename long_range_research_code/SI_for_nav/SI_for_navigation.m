function [mean_SI_nav] = SI_for_navigation(W)

    % input: SC matrix

    ED = table2array(readtable("stand_ED_BN246.csv")); % Euclidean distance of BNA246 atlas (in MNI space)
    
    L_matrix = W;

    n = size(W, 1);
    
    L_matrix(L_matrix ~= 0) = 1./L_matrix(L_matrix~=0);
    
    [~, ~, ~, ~, paths_real] = navigation_wu_Modified(L_matrix, ED); 
    
    SI_nav = search_information_for_nav(W, paths_real, []); 

    SI_nav(1:n+1:end) = inf;

    mean_SI_nav = mean(SI_nav(~isinf(SI_nav)));

end