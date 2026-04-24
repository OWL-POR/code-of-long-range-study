clc;
clear;

%% calculate connection degeneracy of permu long-range connections

%% real network

matrix_con = table2array(readtable("binary_all.csv", 'ReadVariableNames',false)); % mask 30% density

strength_matrix = table2array(readtable("strength_mean.csv", 'ReadVariableNames',false)); % Non-thresholded average connection strength matrix (nearly 100% density)

strength_matrix(matrix_con == 0) = 0;

long_range = table2array(readtable("long_25.csv", 'ReadVariableNames', false)); % long-range connection matrix

node_each_cos_sim = [];
 
for node_i = 1:size(matrix_con,1)

    edge_node_i = long_range(:, node_i);

    node_i_neigh = find(edge_node_i);

    if length(node_i_neigh) >= 2
    
        pair_index = nchoosek(node_i_neigh, 2);
    
        node_cos_sim = zeros(nchoosek(length(node_i_neigh), 2), 1);
    
        for node_neigh_pair_i = 1:nchoosek(length(node_i_neigh), 2)
    
            vector_neigh_a = strength_matrix(:, pair_index(node_neigh_pair_i, 1));
    
            vector_neigh_b = strength_matrix(:, pair_index(node_neigh_pair_i, 2));
    
            node_cos_sim(node_neigh_pair_i, 1) = dot(vector_neigh_a, vector_neigh_b)/(norm(vector_neigh_a)*norm(vector_neigh_b));
    
        end

        node_each_cos_sim = [node_each_cos_sim, mean(node_cos_sim)]; 

    end

end

real_cos_sim = mean(node_each_cos_sim); 

net_num = 1000;

strength_matrix_short = strength_matrix;

strength_matrix_short(long_range == 1) = 0; 

cos_sim_vector_long2short_net = zeros(net_num, 1);

for net_long2short_i = 1:net_num

    random_net_short = table2array(readtable(strcat("reduced_surrogate\\network_25\\random", num2str(net_long2short_i), ".csv"), 'ReadVariableNames',false)); % reduced surrogates

    random_net_long2short = random_net_short;
    
    random_net_long2short(strength_matrix_short > 0) = 0; % only shorten long-range connections

    long2short_node_each_cos_sim = [];
     
    for node_i = 1:size(matrix_con,1)
    
        long2short_edge_node_i = random_net_long2short(:, node_i);
    
        long2short_node_i_neigh = find(long2short_edge_node_i);
    
        if length(long2short_node_i_neigh) >= 2
        
            pair_index = nchoosek(long2short_node_i_neigh, 2);
        
            node_cos_sim_long2short = zeros(nchoosek(length(long2short_node_i_neigh), 2), 1);
        
            for node_neigh_pair_i = 1:nchoosek(length(long2short_node_i_neigh), 2)
        
                vector_neigh_a = random_net_short(:, pair_index(node_neigh_pair_i, 1)); 
        
                vector_neigh_b = random_net_short(:, pair_index(node_neigh_pair_i, 2));
        
                node_cos_sim_long2short(node_neigh_pair_i, 1) = dot(vector_neigh_a, vector_neigh_b)/(norm(vector_neigh_a)*norm(vector_neigh_b));
        
            end
    
            long2short_node_each_cos_sim = [long2short_node_each_cos_sim, mean(node_cos_sim_long2short)];
    
        end
    
    end

    cos_sim_vector_long2short_net(net_long2short_i, 1) = mean(long2short_node_each_cos_sim); 

end

save("redundant_long2short.mat", 'real_cos_sim', 'cos_sim_vector_long2short_net')