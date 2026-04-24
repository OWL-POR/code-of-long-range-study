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
    
            node_cos_sim(node_neigh_pair_i) = dot(vector_neigh_a, vector_neigh_b)/(norm(vector_neigh_a)*norm(vector_neigh_b));
    
        end

        node_each_cos_sim = [node_each_cos_sim, mean(node_cos_sim)];

    end

end

real_cos_sim = mean(node_each_cos_sim); 

strength_matrix(long_range == 1) = 0; 

net_num = 1000;

cos_sim_vector_long_ran_net = zeros(net_num, 1);

for net_long_ran_i = 1:net_num
    
    random_net_long_ran = table2array(readtable(strcat("spatial_surrogate\\random", num2str(net_long_ran_i), ".csv"), 'ReadVariableNames',false)); % spatial surrogates (only long-range connections)

    strength_matrix_longarr = strength_matrix + random_net_long_ran;

    ran_long_node_each_cos_sim = [];
     
    for node_i = 1:size(matrix_con,1)
    
        ran_long_edge_node_i = random_net_long_ran(:, node_i);
    
        ran_long_node_i_neigh = find(ran_long_edge_node_i);
    
        if length(ran_long_node_i_neigh) >= 2
        
            pair_index = nchoosek(ran_long_node_i_neigh, 2);
        
            node_cos_sim_ran_long = zeros(nchoosek(length(ran_long_node_i_neigh), 2), 1);
        
            for node_neigh_pair_i = 1:nchoosek(length(ran_long_node_i_neigh), 2)
        
                vector_neigh_a = strength_matrix_longarr(:, pair_index(node_neigh_pair_i, 1));
        
                vector_neigh_b = strength_matrix_longarr(:, pair_index(node_neigh_pair_i, 2));
        
                node_cos_sim_ran_long(node_neigh_pair_i) = dot(vector_neigh_a, vector_neigh_b)/(norm(vector_neigh_a)*norm(vector_neigh_b));
        
            end
    
            ran_long_node_each_cos_sim = [ran_long_node_each_cos_sim, mean(node_cos_sim_ran_long)];
    
        end
    
    end

    cos_sim_vector_long_ran_net(net_long_ran_i, 1) = mean(ran_long_node_each_cos_sim);

end

save("redundant_arr.mat", 'real_cos_sim', 'cos_sim_vector_long_ran_net')