clc;
clear;

%% calculate diversity

matrix_con = table2array(readtable("binary_all.csv", 'ReadVariableNames',false)); % mask 30% density

strength_matrix = table2array(readtable("strength_mean.csv", 'ReadVariableNames',false)); % Non-thresholded average connection strength matrix (nearly 100% density)

strength_matrix(matrix_con == 0) = 0;

long_range = table2array(readtable("long_25.csv", 'ReadVariableNames',false)); % long-range mask

E = long_range;

[edges_u, edges_v] = find(triu(E,1)); 

long_edges = [edges_u, edges_v];

cos_sim_vector = zeros(size(long_edges,1),1);

% calculate network diversity for real network

for node_pair_i = 1:size(long_edges,1)

    vector_node_a = strength_matrix(:, long_edges(node_pair_i, 1));

    vector_node_b = strength_matrix(:, long_edges(node_pair_i, 2));

    cos_sim_vector(node_pair_i, 1) = dot(vector_node_a, vector_node_b)/(norm(vector_node_a)*norm(vector_node_b));

end

mean_real_cos_sim = mean(cos_sim_vector);

%% calculate diversity for shorten connections

strength_matrix_short = strength_matrix;

strength_matrix_short(long_range == 1) = 0;

net_num = 1000;

cos_sim_vector_shorten = zeros(net_num, 1);

for net_short_i = 1:net_num
    
    random_net_short = table2array(readtable(strcat("random_long2short\\network_25\\random", num2str(net_short_i), ".csv"), 'ReadVariableNames',false));

    random_net_long2short = random_net_short;
    
    random_net_long2short(strength_matrix_short > 0) = 0;

    [edges_short_u, edges_short_v] = find(triu(random_net_long2short,1));
    
    short_edges = [edges_short_u, edges_short_v];
    
    cos_sim_vector_short = zeros(size(short_edges,1),1);
    
    for node_pair_i = 1:size(short_edges,1)

        vector_node_a_short = random_net_short(:, short_edges(node_pair_i, 1));
    
        vector_node_b_short = random_net_short(:, short_edges(node_pair_i, 2));
    
        cos_sim_vector_short(node_pair_i, 1) = dot(vector_node_a_short, vector_node_b_short)/(norm(vector_node_a_short)*norm(vector_node_b_short));
    
    end

    cos_sim_vector_shorten(net_short_i, 1) = mean(cos_sim_vector_short);

end

%% calculate diversity for permu long connections

cos_sim_vector_long_permu = zeros(net_num, 1);

for net_long_ran_i = 1:net_num
    
    random_net_long_ran = table2array(readtable(strcat("random_arrange_weight\\network_25\\random", num2str(net_long_ran_i), ".csv"), 'ReadVariableNames',false));

    [edges_long_ran_u, edges_long_ran_v] = find(triu(random_net_long_ran,1)); 
    
    long_ran_edges = [edges_long_ran_u, edges_long_ran_v];
    
    cos_sim_vector_long_ran = zeros(size(long_ran_edges,1),1);

    ran_long_matrix = random_net_long_ran + strength_matrix_short;
    
    for node_pair_i = 1:size(long_ran_edges,1)

        vector_node_a_long_ran = ran_long_matrix(:, long_ran_edges(node_pair_i, 1));
    
        vector_node_b_long_ran = ran_long_matrix(:, long_ran_edges(node_pair_i, 2));
    
        cos_sim_vector_long_ran(node_pair_i, 1) = dot(vector_node_a_long_ran, vector_node_b_long_ran)/(norm(vector_node_a_long_ran)*norm(vector_node_b_long_ran));
    
    end

    cos_sim_vector_long_permu(net_long_ran_i, 1) = mean(cos_sim_vector_long_ran);

end

save("density.mat", "mean_real_cos_sim", "cos_sim_vector_shorten", "cos_sim_vector_long_permu");