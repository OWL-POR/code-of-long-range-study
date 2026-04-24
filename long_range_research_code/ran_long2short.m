clc;
clear;

% create reduced surrogates

matrix_con = table2array(readtable("thr_30\binary_all.csv")); % group-level SC mask (density 30%)

long_range = table2array(readtable("thr_30\long_25.csv")); % long range connections

origin_mean_length = table2array(readtable("SL\streamline_length_mean.csv")); % Non-thresholded average streamline-length matrix (nearly 100% density)

length_cal = origin_mean_length;

length_cal(long_range == 0) = NaN;

min_long_connection = min(length_cal(~isnan(length_cal)));

cal_index_matrix = matrix_con;

cal_index_matrix(origin_mean_length == 0) = NaN;

long_index = find(triu(origin_mean_length >= min_long_connection,1));

cal_index_matrix(long_index) = NaN;

short_index = find(triu(cal_index_matrix == 0,1));

strength_matrix = table2array(readtable("consense\strength_mean.csv")); % Non-thresholded average connection strength matrix (nearly 100% density)

strength_matrix_used = strength_matrix;

strength_matrix_triu  = triu(strength_matrix_used,1);

long_range_strength = strength_matrix_triu(long_range == 1);

long_range_strength = long_range_strength(long_range_strength ~= 0);

strength_matrix_used(matrix_con == 0) = 0;

strength_matrix_used(long_range == 1) = 0; 

for net_i = 1:ran_num

    random_index = short_index(randperm(length(short_index), long_all));

    wei_random_long2short = zeros(246,246);

    wei_random_long2short(random_index) = long_range_strength;

    wei_random_long2short = wei_random_long2short + wei_random_long2short';

    matrix_cal_strength  = strength_matrix_used;

    matrix_cal_strength = matrix_cal_strength + wei_random_long2short;

    writematrix(matrix_cal_strength, strcat("thr_30\random_long2short_weighted\network\random", num2str(net_i),".csv")); % matrix contains shorten long-range connections

end

