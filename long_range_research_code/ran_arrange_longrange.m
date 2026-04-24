clc;
clear;

% create spatial surrogates

bin_num = 1; % Stepwise increase the number of bins, until until the average streamline length of long-range connections in the real group-level matrix fell between the 45th and 55th percentiles of 1000 surrogate networks

origin_mean_length = table2array(readtable("SL\streamline_length_mean.csv")); % Non-thresholded average streamline-length matrix (nearly 100% density)

long_range_pick = table2array(readtable(strcat("long_25.csv"))); % longrange connections

origin_mean_weight = table2array(readtable("consense\strength_mean.csv")); % Non-thresholded average connection strength matrix (nearly 100% density)

origin_mean_weight(long_range_pick == 0) = 0;

length_cal = origin_mean_length;

length_cal(long_range_pick == 0) = NaN; 

min_long_connection = min(length_cal(~isnan(length_cal)));

max_long_connection = max(length_cal(~isnan(length_cal)));

random_num = 1000;

distbins = linspace(min_long_connection,max_long_connection,bin_num + 1); 

mean_length = zeros(1000,1);

long_range_all = origin_mean_length;

real_long_range = origin_mean_weight;

for net_i = 1:random_num
    matrix_all = zeros(246,246);
    for bin = 1:bin_num

        if bin ~= bin_num

            bin_index = find(triu(long_range_all >= distbins(bin) & long_range_all < distbins(bin + 1),1));

            long_range_strenth = real_long_range(bin_index);

            long_range_strenth = long_range_strenth(long_range_strenth ~= 0);

            bin_random = zeros(246,246);

            random_index = bin_index(randperm(length(bin_index), length(long_range_strenth)));

            bin_random(random_index) = long_range_strenth;

            matrix_all = matrix_all + bin_random;
        else
            bin_index = find(triu(long_range_all >= distbins(bin) & long_range_all <= distbins(bin + 1),1));

            long_range_strenth = real_long_range(bin_index);

            long_range_strenth = long_range_strenth(long_range_strenth ~= 0);

            bin_random = zeros(246,246);

            random_index = bin_index(randperm(length(bin_index), length(long_range_strenth)));

            bin_random(random_index) = long_range_strenth;

            matrix_all = matrix_all + bin_random;
        end

    end
    matrix_all = matrix_all + matrix_all';

    writematrix(matrix_all, strcat("random_arrange_weight\network_permu_longrange\random", num2str(net_i), ".csv")); % This matrix contains permuted long-range connections and should be added to the matrix with only short-range connections to construct the surrogate model.
   
    mean_length(net_i) = mean(long_range_all(matrix_all > 0));

end

real_mean_length = mean(origin_mean_length(long_range_pick == 1));

disp(real_mean_length < prctile(mean_length, 55) & real_mean_length > prctile(mean_length, 45))
