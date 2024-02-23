% Computes statistics from an n x m matrix, where n is number of
% observations and m is the number of metrics measured. Generates an m x 7
% matrix with the following statistics: Min, LQ, Median, UQ, Max, Mean,
% Stdev

function stats_summary = compute_stats_sum(input_matrix)
    stats_summary = [];
    stats_summary = [stats_summary, min(input_matrix)'];
    stats_summary = [stats_summary, quantile(input_matrix, 0.25)'];
    stats_summary = [stats_summary, quantile(input_matrix, 0.5)'];
    stats_summary = [stats_summary, quantile(input_matrix, 0.75)'];
    stats_summary = [stats_summary, max(input_matrix)'];
    stats_summary = [stats_summary, mean(input_matrix)'];
    stats_summary = [stats_summary, std(input_matrix)'];
end