function f = initialize_variables(N, V, min_range, max_range)
% ≥ı ºªØ
% initializing the chromosomes = [decision variables]
% 
% N - Population size
% V - Number of decision variables
% U - the numebr of uncertainty variables
% min_range & max_range - the maximum and minimum range for each decision variable, should be row vector.

%% Initialize each chromosome
f = rand(N,V).*(max_range-min_range)+min_range;
