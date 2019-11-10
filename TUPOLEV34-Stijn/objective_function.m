%input is fuel weight and wing structural weight
function [W_TO_max] = objective_function(W_AW, W_fuel, W_str)
    W_TO_max = W_AW + W_fuel + W_str;
end