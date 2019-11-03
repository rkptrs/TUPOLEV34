%input = W_to_max and S
%output is in kg/m2!!!!!
% original Wingloading =  558.72 kg/m2
function [Wing_Loading_max] = Wing_Loading_function(W_TO_max, S)
    Wing_Loading_max = W_TO_max/S/9.80665;
end