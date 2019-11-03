%input = W_AW, W_fuel, W_str, S, n_factor
% altitude = 36500 ft = 11125.2 m
% rho_cruise = 0.356804 kg/m3  (ISA)
% V = 432 knots = 222.24 m/s
function [CL_design] = DPC_function(W_AW, W_fuel, W_str, S, n_factor)
    W_TO_max = W_AW + W_fuel + W_str;
    W_design = sqrt(W_TO_max*(W_TO_max - W_fuel));
    CLdes = (2*W_design)/(S* 0.356804 *222.24^2);
    CL_design = CLdes * n_factor;
end