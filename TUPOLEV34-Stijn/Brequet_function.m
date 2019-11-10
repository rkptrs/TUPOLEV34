%input = CLwing, CDwing, W_AW, W_fuel_old, W_str
% R = 1550 nautical mile = 2870600 m
% V = 432 knots = 222.24 m/s
% CT = 1.8639*10^-4
function [W_fuel] = Brequet_function(CLwing, CDwing, W_AW, W_fuel_old, W_str, D_AW, S, rho, V)
    C_D_AW = D_AW/(1/2*rho*S*V^2);
    CLCD = CLwing/(CDwing+C_D_AW);
    W_fraction_cruise = exp(2870600/((222.24/(1.8639*10^-4))*(CLCD)));
    W_fuel = (1 - 0.938/W_fraction_cruise)*(W_AW + W_fuel_old + W_str); 
end