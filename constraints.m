function [c, ceq] = constraints(x)

global data;
y_fte = data.y_85;
c_fte = data.c_85;
Xtur = data.Xtur;
Xtut = data.Xtut;
Xtuk = data.Xtuk;
Xtu85 = data.Xtu85;
Xtlr = data.Xtlr;
Xtlt = data.Xtlt;
Xtlk = data.Xtlk;
Xtl85 = data.Xtl85;
W_TO_max = data.W_TO_max;
W_fuel = data.W_fuel;
S = data.S;
y_section = data.y_section;
c_section = data.c_section;

f_tank = 0.93;
rho_fuel = 0.81715*10^3;
V_fuel = W_fuel/9.80665/rho_fuel;
WS_original = 558.72;


[V_tank] = fuel_volume(Xtur, Xtuk, Xtu85, Xtlr, Xtlk, Xtl85, y_fte, c_fte, y_section, c_section);
[WS_new] = Wing_Loading_function(W_TO_max, S);

c = [(V_fuel - V_tank), (WS_new - WS_original)];
ceq = [];
