function [c, ceq] = constraints(x)

global data;
% y_fte = data.y_85;
% c_fte = data.c_85;
% Xtur = data.Xtur;
% Xtut = data.Xtut;
% Xtuk = data.Xtuk;
% Xtu85 = data.Xtu85;
% Xtlr = data.Xtlr;
% Xtlt = data.Xtlt;
% Xtlk = data.Xtlk;
% Xtl85 = data.Xtl85;
% W_TO_max = data.W_TO_max;
% W_fuel = data.W_fuel;
% S = data.S;
% y_section = data.y_section;
% c_section = data.c_section;

x = x.*(data.ub-data.lb)+data.lb;

W_fuel = data.W_fuel_mda;
W_str = data.W_str_mda;
[x_section, y_section, z_section, c_section, twist_section, S, y_fte, c_fte] = geometry_function(x(1:7), data.y_kink, data.dihedral);
W_TO_max = data.W_TO_max_mda;
[Xtur,Xtlr,Xtut,Xtlt,Xtuk,Xtlk,Xtu85,Xtl85] = CSTtoGeom(x(8:13), x(14:19), x(20:25), x(26:31), data.y_kink/x(3),c_section,c_fte);

f_tank = 0.93;
rho_fuel = 0.81715*10^3;
V_fuel = W_fuel/9.80665/rho_fuel;
WS_original = 558.72;

[V_tank] = fuel_volume(Xtur, Xtuk, Xtu85, Xtlr, Xtlk, Xtl85, y_fte, c_fte, y_section, c_section);
[WS_new] = Wing_Loading_function(W_TO_max, S);

no_convergence_value = data.no_convergence_value;

c = [(V_fuel - V_tank*f_tank)/abs(data.c_volume_ref), (WS_new - WS_original)/abs(data.c_wingloading_ref)];
%ceq = [(no_convergence_value - 1)];
ceq = [];