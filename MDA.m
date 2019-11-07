function [W_TO_max] = MDA(x)

s = load('reffile.mat');
W_fuel = s.W_fuel_ref;
W_str = s.W_str_ref;
[x_section, y_section, z_section, c_section, twist_section, S, y_85, c_85] = geometry_function(x(1:7), s.y_kink, s.dihedral);

W_TO_list = [s.W_TO_max, 0];
diff_W_TO_max = abs(W_TO_list(1) - W_TO_list(2));

while diff_W_TO_max > 0.05
W_TO_list(1) = W_TO_max;
[CL_design_max] = DPC_function(s.W_AW, W_fuel, W_str, S, s.n_max);
[cl_distribution, cm_distribution, Y_distribution, chord_distribution, CLwing, CDwing] = Q3D_function(x_section, y_section, z_section, c_section, twist_section, x(8:19), x(20:31), 0, 150, CL_design_max, s.V, s.rho, s.alt, s.Re, s.M);
[W_str] = EMWET_function(x_1, cl_distribution, cm_distribution, Y_distribution, chord_distribution, W_TO_max, W_fuel, s.n_max, S, x_section, y_section, z_section, c_section, s.rho, s.V);
[CL_design_cruise] = DPC_function(s.W_AW, W_fuel, W_str, S, s.n_max);
[cl_distribution, cm_distribution, Y_distribution, chord_distribution, CLwing, CDwing] = Q3D_function(x_section, y_section, z_section, c_section, twist_section, x(8:19), x(20:31), 1, 150, CL_design_cruise, s.V, s.rho, s.alt, s.Re, s.M);
[W_fuel] = Brequet_function(CLwing, CDwing, W_AW, W_fuel_old, W_str, D_AW, S, rho, V);
[W_TO_max] = objective_function(W_AW, W_fuel, W_str);
W_TO_list(2) = W_TO_max;
diff_W_TO_max = abs(W_TO_list(1) - W_TO_list(2));
end

[Xtur,Xtlr,Xtut,Xtlt,Xtuk,Xtlk,Xtu85,Xtl85] = CSTtoGeom(x(8:13), x(14:19), x(20:25), x(26:31), s.y_kink/x(3));
global data;
data.y_85 = y_85;
data.c_85 = c_85;
data.Xtur = Xtur;
data.Xtut = Xtut;
data.Xtuk = Xtuk;
data.Xtu85 = Xtu85;
data.Xtlr = Xtlr;
data.Xtlt = Xtlt;
data.Xtlk = Xtlk;
data.Xtl85 = Xtl85;
data.W_TO_max = W_TO_max;
data.W_fuel = W_fuel;
data.S = S;
data.y_section = y_section;
data.c_section = c_section;
