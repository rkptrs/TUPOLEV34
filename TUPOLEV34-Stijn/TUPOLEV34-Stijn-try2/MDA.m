function [W_TO_max_out] = MDA(x)

global data;

x = x.*(data.ub-data.lb)+data.lb;
x(1:7)  %check the input of the loop

W_fuel = data.W_fuel_mda;
W_str = data.W_str_mda;
[x_section, y_section, z_section, c_section, twist_section, S, y_85, c_85] = geometry_function(x(1:7), data.y_kink, data.dihedral);

W_TO_max = data.W_TO_max_mda;
W_TO_list = [W_TO_max, 0];
diff_W_TO_max = abs(W_TO_list(1) - W_TO_list(2));

[Xtur,Xtlr,Xtut,Xtlt,Xtuk,Xtlk,Xtu85,Xtl85] = CSTtoGeom(x(8:13), x(14:19), x(20:25), x(26:31), data.y_kink/x(3),c_section,c_85);

taper = x(2)/x(1);
mac = x(1)*2/3*(1 + taper + taper^2)/(1 + taper);
Re = data.rho*data.V*mac/(1.422*10^-5);

no_convergence_value = 1;
iteration_mda = 0;
while diff_W_TO_max > 9.80665*0.1
W_TO_list(1) = W_TO_max;
[CL_design_max] = DPC_function(data.W_AW, W_fuel, W_str, S, data.n_max);
[cl_distribution, cm_distribution, Y_distribution, chord_distribution, CLwing, CDwing] = Q3D_function(x_section, y_section, z_section, c_section, twist_section, x(8:19), x(20:31), 0, 150, CL_design_max, data.V, data.rho, data.alt, Re, data.M);
EMWET_function(x(1:7), cl_distribution, cm_distribution, Y_distribution, chord_distribution, W_TO_max, W_fuel, data.n_max, S, x_section, y_section, z_section, c_section, data.rho, data.V);
EMWET TUPOLEVaircraft
filetext = fileread('TUPOLEVaircraft.weight');
W_str_kg = str2double(regexp(filetext, '(?<=Wing total weight[^0-9]*)[0-9]*\.?[0-9]+', 'match'));
W_str = W_str_kg * 9.80665;    
[CL_design_cruise] = DPC_function(data.W_AW, W_fuel, W_str, S, 1);
% x
% W_TO_max, W_fuel, W_str
[cl_distribution, cm_distribution, Y_distribution, chord_distribution, CLwing, CDwing] = Q3D_function(x_section, y_section, z_section, c_section, twist_section, x(8:19), x(20:31), 1, 100, CL_design_cruise, data.V, data.rho, data.alt, Re, data.M);
if isnan(CDwing)
    CDwing = 0.06;
    no_convergence_value = 0;  
end
[W_fuel] = Brequet_function(CLwing, CDwing, data.W_AW, W_fuel, W_str, data.D_AW, S, data.rho, data.V);
[W_TO_max] = objective_function(data.W_AW, W_fuel, W_str);
% [W_TO_max/9.80665/1000; W_str/9.80665/1000]
W_TO_list(2) = W_TO_max;
diff_W_TO_max = abs(W_TO_list(1) - W_TO_list(2));
% W_TO_list(1) - W_TO_list(2)
iteration_mda = iteration_mda + 1
if iteration_mda == 10
    diff_W_TO_max = 0;  %to end mda loop if it wont converge further
end
end

W_TO_max_out = W_TO_max / data.W_TO_max_ref;


%save new guess for next loop
W_str_mda = W_str;
W_fuel_mda = W_fuel;
W_TO_max_mda = W_TO_max;
% abs(sum(data.xold) - sum(x))*100000    %check to see if x changes
xold = x;
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
data.W_str_mda = W_str_mda;
data.W_fuel_mda = W_fuel_mda;
data.W_TO_max_mda = W_TO_max_mda;
data.xold = xold;
data.no_convergence_value = no_convergence_value;

