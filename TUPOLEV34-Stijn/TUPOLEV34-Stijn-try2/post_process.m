global data;
x_optim_n= [0.0043
    0.0138
    0.8085
    0.9424
    0.6271
    0.6137
    0.2760
    0.4212
    0.4643
    0.5991
    0.5351
    0.5354
    0.4591
    0.4288
    0.5309
    0.5311
    0.4554
    0.4094
    0.4413
    0.1273
    0.6263
    0.7717
    0.5730
    0.5188
    0.6154
    0.6492
    0.6536
    0.6285
    0.6971
    0.5303
    0.5520];

x_optim = x_optim_n.*(data.ub-data.lb)+data.lb;
fval =  0.963479;

W_TO_max_optim = data.W_TO_max_ref*fval;


taper = x_optim(2)/x_optim(1);
mac = x_optim(1)*2/3*(1 + taper + taper^2)/(1 + taper);
Re = data.rho*data.V*mac/(1.422*10^-5);

W_fuel_optim = 7569* 9.80665;
W_str_optim = 768 * 9.80665;

p = 0;
while p < 3
[x_section_optim, y_section_optim, z_section_optim, c_section_optim, twist_section_optim, S_optim, y_85_optim, c_85_optim] = geometry_function(x_optim,  data.y_kink,  data.dihedral);
[Aur_optim, Alr_optim, Aut_optim, Alt_optim] = GeomtoCST();
[Xtur,Xtlr,Xtut,Xtlt,Xtuk,Xtlk,Xtu85,Xtl85] = CSTtoGeom(Aur_optim, Alr_optim, Aut_optim, Alt_optim,  data.y_kink/x_optim(3),c_section_optim,c_85_optim);
[CL_max_optim] = DPC_function(data.W_AW, W_fuel_optim,  W_str_optim, S_optim,  n_max);
[cl_distribution_optim_max, cm_distribution_optim_max, Y_distribution_optim_max, chord_distribution_optim_max, CLwing_optim_max, CDwing_optim_max] = Q3D_function(x_section_optim, y_section_optim, z_section_optim, c_section_optim, twist_section_optim, [Aur_optim; Alr_optim], [Aut_optim, Alt_optim], 0, 150, CL_max_optim,  data.V,  data.rho,  data.alt,  Re,  data.M);
EMWET_function(x_optim, cl_distribution_optim_max, cm_distribution_optim_max, Y_distribution_optim_max, chord_distribution_optim_max,  W_TO_max_optim,  W_fuel_optim,  n_max, S_optim, x_section_optim, y_section_optim, z_section_optim, c_section_optim,  data.rho,  data.V);
EMWET TUPOLEVaircraft;
filetext = fileread('TUPOLEVaircraft.weight');
W_str_kg = str2double(regexp(filetext, '(?<=Wing total weight[^0-9]*)[0-9]*\.?[0-9]+', 'match'));
W_str_optim = W_str_kg*9.80665;
[CL_design_optim] = DPC_function(data.W_AW, W_fuel_optim,  W_str_optim, S_optim, 1);
[cl_distribution_optim, cm_distribution_optim, Y_distribution_optim, chord_distribution_optim, CLwing_optim, CDwing_optim] = Q3D_function(x_section_optim, y_section_optim, z_section_optim, c_section_optim, twist_section_optim, [Aur_optim; Alr_optim], [Aut_optim, Alt_optim], 1, 1000, CL_design_optim,  data.V,  data.rho,  data.alt,  Re,  data.M);
[W_fuel_optim] = Brequet_function(CLwing_optim, CDwing_optim,  data.W_AW, W_fuel_optim, W_str_optim, data.D_AW, S_optim,  data.rho,  data.V);
p = p + 1;
end

[data.W_AW/9.80665, W_str_kg, W_fuel_optim/9.80665, W_TO_max_optim/9.80665]
