close all;
clc;
clear all;

%First run the part for the reference aircraft for initial CST coefficients
%and A-W contribution

s = load('reffile.mat');

x_1_ref = [s.C_root, s.C_tip, s.b/2, s.Lambda_1, s.Lambda_2, s.Incidence_root, s.Incidence_tip];
[x_section_ref, y_section_ref, z_section_ref, c_section_ref, twist_section_ref, S_ref, y_85_ref, c_85_ref] = geometry_function(x_1_ref, s.y_kink, s.dihedral);
[Aur_ref, Alr_ref, Aut_ref, Alt_ref] = GeomtoCST();
[Xtur,Xtlr,Xtut,Xtlt,Xtuk,Xtlk,Xtu85,Xtl85] = CSTtoGeom(Aur_ref, Alr_ref, Aut_ref, Alt_ref, s.y_kink/x_1_ref(3));
[CL_design_ref] = DPC_function(0, 0, s.W_TO_max, S_ref, s.n_max);
[cl_distribution_ref, cm_distribution_ref, Y_distribution_ref, chord_distribution_ref, CLwing_ref, CDwing_ref] = Q3D_function(x_section_ref, y_section_ref, z_section_ref, c_section_ref, twist_section_ref, [Aur_ref; Alr_ref], [Aut_ref, Alt_ref], 0, 150, CL_design_ref, s.V, s.rho, s.alt, s.Re, s.M);
[W_str_ref] = EMWET_function(x_1_ref, cl_distribution_ref, cm_distribution_ref, Y_distribution_ref, chord_distribution_ref, s.W_TO_max, s.W_TO_max - s.W_ZFW, s.n_max, S_ref, x_section_ref, y_section_ref, z_section_ref, c_section_ref, s.rho, s.V);
[CL_design_ref] = DPC_function(0, 0, s.W_TO_max, S_ref, 1);
[cl_distribution_ref, cm_distribution_ref, Y_distribution_ref, chord_distribution_ref, CLwing_ref, CDwing_ref] = Q3D_function(x_section_ref, y_section_ref, z_section_ref, c_section_ref, twist_section_ref, [Aur_ref; Alr_ref], [Aut_ref, Alt_ref], 1, 150, CL_design_ref, s.V, s.rho, s.alt, s.Re, s.M);
[W_fuel_ref] = Brequet_function(16, 1, s.W_TO_max, 0, 0, 0, S_ref, s.rho, s.V);
W_AW = s.W_TO_max - W_str_ref - W_fuel_ref;   %Newton
C_D_AW = CLwing_ref/16 - CDwing_ref;
D_AW = 0.5*s.rho*S_ref*s.V^2*C_D_AW;
save('reffile.mat', 'W_AW', 'D_AW', 'W_str_ref', 'W_fuel_ref', '-append')
x0 = [x_1_ref'; Aur_ref; Alr_ref; Aut_ref; Alt_ref];

[x, fval] = fmincon(@(x) MDA(x), x0, [], [], [], [], x0*0.9, x0*1.1, @(x) constraints(x)); 

