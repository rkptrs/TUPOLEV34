

x_1_ref = [s.C_root, s.C_tip, s.b/2, s.Lambda_1, s.Lambda_2, s.Incidence_root, s.Incidence_tip];
[x_section_ref, y_section_ref, z_section_ref, c_section_ref, twist_section_ref, S_ref, y_85_ref, c_85_ref] = geometry_function(x_1_ref, s.y_kink, s.dihedral);
[Aur_ref, Alr_ref, Aut_ref, Alt_ref] = GeomtoCST();
[Xtur,Xtlr,Xtut,Xtlt,Xtuk,Xtlk,Xtu85,Xtl85] = CSTtoGeom(Aur_ref, Alr_ref, Aut_ref, Alt_ref, s.y_kink/x_1_ref(3),c_section_ref,c_85_ref);
[CL_design_ref] = DPC_function(0, 0, s.W_TO_max, S_ref, s.n_max);
[cl_distribution_ref, cm_distribution_ref, Y_distribution_ref, chord_distribution_ref, CLwing_ref, CDwing_ref] = Q3D_function(x_section_ref, y_section_ref, z_section_ref, c_section_ref, twist_section_ref, [Aur_ref; Alr_ref], [Aut_ref, Alt_ref], 0, 150, CL_design_ref, s.V, s.rho, s.alt, s.Re, s.M);
[W_str_ref] = EMWET_function(x_1_ref, cl_distribution_ref, cm_distribution_ref, Y_distribution_ref, chord_distribution_ref, s.W_TO_max, s.W_TO_max - s.W_ZFW, s.n_max, S_ref, x_section_ref, y_section_ref, z_section_ref, c_section_ref, s.rho, s.V);
[CL_design_ref] = DPC_function(0, 0, s.W_TO_max, S_ref, 1);
[cl_distribution_ref, cm_distribution_ref, Y_distribution_ref, chord_distribution_ref, CLwing_ref, CDwing_ref] = Q3D_function(x_section_ref, y_section_ref, z_section_ref, c_section_ref, twist_section_ref, [Aur_ref; Alr_ref], [Aut_ref, Alt_ref], 1, 150, CL_design_ref, s.V, s.rho, s.alt, s.Re, s.M);
[W_fuel_ref] = Brequet_function(16, 1, s.W_TO_max, 0, 0, 0, S_ref, s.rho, s.V);
W_AW = s.W_TO_max - W_str_ref - W_fuel_ref;   %Newton
C_D_AW = CLwing_ref/16 - CDwing_ref;
D_AW = 0.5*s.rho*S_ref*s.V^2*C_D_AW;
x0 = [x_1_ref'; Aur_ref; Alr_ref; Aut_ref; Alt_ref];
lb = [0.5,0.1,5,0.5,0.5,-5,-5, Aur_ref' - 0.5*abs(Aur_ref)',Alr_ref' - 0.5*abs(Alr_ref)',Aut_ref' - 0.5*abs(Aut_ref)',Alt_ref' - 0.5*abs(Alt_ref)']';
ub = [20,20,30,60,60,5,5,Aur_ref' + 0.5*abs(Aur_ref)',Alr_ref' + 0.5*abs(Alr_ref)',Aut_ref' + 0.5*abs(Aut_ref)',Alt_ref' + 0.5*abs(Alt_ref)']';

W_str_mda = W_str_ref;
W_fuel_mda = W_fuel_ref;
W_TO_max_mda = s.W_TO_max;

save('reffile.mat', 'W_AW', 'D_AW', 'W_str_mda', 'W_fuel_mda', 'W_TO_max_mda', 'lb', 'ub', 'x0', '-append');