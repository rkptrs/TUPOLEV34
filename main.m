
%First run the part for the reference aircraft for initial CST coefficients
%and A-W contribution

s = load('reffile.mat');

x_1_ref = [s.C_root, s.C_tip, s.b/2, s.Lambda_1, s.Lambda_2, s.Incidence_root, s.Incidence_tip];
[x_section, y_section, z_section, c_section, twist_section, S, y_85, c_85] = geometry_function(x_1_ref, s.y_kink, s.dihedral);
