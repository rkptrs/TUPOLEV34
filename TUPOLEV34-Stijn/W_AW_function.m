function [W_AW] = W_AW_function(x_1, y_kink, dihedral, )

[x_section, y_section, z_section, c_section, twist_section, S, y_85, c_85] = geometry_function(x_1, y_kink, dihedral);

[cl_distribution, cm_distribution, Y_distribution, chord_distribution, CLwing, CDwing] = Q3D_function(x_section, y_section, z_section, c_section, twist_section, CST_root, CST_kink, CST_tip, Viscosity_coefficient, MaxIterIndex, CL_design, V, rho, alt, Re, M)
