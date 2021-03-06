%input = design vector x_1 (vector), kink location y_kink and dihedral angle.
function [x_section, y_section, z_section, c_section, twist_section, S, y_85, c_85] = geometry_function(x_1, y_kink, dihedral)
    %input for q3d, emwet, wing loading, fuel volume
    kink_percentage = y_kink/(x_1(3));
    c_kink = x_1(1) - kink_percentage*(x_1(1) - x_1(2));
    twist_kink = x_1(6) - kink_percentage*(x_1(6) - x_1(7));
    x_kink = sind(x_1(4))*y_kink;
    z_kink = sind(dihedral)*y_kink;
    
    y_tip = x_1(3);
    x_tip = x_kink + sind(x_1(5))*(y_tip - y_kink);
    z_tip = sind(dihedral)*y_tip;
    
    x_section = [0, x_kink, x_tip];
    y_section = [0, y_kink, y_tip];
    z_section = [0, z_kink, z_tip];
    c_section = [x_1(1), c_kink, x_1(2)];
    twist_section = [x_1(6), twist_kink, x_1(7)];
    
    S1 = y_kink*(x_1(1) + c_kink)/2;
    S2 = (y_tip - y_kink)*(x_1(2) + c_kink)/2;
    S = 2*(S1 + S2);
    
    %extra for fuel volume
    y_85 = x_1(3)*0.85;
    c_85 = x_1(1) - 0.85*(x_1(1) - x_1(2));
end