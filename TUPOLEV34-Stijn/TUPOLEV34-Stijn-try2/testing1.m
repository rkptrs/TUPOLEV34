
global data;
%x = x0n.*(data.ub-data.lb)+data.lb;
x(1:7)  = [2;  0.2;  10;  30;   4;   -5;    -5];%check the input of the loop

% lb = x(1:7)  = [3;  0.2;  10;  0.5;   0.5;   -5;    -5];%check the input of the loop
% ub = x(1:7)  = [20;  20;  30;  45;   45;   5;    5];%check the input of the loop
x = data.lb;
% 1.0000
%     0.2000
%     7.5000
%     1.0000
%    50.0000
%     4.5957
%    -5.0000

x = [4.3260
    0.5774
   15.8125
   34.8000
   34.8000
    1.2495
    1.9400
    0.1765
    0.0601
    0.3338
    0.0671
    0.2104
    0.2878
   -0.2806
   -0.2035
   -0.0584
   -0.3603
    0.0555
    0.2458
    0.1058
    0.0598
    0.1998
    0.0669
    0.2079
    0.2848
   -0.1021
   -0.0742
   -0.0210
   -0.2164
    0.0552
    0.2431];

W_fuel = 7.8384e+04;

   
W_str = data.W_str_mda;
[x_section, y_section, z_section, c_section, twist_section, S, y_85, c_85] = geometry_function(x(1:7), data.y_kink, data.dihedral);

W_TO_max = 4.3860e+05;

[Xtur,Xtlr,Xtut,Xtlt,Xtuk,Xtlk,Xtu85,Xtl85] = CSTtoGeom(x(8:13), x(14:19), x(20:25), x(26:31), data.y_kink/x(3),c_section,c_85);

[CL_design_max] = DPC_function(data.W_AW, W_fuel, W_str, S, data.n_max);
[cl_distribution, cm_distribution, Y_distribution, chord_distribution, CLwing, CDwing] = Q3D_function(x_section, y_section, z_section, c_section, twist_section, x(8:19), x(20:31), 1, 1500, CL_design_max, data.V, data.rho, data.alt, data.Re, data.M);

EMWET_function(x(1:7), cl_distribution, cm_distribution, Y_distribution, chord_distribution, W_TO_max, W_fuel, data.n_max, S, x_section, y_section, z_section, c_section, data.rho, data.V);
EMWET TUPOLEVaircraft
