close all;
clc;
clear all;

C_root = 5.6;
C_tip = 0.94;
b = 28;
Lambda_1 = 25;
Lambda_2 = 25;
Incidence_root = 0;
Incidence_tip = -1;

y_kink = 5;
dihedral = 6;
W_TO_max = 46500*9.80665;
W_ZFW = 41050*9.80665;
MAC = 3.53;
n_max = 2.5;
V = 222.24;
rho = 0.356804;
alt = 11125.2;
Re = rho*V*3.53/(1.422*10^-5);
speedofsound = 295;
M = V/speedofsound;

xold = [0;0];

no_convergence_value = 1;

x_1_ref = [ C_root,  C_tip,  b/2,  Lambda_1,  Lambda_2,  Incidence_root,  Incidence_tip];
[x_section_ref, y_section_ref, z_section_ref, c_section_ref, twist_section_ref, S_ref, y_85_ref, c_85_ref] = geometry_function(x_1_ref,  y_kink,  dihedral);
[Aur_ref, Alr_ref, Aut_ref, Alt_ref] = GeomtoCST();
[Xtur,Xtlr,Xtut,Xtlt,Xtuk,Xtlk,Xtu85,Xtl85] = CSTtoGeom(Aur_ref, Alr_ref, Aut_ref, Alt_ref,  y_kink/x_1_ref(3),c_section_ref,c_85_ref);
[CL_max_ref] = DPC_function(0, 0,  W_TO_max, S_ref,  n_max);
[cl_distribution_ref, cm_distribution_ref, Y_distribution_ref, chord_distribution_ref, CLwing_ref, CDwing_ref] = Q3D_function(x_section_ref, y_section_ref, z_section_ref, c_section_ref, twist_section_ref, [Aur_ref; Alr_ref], [Aut_ref, Alt_ref], 0, 150, CL_max_ref,  V,  rho,  alt,  Re,  M);
EMWET_function(x_1_ref, cl_distribution_ref, cm_distribution_ref, Y_distribution_ref, chord_distribution_ref,  W_TO_max,  W_TO_max -  W_ZFW,  n_max, S_ref, x_section_ref, y_section_ref, z_section_ref, c_section_ref,  rho,  V);
EMWET TUPOLEVaircraft;
filetext = fileread('TUPOLEVaircraft.weight');
W_str_kg = str2double(regexp(filetext, '(?<=Wing total weight[^0-9]*)[0-9]*\.?[0-9]+', 'match'));
W_str_ref = W_str_kg*9.80665;
[CL_design_ref] = DPC_function(0, 0,  W_TO_max, S_ref, 1);
[cl_distribution_ref, cm_distribution_ref, Y_distribution_ref, chord_distribution_ref, CLwing_ref, CDwing_ref] = Q3D_function(x_section_ref, y_section_ref, z_section_ref, c_section_ref, twist_section_ref, [Aur_ref; Alr_ref], [Aut_ref, Alt_ref], 1, 1000, CL_design_ref,  V,  rho,  alt,  Re,  M);
[W_fuel_ref] = Brequet_function(16, 1,  W_TO_max, 0, 0, 0, S_ref,  rho,  V);
W_AW =  W_TO_max - W_str_ref - W_fuel_ref;   %Newton
C_D_AW = CLwing_ref/16 - CDwing_ref;
D_AW = 0.5* rho*S_ref* V^2*C_D_AW;
x0 = [x_1_ref'; Aur_ref; Alr_ref; Aut_ref; Alt_ref];
lb = [3,  0.2,  10,  0.5,   0.5,   -4,    -4, Aur_ref' - 0.3*abs(Aur_ref)',Alr_ref' - 0.3*abs(Alr_ref)',Aut_ref' - 0.3*abs(Aut_ref)',Alt_ref' - 0.3*abs(Alt_ref)']';
ub = [20, 20, 30, 45, 45, 4, 4, Aur_ref' + 0.3*abs(Aur_ref)',Alr_ref' + 0.3*abs(Alr_ref)',Aut_ref' + 0.3*abs(Aut_ref)',Alt_ref' + 0.3*abs(Alt_ref)']';
lbtest = [5.7,  0.8,  13,  24,   24,   -1.1,    -0.9, Aur_ref' - 0.05*abs(Aur_ref)',Alr_ref' - 0.05*abs(Alr_ref)',Aut_ref' - 0.05*abs(Aut_ref)',Alt_ref' - 0.05*abs(Alt_ref)']';
ubtest = [5.8, 1.1, 15, 26, 26, -0.1, 0.1, Aur_ref' + 0.05*abs(Aur_ref)',Alr_ref' + 0.05*abs(Alr_ref)',Aut_ref' + 0.05*abs(Aut_ref)',Alt_ref' + 0.05*abs(Alt_ref)']';

W_str_mda = W_str_ref;
W_fuel_mda = W_fuel_ref;
W_TO_max_mda =  W_TO_max;


f_tank = 0.93;
rho_fuel = 0.81715*10^3;
V_fuel_ref = W_fuel_ref/9.80665/rho_fuel;
WS_original = 558.72;

[V_tank_ref] = fuel_volume(Xtur, Xtuk, Xtu85, Xtlr, Xtlk, Xtl85, y_85_ref, c_85_ref, y_section_ref, c_section_ref);
[WS_new_ref] = Wing_Loading_function(W_TO_max, S_ref);
c_volume_ref = (V_fuel_ref - V_tank_ref*f_tank);
c_wingloading_ref = (WS_new_ref - WS_original);

global data;
data.C_root = C_root;
data.C_tip = C_tip;
data.b = b;
data.Lambda_1 = Lambda_1;
data.Lambda_2 = Lambda_2;
data.Incidence_root = Incidence_root;
data.Incidence_tip = Incidence_tip;
data.y_kink = y_kink;
data.dihedral = dihedral;
data.W_TO_max = W_TO_max;
data.W_TO_max_ref = W_TO_max;
data.W_ZFW = W_ZFW; 
data.MAC = MAC;
data.n_max = n_max;
data.V = V;
data.rho = rho;
data.alt = alt;
data.Re = Re;
data.speedofsound = speedofsound;
data.M = M;
data.W_AW = W_AW;
data.D_AW = D_AW;
data.W_str_mda = W_str_mda;
data.W_fuel_mda = W_fuel_mda;
data.W_TO_max_mda = W_TO_max_mda;
data.lb = lb;
data.ub = ub;
data.lbtest = lbtest;
data.ubtest = ubtest;
data.x0 = x0;
data.xold = xold;
data.c_volume_ref = c_volume_ref;
data.c_wingloading_ref = c_wingloading_ref;
data.no_convergence_value = no_convergence_value;
data.clist = [];

