%helicopters

mass = 17650 * 0.45359237; %lbs to kg
weight = mass * 9.80665;
blade_radius = 24 * 0.3048; %ft to m
blade_radius_tail = 4.6*0.3048;
omega_r = 221.2848;
omega_r_tail = 677*0.3048;
rho = 1.225;
solidity = 0.092;
solidity_tail = 0.231;

l_tail = 9;
%% hover
FM = 0.6;
CL_medium = 0.55;
CD_p_medium = 0.02;
CD_p_medium_tail = 0.02;

k = 1.15;
k_tail =1.4;
P_id = weight * sqrt(weight/(2*rho*pi*blade_radius^2));
v_i = sqrt(weight/(2*rho*pi*blade_radius^2));

%ACT
P_hov_ACT = P_id / FM;

%BEM
CT_BEM = solidity*CL_medium/6.6;
T_BEM = CT_BEM*(rho*pi*blade_radius^2*omega_r^2);
P_i_BEM = k*T_BEM*v_i;
P_p_BEM = solidity*CD_p_medium/8*rho*omega_r^3*pi*blade_radius^2;
P_p_BEM_tail = solidity_tail*CD_p_medium_tail/8*rho*omega_r_tail^3*pi*blade_radius_tail^2;
P_hov_BEM = P_i_BEM + P_p_BEM;
FM_BEM = P_i_BEM/P_hov_BEM;

T_tail = P_hov_BEM/(omega_r/blade_radius)/l_tail;
v_i_tail = sqrt(T_tail/(2*rho*pi*blade_radius_tail^2));
P_i_tail = 1.1*k_tail*T_tail*v_i_tail;
%% forward flight
V = 0:100;

mu = V/omega_r;
mu_tail = V/omega_r_tail;

%Main rotor:
P_ind_fw_rotor = P_i_BEM - (V)/V(end)*(0.75*P_i_BEM);
P_profiledrag_fw_rotor = P_p_BEM.*(1+mu.^2);
P_drag_fw_rotor = P_p_BEM.*mu.^2.*2;

P_p_fw_rotor = P_p_BEM.*(1+4.65*mu.^2);

%Tail rotor
P_ind_fw_tail = P_i_tail - (V)/V(end)*(0.75*P_i_tail);
P_profiledrag_fw_tail = P_p_BEM_tail.*(1+mu_tail.^2);
P_drag_fw_tail = P_p_BEM_tail.*mu_tail.^2.*2;

P_p_fw_tail = P_p_BEM_tail.*(1+4.65*mu_tail.^2);
P_total_tail = P_ind_fw_tail + P_profiledrag_fw_tail + P_drag_fw_tail;

P_total = P_ind_fw_rotor + P_profiledrag_fw_rotor + P_drag_fw_rotor + P_total_tail;

plot(V, P_ind_fw_rotor); hold on; plot(V, P_profiledrag_fw_rotor); hold on; plot(V, P_drag_fw_rotor); hold on; plot(V, P_total);hold on; plot(V, P_total_tail)

