clear all;
clc;
close all;

%% Helicopter parameters

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


%% Hover
FM = 0.6;
CL_medium = 0.55;
CD_p_medium = 0.02;
CD_p_medium_tail = 0.02;

k = 1.15;
k_tail =1.4;
P_id = weight * sqrt(weight/(2*rho*pi*blade_radius^2));
v_i_hov = sqrt(weight/(2*rho*pi*blade_radius^2));

%ACT
P_hov_ACT = P_id / FM;

%BEM
CT_BEM = solidity*CL_medium/6.6;
T_BEM = CT_BEM*(rho*pi*blade_radius^2*omega_r^2);
P_i_BEM = k*T_BEM*v_i_hov;
P_p_BEM = solidity*CD_p_medium/8*rho*omega_r^3*pi*blade_radius^2;
P_p_BEM_tail = solidity_tail*CD_p_medium_tail/8*rho*omega_r_tail^3*pi*blade_radius_tail^2;
P_hov_BEM = P_i_BEM + P_p_BEM;
FM_BEM = P_i_BEM/P_hov_BEM;


%% Forward flight
V = 0:0.01:100;

mu = V/omega_r;
mu_tail = V/omega_r_tail;

%Main rotor:
%Induced power
V_bar = V./v_i_hov;
v_i_bar = sqrt(-V_bar.^2/2 + sqrt(V_bar.^4/4 + 1));
P_ind_fw_rotor = k*T_BEM*v_i_bar*v_i_hov;

%Profile drag power
H_0 = solidity*CD_p_medium/4*rho*omega_r^2*pi*blade_radius^2*mu;
P_d = H_0.*omega_r.*mu;
P_p = solidity*CD_p_medium/8*rho*omega_r^3*pi*blade_radius^2*(1+mu.^2);
P_profiledrag_fw_rotor = P_d + P_p;
P_profiledrag_fw_rotor_Bennet = solidity*CD_p_medium/8*rho*omega_r^3*pi*blade_radius^2*(1+4.65*mu.^2);

%Parasite drag power
sum_CD_S = 3;   %assuming the apache is not aerodynamically refined (Figure 5.11)
P_par = sum_CD_S*0.5*rho*V.^3;

%Total main rotor power
P_total_rotor = P_ind_fw_rotor + P_profiledrag_fw_rotor_Bennet + P_par;

%Tail rotor power
T_tail = P_hov_BEM/(omega_r/blade_radius)/l_tail;
v_i_tail = sqrt(T_tail/(2*rho*pi*blade_radius_tail^2));
P_i_tail = 1.1*k_tail*T_tail*v_i_tail;

V_bar_tail = V./v_i_tail;
v_i_bar_tail = sqrt(-V_bar_tail.^2/2 + sqrt(V_bar_tail.^4/4 + 1));
P_ind_fw_tail = P_i_tail*v_i_bar_tail;

P_profiledrag_fw_tail = P_p_BEM_tail.*(1+mu_tail.^2);
P_drag_fw_tail = P_p_BEM_tail.*mu_tail.^2.*2;
P_p_fw_tail = P_p_BEM_tail.*(1+4.65*mu_tail.^2);
% P_total_tail = P_ind_fw_tail + P_profiledrag_fw_tail + P_drag_fw_tail;
P_total_tail = P_ind_fw_tail + P_p_fw_tail;


%Total required power
P_total = P_total_rotor + P_total_tail;


% %% plotting total required power
figure()
plot(V, P_ind_fw_rotor/1000); hold on; plot(V, P_profiledrag_fw_rotor/1000);hold on; plot(V, P_par/1000);hold on; plot(V, P_total_tail/1000); hold on; plot(V, P_total/1000)
grid on;
legend('Induced power', 'Profile drag power', 'Parasite drag power', 'Tail power', 'Total required power', 'Location', 'best')
xlabel('Helicopter velocity [m/s]')
ylabel('Power [kW]')



%% Power speed calculations
%Speed for maximum endurance
P_r = P_total;
index_endurance = find(P_r == min(P_r));
V_endurance = V(index_endurance);

%Speed for maximum range
P_r_slope = [];
for i = 2:length(P_r)
    P_r_slope_i = (P_r(i)-P_r(i-1))/(V(i)-V(i-1));
    P_r_slope = [P_r_slope; P_r_slope_i];
end
P_slope0 = P_r./V;

point = InterX([V(2:end); P_slope0(2:end)], [V(2:end); P_r_slope']);
V_range = point(1);
P_r_range = P_r(find(V == round(V_range,2)));

%plotting
figure()
plot(V, P_r/1000); hold on; plot([0, V_range], [0, P_r_range/1000]); hold on; scatter(V_endurance, min(P_r)/1000); hold on; scatter(V_range, P_r_range/1000); hold on;
grid on;
legend('Total required power', 'Tangency line', 'Point of maximum endurance', 'Point of maximum range', 'Location', 'best')
xlabel('Helicopter velocity [m/s]')
ylabel('Power [kW]')

