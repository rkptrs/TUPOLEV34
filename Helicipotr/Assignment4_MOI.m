%% Helicopter parameters

mass = 17650 * 0.45359237; %kg

% geometric parameters in m
l_fus = 7.7;
D_fus = 2.3;
l_eng = 2.3;
D_eng = 0.77;
blade_radius = 24 * 0.3048; 
blade_radius_tail = 4.6*0.3048;
l_fft_shaft = 5.4;
D_fft_shaft = 1.1;
l_vt = 2.8;
w_vt = 1.0;
l_ht = 3.4;
w_ht = 0.92;
l_wing = 1.8;
w_wing = 1.0;

%% Center of gravity
% component lists are: [fuselage; engines; main rotor; tail rotor; ...
%                       fuselage-to-tail-shaft; vertical tail; ...
%                       horizontal tail; wings]

% masses as percentage of MTOM
masses = [55; 15; 8; 2; 8; 3; 3; 6];
masses_kg = masses./100.*mass;

% cog positions of each of the components as [x, z]
cog_positions = [4.7614 1.7782;
                6.4546 2.1872;
                5.0568 3.6221;
                14.068 2.8698;
                10.284 1.0995;
                13.921 2.1455;
                14.136 1.0400;
                5.4887 1.4338];

% cog of the helicopter
x_cg = sum(cog_positions(:,1).*masses) / sum(masses);
z_cg = sum(cog_positions(:,2).*masses) / sum(masses);
cog_heli = [x_cg, z_cg];

%% Mass moment of inertia about y

% MMOI due to geometric shapes
mmoi_fus = (masses_kg(1) * (D_fus / 2) / 5) * (1 + ((l_fus/2)/(D_fus/2))^2);
mmoi_engines = 1 / 12 * masses_kg(2) * (3 * (D_eng/2)^2 + l_eng^2);
mmoi_main_rotor = 1 / 4 * masses_kg(3) * blade_radius^2;
mmoi_tail_rotor = 1 / 2 * masses_kg(4) * blade_radius_tail^2;
mmoi_ftt_shaft = 1 / 12 * masses_kg(5) * (3 * (D_fft_shaft/2)^2 + l_fft_shaft^2);
mmoi_vt = 1 / 12 * masses_kg(6) * (w_vt^2 + l_vt^2);
mmoi_ht = 1 / 12 * masses_kg(7) * w_ht^2;
mmoi_wings = 1 / 12 * masses_kg(8) * w_wing^2;

MMOI_geom = [mmoi_fus; mmoi_engines; mmoi_main_rotor; mmoi_tail_rotor; mmoi_ftt_shaft; mmoi_vt; mmoi_ht; mmoi_wings];

% MMOI due to parallel axis theorem I = m*d^2
distances = sqrt((cog_positions(:, 1) - x_cg).^2 + (cog_positions(:, 2) - z_cg).^2);
MMOI_PA = masses_kg .* distances.^2;

% total MMOI
MMOI_total = MMOI_geom + MMOI_PA;
MMOI_yy = sum(MMOI_total);


