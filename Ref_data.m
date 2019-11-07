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

save('reffile.mat','C_root','C_tip', 'b', 'Lambda_1', 'Lambda_2', 'Incidence_root', 'Incidence_tip', 'y_kink', 'dihedral', 'W_TO_max', 'MAC', 'n_max', 'V', 'rho', 'alt', 'Re', 'speedofsound', 'M' )
