

%% Helicopter parameters

mass = 17650 * 0.45359237; %lbs to kg
weight = mass * 9.80665;
blade_radius = 24 * 0.3048; %ft to m
omega_r = 221.2848;
rho = 1.225;
CD_S = 3.0;
solidity = 0.092;
cl_alpha = 0.115;



%% Other variables
max_speed = 100;
V = 0:1:max_speed;

mu = V./(omega_r);
D_fus = 0.5*rho*CD_S.*V.^2;
W = weight;
T = sqrt(W^2 + D_fus.^2);
CT = T ./ (rho * (omega_r)^2 * pi * blade_radius^2);



%% Calculate lambda_i
lambda_i = 0:10^-7:1;

CT_glauert_list = [];
lambda_i_list = [];
for i = 1:length(V)
    for j = 1:length(lambda_i)
        CT_glauert = 2 * lambda_i(j) * sqrt((mu(i) * cos(D_fus(i)/W))^2 + (mu(i) * sin(D_fus(i)/W) + lambda_i(j))^2);
        if abs(CT_glauert - CT(i)) < 10^-8
            CT_glauert_list = [CT_glauert_list, CT_glauert];
            lambda_i_list = [lambda_i_list, lambda_i(j)];
            break
        end
    end
end

%% Calculating the trim angles
theta_c_list = [];
theta_0_list = [];
for i = 1:length(V)
    A = [(1 + 3/2*mu(i)^2), (-8/3*mu(i));...
        (-mu(i)), (2/3 + mu(i)^2)];
    B = [(-2 * mu(i)^2 * D_fus(i)/W - 2* mu(i) * lambda_i(i));...
        (4 / solidity * CT(i) / cl_alpha + mu(i) * D_fus(i) / W + lambda_i(i))];
    C = A\B;
    theta_c_list = [theta_c_list C(1)];
    theta_0_list = [theta_0_list C(2)];
end

%% Plotting
figure()
plot(V, theta_c_list); hold on; plot(V, theta_0_list); hold off;
legend('Cyclic pitch angle \theta_c', 'Collective pitch angle \theta_0')
ylabel('Trim angle [deg]')
xlabel('Level flight velocity [m/s]')
grid on
 