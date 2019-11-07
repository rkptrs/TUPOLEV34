function [W_str] = EMWET_function(x_1, cl_distribution, cm_distribution, Y_distribution, chord_distribution, W_TO_max, W_fuel, n_max, S, x_section, y_section, z_section, c_section, rho, V)

W_ZFW = (W_TO_max - W_fuel)/9.80665; %kg

% EMWET 
fid = fopen('tu334.init', 'wt');
    fprintf(fid,'%g %g\n', W_TO_max/9.80665, W_ZFW);
    fprintf(fid,'%g\n', n_max);
    fprintf(fid,'%g %g %g %g\n', S, x_1(3)*2, 3, 2);
    fprintf(fid,'%g %s\n', 0, 'airfoilroot');
    fprintf(fid,'%g %s\n', 1, 'airfoiltip');
    fprintf(fid,'%g %g %g %g %g %g\n', c_section(1), x_section(1), y_section(1), z_section(1), 0.15, 0.8);
    fprintf(fid,'%g %g %g %g %g %g\n', c_section(2), x_section(2), y_section(2), z_section(2), 0.25, 0.65);
    fprintf(fid,'%g %g %g %g %g %g\n', c_section(3), x_section(3), y_section(3), z_section(3), 0.2, 0.6);
    fprintf(fid,'%g %g\n', 0, 0.85);
    fprintf(fid,'%g\n', 0);
    fprintf(fid,'%g %g %g %g\n', 7e+010, 2800, 2.95e+08, 2.95e+08);
    fprintf(fid,'%g %g %g %g\n', 7e+010, 2800, 2.95e+08, 2.95e+08);
    fprintf(fid,'%g %g %g %g\n', 7e+010, 2800, 2.95e+08, 2.95e+08);
    fprintf(fid,'%g %g %g %g\n', 7e+010, 2800, 2.95e+08, 2.95e+08);
    fprintf(fid,'%g %g\n', 0.96, 0.5);
    fprintf(fid,'%g\n', 1);   
fclose(fid);


taper = x_1(2)/x_1(1);
mac = x_1(1)*2/3*(1 + taper + taper^2)/(1 + taper);
q = 0.5*rho*V^2;
chord = chord_distribution;
moment = mac .* chord .* cm_distribution .* q;
lift = n_max*chord_distribution .* cl_distribution .* q;

fid = fopen('tu334.load', 'wt');
    for i = 1:length(lift)
        fprintf(fid,'%g %g %g\n', Y_distribution(i)/14, lift(i), moment(i));
    end
fclose(fid);

EMWET tu334;

filetext = fileread('tu334.weight');
W_str_kg = str2double(regexp(filetext, '(?<=Wing total weight[^0-9]*)[0-9]*\.?[0-9]+', 'match'));

W_str = W_str_kg * 9.80665;    
    