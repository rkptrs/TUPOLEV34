function [W_str] = EMWET_function(x_1, cl_distribution, cm_distribution, Y_distribution, chord_distribution, W_TO_max, W_fuel, n_max, S, x_section, y_section, z_section, c_section, rho, V)

% EMWET 
fid = fopen('f50.init', 'wt');
    fprintf(fid,'%g %g\n', 20820, 18600);
    fprintf(fid,'%g\n', 2.5);
    fprintf(fid,'%g %g %g %g\n', 28*((3.5 + 3.5*0.25)/2), 28, 2, 2);
    fprintf(fid,'%g %s\n', 0, 'e553');
    fprintf(fid,'%g %s\n', 1, 'e553');
    fprintf(fid,'%g %g %g %g %g %g\n', 3.5, 0, 0, 0, 0.2, 0.8);
    fprintf(fid,'%g %g %g %g %g %g\n', 3.5*0.25, 1.22484, 14, 0, 0.2, 0.8);
    fprintf(fid,'%g %g\n', 0.1, 0.7);
    fprintf(fid,'%g\n', 1);
    fprintf(fid,'%g %g\n', 0.25, 1200);
    fprintf(fid,'%g %g %g %g\n', 7.1e+010, 2800, 4.8e+08, 4.6e+08);
    fprintf(fid,'%g %g %g %g\n', 7.1e+010, 2800, 4.8e+08, 4.6e+08);
    fprintf(fid,'%g %g %g %g\n', 7.1e+010, 2800, 4.8e+08, 4.6e+08);
    fprintf(fid,'%g %g %g %g\n', 7.1e+010, 2800, 4.8e+08, 4.6e+08);
    fprintf(fid,'%g %g\n', 0.96, 0.5);
    fprintf(fid,'%g\n', 1);   
fclose(fid);


taper = x_1(2)/x_1(1);
mac = x_1(1)*2/3*(1 + taper + taper^2)/(1 + taper);
q = 0.5*rho*V^2;
chord = chord_distribution;
moment = mac .* chord .* cm_distribution .* q;
lift = n_max*Res.Wing.chord .* cl_distribution .* q;

fid = fopen('f50.load', 'wt');
    for i = 1:length(lift)
        fprintf(fid,'%g %g %g\n', Y_distribution(i)/14, lift(i), moment(i));
    end
fclose(fid);

EMWET f50
    
    