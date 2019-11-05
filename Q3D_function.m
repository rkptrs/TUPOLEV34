%input = output of geometryfunction and CST function
% Viscosity_coefficient = 1 if viscous, 0 if inviscid analysis
% MaxIterIndex = maximum convergence iterations for viscous analysis
% CL_design from DPC function
% flight conditions V, rho, alt, Re, M

function [cl_distribution, cm_distribution, Y_distribution, chord_distribution, CLwing, CDwing] = Q3D_function(x_section, y_section, z_section, c_section, twist_section, CST_root, CST_kink, CST_tip, Viscosity_coefficient, MaxIterIndex, CL_design, V, rho, alt, Re, M)

    % Wing planform geometry 
    %                x    y     z   chord(m)    twist angle (deg) 
    AC.Wing.Geom = [x_section(1)     y_section(1)     z_section(1)     c_section(1)       twist_section(1);
                    x_section(2)     y_section(2)     z_section(2)     c_section(2)       twist_section(2);
                    x_section(3)     y_section(3)     z_section(3)     c_section(3)       twist_section(3)];

    % Wing incidence angle (degree)
    AC.Wing.inc  = 0;   

    % Airfoil coefficients input matrix
    %                    | ->     upper curve coeff.                <-|   | ->       lower curve coeff.       <-| 
    AC.Wing.Airfoils   = [CST_root(1)    CST_root(2)    CST_root(3)    CST_root(4)    CST_root(5)  CST_root(6)   CST_root(7)   CST_root(8)   CST_root(9)    CST_root(10)    CST_root(11)   CST_root(12);
                          CST_kink(1)    CST_kink(2)    CST_kink(3)    CST_kink(4)    CST_kink(5)  CST_kink(6)   CST_kink(7)   CST_kink(8)   CST_kink(9)    CST_kink(10)    CST_kink(11)   CST_kink(12);
                          CST_tip(1)    CST_tip(2)    CST_tip(3)    CST_tip(4)    CST_tip(5)  CST_tip(6)   CST_tip(7)   CST_tip(8)   CST_tip(9)    CST_tip(10)    CST_tip(11)   CST_tip(12)];

    percentage_kink = y_section(2)/y_section(3);
    AC.Wing.eta = [0; percentage_kink; 1];  % Spanwise location of the airfoil sections

    % Viscous vs inviscid
    AC.Visc  = Viscosity_coefficient;              % 0 for inviscid and 1 for viscous analysis
    AC.Aero.MaxIterIndex = MaxIterIndex;    %Maximum number of Iteration for the
                                    %convergence of viscous calculation

    % Flight Condition
    AC.Aero.V     = V;            % flight speed (m/s)
    AC.Aero.rho   = rho;         % air density  (kg/m3)
    AC.Aero.alt   = alt;             % flight altitude (m)
    AC.Aero.Re    = Re;        % reynolds number (bqased on mean aerodynamic chord)
    AC.Aero.M     = M;           % flight Mach number 

    AC.Aero.CL    = CL_design;          % lift coefficient - comment this line to run the code for given alpha%
    %AC.Aero.Alpha = 2;             % angle of attack -  comment this line to run the code for given cl 

    Res = Q3D_solver(AC);
    
    cl_distribution = Res.Wing.cl;  
    cm_distribution = Res.Wing.cm_c4; 
    Y_distribution = Res.Wing.Yst;
    chord_distribution = Res.Wing.chord;
    CLwing = Res.CLwing;
    if AC.Visc == 0
        CDwing = Res.CDiwing;  % if AC.Visc = 0, this value does not matter as it is used for EMWET
    end
    if AC.Visc == 1
        CDwing = Res.CDwing;   % if AC.Visc = 1, this value does matter as it is used for Brequet
    end
        
end 