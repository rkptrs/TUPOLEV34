function [V_tank] = fuel_volume(Xtu_root, Xtu_kink, Xtu_fte, Xtl_root, Xtl_kink, Xtl_fte, y_fte, c_fte, y_section, c_section)
    %calculate root area
    root_indices_low = find(Xtu_root(:,2) > 0.15*c_section(1), 1, 'first');
    root_indices_high = find(Xtu_root(:,2) < 0.8*c_section(1), 1, 'last');
    A_eff_root = 0;
    for i = root_indices_low: root_indices_high-1
        A_eff_root_add = (Xtu_root(i+1,2) - Xtu_root(i,2))*((Xtu_root(i+1,1)-Xtl_root(i+1,1)) + (Xtu_root(i,1)-Xtl_root(i,1)))/2;
        A_eff_root = A_eff_root + A_eff_root_add;
    end
    
    %calculate kink area
    kink_indices_low = find(Xtu_kink(:,2) > 0.25*c_section(2), 1, 'first');
    kink_indices_high = find(Xtu_kink(:,2) < 0.65*c_section(2), 1, 'last');
    A_eff_kink = 0;
    for i = kink_indices_low: kink_indices_high-1
        A_eff_kink_add = (Xtu_kink(i+1,2) - Xtu_kink(i,2))*((Xtu_kink(i+1,1)-Xtl_kink(i+1,1)) + (Xtu_kink(i,1)-Xtl_kink(i,1)))/2;
        A_eff_kink = A_eff_kink + A_eff_kink_add;
    end

    %calculate fuel tank end (fte) area
    fte_indices_low = find(Xtu_fte(:,2) > 0.2*c_fte, 1, 'first');
    fte_indices_high = find(Xtu_fte(:,2) < 0.6*c_fte, 1, 'last');
    A_eff_fte = 0;
    for i = fte_indices_low: fte_indices_high-1
        A_eff_fte_add = (Xtu_fte(i+1,2) - Xtu_fte(i,2))*((Xtu_fte(i+1,1)-Xtl_fte(i+1,1)) + (Xtu_fte(i,1)-Xtl_fte(i,1)))/2;
        A_eff_fte = A_eff_fte + A_eff_fte_add;
    end
    
    %calculate total fuel tank volume
    V_tank = 2*((y_section(2)*(A_eff_root + A_eff_kink)/2) + (y_fte - y_section(2))*(A_eff_kink + A_eff_fte)/2);
end