%% Assign_lim_forbb_ang
% lim_forbb_ang ... the largest forbidden normal angle to Wulff shape defined by parameters in input structure in.intf.params_incl_dep
%       - given in the 1st quadrant
function lim_forbb_ang = Assign_lim_forbb_ang(in)
    lim_forbb_ang = 0;
    if any(in.is_locally_aniso_IE)
        if any([in.intf.params_incl_dep.Omega] >= 1) % i.e. STRONG ANISO (non-convex anisofun)
            clear lim_forbb_ang
            lim_forbb_ang = border_angles(in.intf.params_incl_dep,false);
            lim_forbb_ang = lim_forbb_ang(1,:); % 1st segment suffices
        end
    end
end % Assign_lim_forbb_ang