%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Benchmarking of different strategies to include anisotropy in a curvature-driven multi-phase-field model”, Mendeley Data, V2, doi: 10.17632/5wrv3ky9pp.2> 
% coupled to publication of the same name by Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754).
% 
%% Assign_lim_forbb_ang
% lim_forbb_ang = Assign_lim_forbb_ang(in)
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