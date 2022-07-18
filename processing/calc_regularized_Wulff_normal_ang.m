%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Benchmarking of different strategies to include anisotropy in a curvature-driven multi-phase-field model”, Mendeley Data, V2, doi: 10.17632/5wrv3ky9pp.2> 
% coupled to publication of the same name by Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754).
% 
%%calc_regularized_Wulff_normal_ang
% [phi,th_allwd, phi_rot] = calc_regularized_Wulff_normal_ang(phi,results_entry)
% - to see the interval of allowed angles set in the script  'plot_allowed_angs=true'
% INPUT
%   phi ... linspace(-pi,pi) ... interval of POLAR angles for Wulff plotting
%   results_entry ... results_entry.intf.params_incl_dep.(...)
% OUTPUT
%   phi        ... as phi, only polar angles of sharp corner points added if present
%   th_allwd ... allowed INTERFACE NORMAL angles, rotated (to plot regularized Wulff shape)
%   phi_rot ... POLAR angles rotated by offset_ang (to plot rotated Wulff shape)
% 
function [phi,th_allwd, phi_rot] = calc_regularized_Wulff_normal_ang(phi,results_entry)

    offset_ang = results_entry.intf.params_incl_dep.offset_ang;

    if results_entry.intf.params_incl_dep.Omega>1
        % __ adds exact angles under which there should be the vortices and sorts ascendingly, if
        % they are present it does nothing
        phi = AddVorticesAngleCoord(results_entry.intf.params_incl_dep,phi); 
        
        % __ interface normal angles 
        th_m = border_angles(results_entry.intf.params_incl_dep,false);
        % __ apply rotation by offset_ang and bring all angles to (-pi,pi) interval
        th_m = rotate_to_first_inetrval(th_m,offset_ang);
        
        % __ get vector of allowed angles on the nfold Wulff edges of the same length as phi
        numel_per_int = floor(length(phi)/results_entry.intf.params_incl_dep.nfold);
        plot_allowed_angs = false;
        [th_allwd, ~]  = GetAllowedAngles(th_m,numel_per_int,plot_allowed_angs);
    else
        th_allwd = phi;
    end% if input.Omega>1
    
    % __ rotate the wulff plot 
        phi_rot = rotate_to_first_inetrval(phi,-offset_ang); 
end

%% AddVorticesAngleCoord
function th = AddVorticesAngleCoord(input,th)
    
    for k = 1: input.nfold
        la = length(th);
        seg_width = 2*pi/input.nfold;
        th(la+1) = (k-1)*seg_width + input.offset_ang;
%         rw(la+1) = fIEijfun(lim_ang,input.soaIE)/cos(lim_ang);
    end % for
    th = rotate_to_first_inetrval(th,0);
    th = sort(th,'ascend');
end % func
