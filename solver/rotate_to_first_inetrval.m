%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Benchmarking of different strategies to include anisotropy in a curvature-driven multi-phase-field model”, Mendeley Data, V2, doi: 10.17632/5wrv3ky9pp.2> 
% coupled to publication of the same name by Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754).
% 
%% alph_rot = rotate_to_first_inetrval(alph,offset_ang)
% bring angles outside the interval (-pi,pi) to that interval
function alph_rot = rotate_to_first_inetrval(alph,offset_ang)
    alph_rot = alph + offset_ang;
    
    in_first_interval = (alph_rot >= -pi) & (alph_rot <= pi);
    all_in_first_interval = all(in_first_interval);
    while ~all_in_first_interval
        % bringing angles around the jump to the (-pi,pi) interval
        alph_rot(alph_rot > pi) = alph_rot(alph_rot > pi) - 2*pi;
        alph_rot(alph_rot < -pi) = alph_rot(alph_rot < -pi) + 2*pi;
        % test for this is in comment behind the function
        in_first_interval = (alph_rot >= -pi) & (alph_rot <= pi);
        all_in_first_interval = all(in_first_interval);
    end

end % func


%% conceptual TEST for bringing the values out of (-pi,pi) back
% th = 170:179;
%     polarplot(th*pi/180,ones(size(th)),'.'), hold on
% offset_ang = -30;
% th = th - offset_ang;
%     polarplot(th*pi/180,ones(size(th)),'.')
%     title(['offset ang = ' num2str(offset_ang)])
% %
% th(th > 180) = th(th > 180) - 360;
%     polarplot(th*pi/180,ones(size(th)),'o')
%     legend('given angle','angle in 0deg anisofun','the same angle in (-pi,pi)')
%     hold off
%     assert(all(th>=-180 & th<=180))
% 
% th = -179:-170;
%     polarplot(th*pi/180,ones(size(th)),'.'), hold on
% offset_ang = 30;
% th = th - offset_ang;
%     polarplot(th*pi/180,ones(size(th)),'.')
%     title(['offset ang = ' num2str(offset_ang)])
% % 
% th(th < -180) = th(th < -180) + 360;
%     polarplot(th*pi/180,ones(size(th)),'o')
%     legend('given angle','angle in 0deg anisofun','the same angle in (-pi,pi)')
%     hold off
%     assert(all(th>=-180 & th<=180))