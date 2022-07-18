%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Benchmarking of different strategies to include anisotropy in a curvature-driven multi-phase-field model”, Mendeley Data, V2, doi: 10.17632/5wrv3ky9pp.2> 
% coupled to publication of the same name by Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754).
% 
%% r_interp = CalcRadiusFromXYcontour(XYcontour,th)
% - Transfers contour matrix from Cartesian to polar coordinates and interpolates the radius
% in values th. Assumes the contour is already centered.
% - resulting contour to be plotted in polar coordinates (r_interp,th)
% - Occurrence of NaNs from the interpolation at the angle jump is prevented
% by extension of both curve ends by points from the other end.
% INPUT
%   XYcontour - contour matrix with n points, size(XYcontour)=[n,2]
%       - contour to be a closed curve
%   th - polar angles (rad) for the 
% OUTPUT
% r_interp - radius interpolated in angles th

function r_interp = CalcRadiusFromXYcontour(XYcontour,th)
% __ XYcontour ... size(XYcontour) = [ptscount 2] = [ xcoord ycoord ]
    
    ang2 = atan2(XYcontour(:,2),XYcontour(:,1));

    [~ , indsorted] = sort(ang2,'ascend');
    % __ jump in angles removed from middle
    XYcontour = XYcontour(indsorted,:);
    % __add some more poits from back to the front and vice versa to avoid NaNs in interpolation
    addpts = 10;
    numpts = size(XYcontour,1);
    ind_front_to_back = 1:addpts;
    ind_back_to_front = (numpts-addpts+1):numpts;
    ind_front = ind_front_to_back;
    ind_back = (numpts+addpts+1):(numpts+2*addpts);
    % __ extended contour points
    XYcontour_ext = zeros(numpts+2*addpts,2);
    XYcontour_ext(ind_front,:) = XYcontour(ind_back_to_front,:);
    XYcontour_ext(ind_back,:) = XYcontour(ind_front_to_back,:);
    XYcontour_ext(addpts+(1:numpts),:) = XYcontour;

    ang1 = atan2(XYcontour_ext(:,2),XYcontour_ext(:,1));
    r1 = sqrt(XYcontour_ext(:,1).^2+XYcontour_ext(:,2).^2);
    
    % __ angles at front to be negative thanks to sorting before
    ang1(ind_front)= ang1(ind_front)-2*pi;
    % __ angles at end to be positive thanks to sorting before
    ang1(ind_back)= ang1(ind_back)+2*pi;
    
    [ang1, i1] = unique(ang1);
    r1 =r1(i1);
    
    r_interp(:,1) = interp1(ang1,r1,th,'linear');
    
end