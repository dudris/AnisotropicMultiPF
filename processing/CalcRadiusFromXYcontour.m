%% r_interp = CalcRadiusFromXYcontour(XYcontour,th)
% - Transfers contour matrix from Cartesian to polar coordinates and interpolates the radius
% in values th. Assumes the contour is already centered.
% - Occurrence of NaNs from the interpolation at the angle jump is prevented
% by extension of both curve ends by points from the other end.
% INPUT
% XYcontour - contour matrix with n points, size(XYcontour)=[n,2]
% - contour to be a closed curve
% OUTPUT
% r_interp - radius interpolated in angles th

function r_interp = CalcRadiusFromXYcontour(XYcontour,th)
% XYcontour ... size(XYcontour) = [ptscount 2] = [ xcoord ycoord ]
    
    ang2 = atan2(XYcontour(:,2),XYcontour(:,1));
%     plot(ang2,'o')
%     r2 = sqrt(XYcontour(:,1).^2+XYcontour(:,2).^2);
    [~ , indsorted] = sort(ang2,'ascend');
    % jump in angles removed from middle
    XYcontour = XYcontour(indsorted,:);
    % add some more poits from back to the front and vice versa to avoid
    % NaNs in interpolation
    addpts = 10;
    numpts = size(XYcontour,1);
    ind_front_to_back = 1:addpts;
    ind_back_to_front = (numpts-addpts+1):numpts;
    ind_front = ind_front_to_back;
    ind_back = (numpts+addpts+1):(numpts+2*addpts);
    % extended contour points
    XYcontour_ext = zeros(numpts+2*addpts,2);
    XYcontour_ext(ind_front,:) = XYcontour(ind_back_to_front,:);
    XYcontour_ext(ind_back,:) = XYcontour(ind_front_to_back,:);
    XYcontour_ext(addpts+(1:numpts),:) = XYcontour;
    
%     plot(XYcontour(:,1),XYcontour(:,2),'.'), axis equal, hold on
%     plot(XYcontour(ind_back_to_front,1),XYcontour(ind_back_to_front,2),'or')
%     plot(XYcontour(ind_front_to_back,1),XYcontour(ind_front_to_back,2),'ob')
%     plot(XYcontour_ext(ind_front,1),XYcontour_ext(ind_front,2),'xr')
%     plot(XYcontour_ext(ind_back,1),XYcontour_ext(ind_back,2),'xb')
    
    ang1 = atan2(XYcontour_ext(:,2),XYcontour_ext(:,1));
    r1 = sqrt(XYcontour_ext(:,1).^2+XYcontour_ext(:,2).^2);
    
    % angles at front to be negative thanks to sorting before
    ang1(ind_front)= ang1(ind_front)-2*pi;
    % angles at end to be positive thanks to sorting before
    ang1(ind_back)= ang1(ind_back)+2*pi;
    
%     plot(ang1,'o')
    
%     % move added points out of (-pi,pi)
%     jumpinds = find(abs(diff(ang1))>6);
%     if ~isempty(jumpinds) && length(jumpinds)==2
%         front = 1:jumpinds(1);
%         back = (jumpinds(2)+1):size(ang1,1);
%         front_to_shift = ang1(front)<0;
%         back_to_shift = ang1(back)>0;
%         ang1(front(front_to_shift)) = ang1(front(front_to_shift)) + 2*pi;
%         ang1(back(back_to_shift)) = ang1(back(back_to_shift)) - 2*pi;
%     elseif length(jumpinds)==1
%         front = 1:jumpinds;
%         back = (jumpinds+1):size(ang1,1);
%         if length(front)>length(back)
%             back_to_shift = ang1(back)>0;
%             ang1(back(back_to_shift)) = ang1(back(back_to_shift)) - 2*pi;
%         elseif length(front)<length(back)
%             front_to_shift = ang1(front)<0;
%             ang1(front(front_to_shift)) = ang1(front(front_to_shift)) + 2*pi;
%         else % length(front)==length(back)
%             warning('CalcRadiusFromXYcontour: unexpected situation in prevention of NaNs in radius interpolation. ')
%         end
%     end
    
    [ang1, i1] = unique(ang1);
    r1 =r1(i1);
    
    r_interp(:,1) = interp1(ang1,r1,th,'linear');
% plot(th,r_interp,'o')    
%     plot(ang1,r1,'o')    
%     plot(ang1,'o')    
%     plot(diff(ang1),'o')    
    
end