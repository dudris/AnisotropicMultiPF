%%calc_regularized_Wulff_normal_ang
% [th,th_mod, th_wulff] = calc_regularized_Wulff_normal_ang(phi,results_entry)
% - to see the interval of allowed angles set in the script  'plot_allowed_angs=true'
% INPUT
%   phi ... linspace(-pi,pi) ... interval of POLAR angles for Wulff plotting
%   results_entry ... results_entry.intf.params_incl_dep.(...)
% OUTPUT
%   phi        ... polar angles of sharp corner points added if present
%   th_mod ... allowed INTERFACE NORMAL angles, rotated (to plot regularized Wulff shape)
%   phi_rot ... POLAR angles rotated by offset_ang (to plot rotated Wulff shape)
% 
function [phi,th_mod, phi_rot] = calc_regularized_Wulff_normal_ang(phi,results_entry)

    offset_ang = results_entry.intf.params_incl_dep.offset_ang;

    if results_entry.intf.params_incl_dep.Omega>1
        % adds exact angles under which there should be the vortices and sorts ascendingly, if
        % they are present it does nothing
        phi = AddVorticesAngleCoord(results_entry.intf.params_incl_dep,phi); 
%         th = unique(th);
%         [th,ind_th_sort] = sort(th,'ascend');
        th_m = border_angles(results_entry.intf.params_incl_dep,false);
        th_m = rotate_to_first_inetrval(th_m,offset_ang);
%         !!! 
        numel_per_int = floor(length(phi)/results_entry.intf.params_incl_dep.nfold);
        plot_allowed_angs = false;
        [th_mod, ~]  = GetAllowedAngles(th_m,numel_per_int,plot_allowed_angs);
%         th_for_forbb = linspace(-pi,pi,15e3)';
%         cond_forb_ang = find_forbidden_angles(th_for_forbb,th_m,results_entry.intf.params_incl_dep);
%         th_mod = th_for_forbb(~cond_forb_ang);
    else
        th_mod = phi;
    end% if input.Omega>1
    
    % rotate the wulff plot 
        phi_rot = rotate_to_first_inetrval(phi,-offset_ang); % not sure why this takes the negative of offset ang
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

% %% GetAllowedAngles
% % use linspace function to obtain allowed angles between the border points.
% function [th_mod, ba_allwd]= GetAllowedAngles(border_angles,numpts,plotting)
% 
%     nfold = size(border_angles,1);
% %     assert(mod(nfold,2)==0,'GetAllowedAngles<calc_regularized_Wulff_normal_ang not optimized for use with odd nfold')
%     
%     % makes sure  1st column is always smaller than 2nd
%     ba_forb(:,1) = min(border_angles,[],2);
%     ba_forb(:,2) = max(border_angles,[],2);
%     
% %     rearrange to identify interval of ALLOWED angles
%     % the below fails when there should be jump in ALLOWED angles interval
%     % then the rearranged still marks FORBIDDEN intervals
%     ba_allwd = reshape(sort(ba_forb(:)),[2,nfold])';
%     negcheck = all(sort(ba_allwd(:,1))==sort(ba_forb(:,1))); % true when failed
% 
%     if negcheck 
%         % assume that allowed interval must have a jump
%         clear ba_allwd
%         ba_sorted = sort(ba_forb(:));
%         % largest and smallest value form the interval with jump
%         ba_allwd(1,[1,2]) = [max(ba_sorted), (min(ba_sorted)+2*pi) ];
%         %the rest to be processed as before
%         ba_temp = ba_sorted(2:(end-1));
%         ba_allwd(2:nfold,:) = reshape(ba_temp,[2,(nfold-1)])';
% %         ba_allwd = rotate_to_first_inetrval(ba_allwd,0);
%         negcheck = all(sort(ba_allwd(:,1))==sort(ba_forb(:,1)));
%         assert(~negcheck,'GetAllowedAngles in calc_regularized_Wulff_normal_ang failed')
%     end
%         
%     % assumes 1st column is always smaller than 2nd
%     int_width = diff(ba_allwd,[],2);
%     jump_segment_ind = find((int_width-max(int_width))>pi);
%     th_mod =  nan(nfold*numpts,1);
% 
%     indkk = 1:numpts;
%     for kk = 1:nfold
%         indkk = (1:numpts) + numpts*(kk-1);
% %         disp(num2str(indkk))
%         
%         if isempty(jump_segment_ind) || (kk~=jump_segment_ind)
%             th_mod(indkk,1) = linspace(ba_allwd(kk,1),ba_allwd(kk,2),numpts);
%         else
%             % segment with jump 
%             negind = 1:floor(numpts/2);
%             posind = (floor(numpts/2)+1):numpts;
%             th_mod(indkk(negind),1) = linspace(ba_allwd(kk,1),-pi,numel(negind));
%             th_mod(indkk(posind),1) = linspace(ba_allwd(kk,2),pi,numel(posind));
%         end% if
%     end% for
%     
%     th_mod = rotate_to_first_inetrval(th_mod,0);
%     th_mod = unique(th_mod); % to throuw away duplicaets and sort
%     
% % %     the below samples the forbidden regions, not the allowed ones
% % %     % assumes 1st column is always smaller than 2nd
% %     jump_segment_ind = find(diff(border_angles,[],2)==max(diff(border_angles,[],2)));
% %     th_mod =  nan(nfold*numpts,1);
% % % 
% %     for kk = 1:nfold
% %         indkk = (1:numpts) + numpts*(kk-1);
% % %         disp(num2str(indkk))
% %         
% %         if kk~=jump_segment_ind
% %             th_mod(indkk,1) = linspace(border_angles(kk,1),border_angles(kk,2),numpts);
% %         else
% %             % segment with jump 
% %             negind = 1:floor(numpts/2);
% %             posind = (floor(numpts/2)+1):numpts;
% %             th_mod(indkk(negind),1) = linspace(border_angles(kk,1),-pi,numel(negind));
% %             th_mod(indkk(posind),1) = linspace(border_angles(kk,2),pi,numel(posind));
% %         end% if
% %     end% for
% 
%     if plotting
%         figure(66)
%         x = 1:numel(th_mod); plot(x,th_mod,'o',x(isnan(th_mod)),th_mod(isnan(th_mod)),'rx'), 
%         polarplot(th_mod,ones(size(th_mod)),'.')
%         title('allowed angles intervals')
%     end
% % sum(isnan(th_mod))
% 
% end 
