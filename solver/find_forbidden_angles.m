% cond_forb_ang = find_forbidden_angles(th,th_m,input)
% th_m ... limiting angles in individual segments, size(th_m) = [nfold, 2]
% input ... in.intf.params_incl_dep
function cond_forbb_ang = find_forbidden_angles(th,th_m,input)
    
    cond_forbb_ang = false(size(th));

    forb_interv_width = diff(th_m,1,2);
    segment_width = 2*pi/input.nfold;
    is_not_segment_w_jump = (forb_interv_width - segment_width)<0 ;
    
    for k = 1:size(th_m,1)
        if is_not_segment_w_jump(k)
            cond_forb_ang_segm = th > th_m(k,1)  & th < th_m(k,2);
            cond_forbb_ang = cond_forbb_ang | cond_forb_ang_segm;
            
        else % it is a segment with jump in angles
            cond_forb_ang_segm = ( th < th_m(k,1) & th >= -pi ) | ( th > th_m(k,2) & th <= pi ); 
            cond_forbb_ang = cond_forbb_ang | cond_forb_ang_segm;
        end
    end
    
end




%
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     cond_forb_ang = false(size(th));
%     nfold = input.nfold;
%     
%     is_nfold_odd = mod(nfold ,2)~=0;
% 
%     segment_width = 2*pi/nfold;
% 
%     sign = [1 -1];
%     
%     
%     
%     if is_nfold_odd
%             is_in_segment = th>-th_m & th<th_m;
%             segments_in_halfspace = floor(nfold/2);
%             for k = 1:segments_in_halfspace
%                 for halfspace = 1:2
%                     th_temp = th + sign(halfspace)*k*segment_width;
%                     is_in_segment = th_temp>lim_bott & th_temp<lim_upp;
%                     aniso_reg(is_in_segment) = added_arc_diam*cos(th_temp(is_in_segment));
%                     daniso_reg(is_in_segment) = - added_arc_diam*sin(th_temp(is_in_segment));
%                 end
%             end
% 
%         else % nfold even
%             segments_in_halfspace = nfold/2 -1;
%             % 1st segment
%             is_in_segment = th>lim_bott & th<lim_upp;
%             aniso_reg(is_in_segment) = added_arc_diam*cos(th(is_in_segment));
%             daniso_reg(is_in_segment) = - added_arc_diam*sin(th(is_in_segment));
%             % opposite segment (the one with jump in angle)
%             is_in_segment = (th+pi<lim_upp & th+pi>=0) ;
%             aniso_reg(is_in_segment) = added_arc_diam*cos(th(is_in_segment)+pi);
%             is_in_segment = (th-pi>lim_bott & th-pi<=0) ;
%             aniso_reg(is_in_segment) = added_arc_diam*cos(th(is_in_segment)-pi);
% 
%                 for k = 1:segments_in_halfspace
%                     for halfspace = 1:2
%                         th_temp = th + sign(halfspace)*k*segment_width;
%                         is_in_segment = th_temp>lim_bott & th_temp<lim_upp;
%                         aniso_reg(is_in_segment) = added_arc_diam*cos(th_temp(is_in_segment));
%                         daniso_reg(is_in_segment) = - added_arc_diam*sin(th_temp(is_in_segment));
%                     end
%                 end
%             %     figure(4)
%         %         subplot(121)
%         %         polarplot(th(is_in_segment),strong_aniso_reg(is_in_segment),'r.')
%         %         subplot(122)
%         %         polarplot(th(is_in_segment),1./strong_aniso_reg(is_in_segment),'r.')
%         %         pause
%         end % if is_nfold_odd
% 
% end