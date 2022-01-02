%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Three strategies of parameters assignment
% in multi phase field model of grain growth with anisotorpic grain boundary properties”, Mendeley Data, 
% v1 http://dx.doi.org/10.17632/5wrv3ky9pp.1>, coupled to publication of the same name by 
% Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% 
%% [aniso_reg , daniso_reg] = calc_strong_reg_anisofun(th, th_m,input,soa,plotting)

function [aniso_reg , daniso_reg] = calc_strong_reg_anisofun(th, th_m,input,soa,plotting)
    
    offset_ang = input.offset_ang;
    
    th = rotate_by_offset_ang_to_seg1(th,offset_ang);
    
    offset_ang_included = false;
    [fIEijfun,dfIEijfun,~] = AssignAnisotropyFunction(input,offset_ang_included);
    aniso = fIEijfun(th,soa);
    aniso_reg = aniso;
    daniso_reg = dfIEijfun(th,soa);
    
    Omg = input.Omega;
    
    if Omg <= 1
        disp('msg calc_strong_reg_anisofun: Anisotropy is weak (Omega <=1). Unmodified anisotropy function written to the output,')
        return
    else %  Omg > 1
%         th_m =  border_angles(input,false);
        nfold = input.nfold;
        is_nfold_odd = mod(nfold ,2)~=0;

        lim_bott = min(th_m(1,:)) ;
        lim_upp = max(th_m(1,:)) ;
        added_arc_diam = fIEijfun(lim_bott,input.soaIE)/cos(lim_bott);
        segment_width = 2*pi/nfold;

        sign = [1 -1];

        if is_nfold_odd
            is_in_segment = th>lim_bott & th<lim_upp;
            aniso_reg(is_in_segment) = added_arc_diam*cos(th(is_in_segment));
            daniso_reg(is_in_segment) = - added_arc_diam*sin(th(is_in_segment));
            segments_in_halfspace = floor(nfold/2);
            for k = 1:segments_in_halfspace
                for halfspace = 1:2
                    th_temp = th + sign(halfspace)*k*segment_width;
                    is_in_segment = th_temp>lim_bott & th_temp<lim_upp;
                    aniso_reg(is_in_segment) = added_arc_diam*cos(th_temp(is_in_segment));
                    daniso_reg(is_in_segment) = - added_arc_diam*sin(th_temp(is_in_segment));
                end
            end

        else % nfold even
            segments_in_halfspace = nfold/2 -1;
            % 1st segment
            is_in_segment = th>lim_bott & th<lim_upp;
            aniso_reg(is_in_segment) = added_arc_diam*cos(th(is_in_segment));
            daniso_reg(is_in_segment) = - added_arc_diam*sin(th(is_in_segment));
            % opposite segment (the one with jump in angle)
            is_in_segment = (th+pi<lim_upp & th+pi>=0) ;
            aniso_reg(is_in_segment) = added_arc_diam*cos(th(is_in_segment)+pi);
            daniso_reg(is_in_segment) = - added_arc_diam*sin(th(is_in_segment)+pi);
            is_in_segment = (th-pi>lim_bott & th-pi<=0) ;
            aniso_reg(is_in_segment) = added_arc_diam*cos(th(is_in_segment)-pi);
            daniso_reg(is_in_segment) = - added_arc_diam*sin(th(is_in_segment)-pi);
            % rest of the segments
            for k = 1:segments_in_halfspace
                for halfspace = 1:2
                    th_temp = th + sign(halfspace)*k*segment_width;
                    is_in_segment = th_temp>lim_bott & th_temp<lim_upp;
                    aniso_reg(is_in_segment) = added_arc_diam*cos(th_temp(is_in_segment));
                    daniso_reg(is_in_segment) = - added_arc_diam*sin(th_temp(is_in_segment));
                end
            end
            %     figure(4)
        %         subplot(121)
        %         polarplot(th(is_in_segment),strong_aniso_reg(is_in_segment),'r.')
        %         subplot(122)
        %         polarplot(th(is_in_segment),1./strong_aniso_reg(is_in_segment),'r.')
        %         pause
        end % if is_nfold_odd
    
    end % if Omg <= 1

    if plotting
        figure(4)
            subplot(121)
                polarplot(th,aniso,'.',th,aniso_reg,'.')
                legend('unmodified','regularized')
                title('F(\theta) plot')
            subplot(122)
                polarplot(th,1./aniso,'.',th,1./aniso_reg,'.')
                legend('unmodified','regularized')
                title('1/F(\theta) plot')
        disp('msg replace_concave_segments: plotting to figure 4.')
    end % if plotting

end% function


%% fnction
function th = rotate_by_offset_ang_to_seg1(th,offset_ang)
    % rotate to get rid of offset, offset = 0 does not do anything
    th = th - offset_ang;
    % bringing angles around the jump to the (-pi,pi) interval
    th(th > pi) = th(th > pi) - 2*pi;
    th(th < -pi) = th(th < -pi) + 2*pi;
    % test for this is in comment behind the function
        
end % function rotate_by_offset_ang_to_seg1

