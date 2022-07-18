%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Benchmarking of different strategies to include anisotropy in a curvature-driven multi-phase-field model”, Mendeley Data, V2, doi: 10.17632/5wrv3ky9pp.2> 
% coupled to publication of the same name by Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754).
% 
%% function th_m = border_angles(params_incl_dep,plotting)
%   - returns borders of interval of forbidden angles
%   - params_incl_dep is a field in input structure 'in': in.intf.params_incl_dep
%   - for odd nfold th_m(:,1)<th_m(:,2)
function th_m = border_angles(params_incl_dep,plotting)
    
    th_m = border_angles_segm1(params_incl_dep,plotting);
    
    cond_single_maxmin = (length(th_m) == 1);
    if cond_single_maxmin
        if th_m<1e-15
%             warning('The anisotropy is probably weak (Omega<=1) or too strong [Omega>=(n^2-1)] (a single extreme at th ==0 was found in cos(th)/F(th) in interval (-pi,pi)/n )')
            th_m = 0;
%             th_m = [];
%             disp('th_m is set empty')
        else
            warning(['A single extreme found at unexpected position th == ' num2str(th_m) ' in cos(th)/F(th) in interval (-pi,pi)/n. Must check.'])
        end % if extreme is around th == 0
        return
    else
        nfold = params_incl_dep.nfold;
%         nfold = nfold(1);
        is_nfold_odd = mod(nfold ,2)~=0;
        segment_width = 2*pi/nfold;
        if is_nfold_odd
            segments_in_halfspace = floor(nfold/2);
            for k = 1:segments_in_halfspace
                th_m(k+1,:) = th_m(1,:) + k*segment_width;
                th_m(segments_in_halfspace+k+1,:) = th_m(1,:) - k*segment_width;
            end
            th_m = sort(th_m,2,'ascend');
        else % nfold is even
            segments_in_halfspace = nfold/2 -1;
            for k = 1:segments_in_halfspace
                th_m(k+1,:) = th_m(1,:) + k*segment_width;
                th_m(segments_in_halfspace+k+1,:) = th_m(1,:) - k*segment_width;
            end
            th_m_neg = th_m(1,th_m(1,:)<0);
            th_m_pos = th_m(1,th_m(1,:)>0);
            % to the last rotation I must add length of the ALLOWED interval 
            allowed_segment_width = 2*(pi/nfold-th_m_pos);
            th_m(nfold,1) = th_m_neg - segment_width*segments_in_halfspace - allowed_segment_width;
            th_m(nfold,2) = th_m_pos + segment_width*segments_in_halfspace + allowed_segment_width;
    %         th_m(nfold,:) = th_m(1,:) + nfold/2*segment_width;
        end % if is_nfold_odd
    end % if cond_single_maxmin
    
%     if input.offset_ang ~= 0
%         % rotate to get rid of offset
%         th_m = th_m + offset_ang;
%         % bringing angles around the jump to the (-pi,pi) interval
%         th(th > pi) = th(th > pi) - 2*pi;
%         th(th < -pi) = th(th < -pi) + 2*pi;
%         % test for this is in comment behind the function
%     end
end % func main





%% border_angles_segm1
% assumes no offset angle
function th_m = border_angles_segm1(params_incl_dep,plotting)
    
    nfold = [params_incl_dep.nfold];
    nfold = nfold(1);
    delta = [params_incl_dep.soaIE];
    delta = delta(1);
    assert(all(([params_incl_dep.nfold]-nfold)==0),'too many inclination dpendent interfaces specified. Check also border_angles>border_angles_segm1.')
    assert(all(([params_incl_dep.soaIE]-delta)==0),'too many inclination dpendent interfaces specified. Check also border_angles>border_angles_segm1.')
    
    syms th tangline(th) dtangline(th)
    tangline = cos(th) / (1+delta*cos(nfold*th) );
    dtangline = diff(tangline);

    g= matlabFunction(tangline);
    dg= matlabFunction(dtangline);

    % fitting of dg by polynomial of polyorder-th order 
    ang_all = linspace(-pi/nfold,pi/nfold,300)';
    
    if plotting 
        figure(22)
            plot(ang_all*180/pi,g(ang_all),'.',ang_all*180/pi,dg(ang_all),'o')
            grid on
            hold on
%             plot(ang*180/pi,polymatrix*coeff,'LineWidth',1.5)
%             plot([th_m;th_m]*180/pi,[g(th_m);dg(th_m)],'kd')
%             legend('cos(th)/f(th)','d/d th ( cos(th)/f(th) )','fit poly 9','maxima','Location','south')
            xlim(180/nfold*[-1,1])
            xticks([-180/nfold,0,180/nfold])
            xticklabels({'-\pi/n','0','\pi/n'})
            hold off
            set(gca,'FontSize',14)
    end
    
    ang = find_regions_of_changing_sign(ang_all,dg); % ang is a cell, size(ang)=[2,1]
    
    for k = 1:length(ang)
        [coeff , polymatrix] = fit_poly_to_data(5,ang{k},dg); %fit_poly_to_data(polyorder,ang,dg)

        % finding roots of the polynomial 
        tth_m = (roots(coeff));
        tth_m = tth_m(imag(tth_m)==0);
        is_in_neighbourhood = ( abs(tth_m) < max(abs(ang{k})) ) & ( abs(tth_m) > min(abs(ang{k})) );
        tth_m = tth_m(is_in_neighbourhood);
        
        num_of_real_roots_in_segment = length(tth_m);
        % sometimes the roots lay in the absolute value of the interval but with other sign
        if num_of_real_roots_in_segment == 2 && sign(tth_m(1))~=sign(tth_m(2))
            wtf1 = mean(ang{k}-tth_m(1)); % mean of the value which lies in the interval must be closer to 0 than the other
            wtf2 = mean(ang{k}-tth_m(2));
            if abs(wtf1) > abs(wtf2)
                tth_m = tth_m(2);
            elseif abs(wtf1) < abs(wtf2)
                tth_m = tth_m(1);
            else
                disp('wtf')
            end
            num_of_real_roots_in_segment = 1;
        end

        if plotting
            figure(22)
            hold on
            plot(ang{k}*180/pi,polymatrix*coeff,'--k','LineWidth',1.5)
            hold off
        end % if plotting
        assert(num_of_real_roots_in_segment == 1,['1 root around maximum expected, observed: ' num2str(num_of_real_roots_in_segment) ])
        th_m(k) = tth_m;
    end% fork = 1:size(ang,2)
    
    if plotting
        figure(22)
        hold on
        plot([th_m;th_m]*180/pi,[g(th_m);dg(th_m)],'kd')
        if length(ang)==1
            legend('cos(th)/f(th)','d/d th ( cos(th)/f(th) )','fit poly 5','maxima','Location','south')
        else
            legend('cos(th)/f(th)','d/d th ( cos(th)/f(th) )','fit poly 5','fit poly 5','maxima','Location','south')
        end
        hold off
    end
    
    
end % func border_angles_segm1

%% fit_poly_to_data
function [coeff , polymatrix] = fit_poly_to_data(polyorder,ang,dg)
    numofpts = length(ang);
    polymatrix = ones(numofpts,polyorder+1); 
    for k = 1:(polyorder+1)
        polymatrix(:,k) = ang.^(10-k);
    end
    coeff = polymatrix\dg(ang); % coefficients [a; b; c; ... ; j] in ax^8 + bx^7 + cx^6+.... + j
end % func fit_poly_to_data

%% find_regions_of_changing_sign
function [ang] = find_regions_of_changing_sign(ang_all,dg)
    sgndg = sign(dg(ang_all));
    ind = find(abs(diff(sgndg))==2); % where sign changes the diff signum is equal to -2 or 2
    
    numof_crossing_pts = length(ind);
    assert(numof_crossing_pts==3 | numof_crossing_pts==1,['number of crossing points is not 1 nor 3, numof_crossing_pts=' num2str(numof_crossing_pts)])
    
    % find indices of the points in the area of crossing points
    width_in_pts = 30;
    num_of_all_pts = length(ang_all);
    
    reglim_low = ind-(width_in_pts/2);
    if any(reglim_low<1)
        reglim_low(reglim_low<1) = 1;
    end
    reglim_high = ind+(width_in_pts/2-1);
    if any(reglim_high>num_of_all_pts)
        reglim_high(reglim_high>num_of_all_pts) = num_of_all_pts;
    end
    
%     assert(mean(ang_all(reglim_low(2):reglim_high(2)))==0,'the middle crossing point not second of three')
    
    if numof_crossing_pts==3
        reg1 = ang_all(reglim_low(1):reglim_high(1));
        reg2 = ang_all(reglim_low(3):reglim_high(3));
        ang = {reg1; reg2};
    elseif numof_crossing_pts==1
        ang = {ang_all(reglim_low(1):reglim_high(1))};
    end
    
end % func find_regions_of_changing_sign

