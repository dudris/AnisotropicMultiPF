%% [XYcentroid, r_contour, norm_sum_squared_diff] = GetContourAndCompareToWulff(S, th, input,plotting)
% INPUT
%   S      ... single matrix of phase field values where the contour is assumed to lie on value 0
              ... recommended to take p{1}-p{2} or reverse
%   th      ... column vector of angles in rad to which the PF contour and Wulff will be displayed
%   input   ... input structure 
%   plotting ... bool to supress or allow results plotting
% OUTPUT
%   XYcentroid             ... size(XYcentroid)=[1,2], geometric center of the inner PF (centroid)
%   r_contour               ... size(r_contour)=size(th), the radial distance of the PF contour interpolated to th
%   norm_sum_squared_diff ... scalar, normalized sum of squared differences between the contour and Wulff

function [XYcentroid, RW ,r_contour, norm_sum_squared_diff,HausdorffD] = GetContourAndCompareToWulff(S, th, input,plotting)
    
    th=th(:);
    
    % workaround to assure backward compatibility with older input structures
    if ~isfield(input.intf.params_incl_dep,'offset_ang') % if offset_ang is not one of fields in params_incl_dep
        input.intf.params_incl_dep.offset_ang = input.misori(1); % implying new input including misorientation
    end
    
    % th_mod ... allowed angles to plot Wulff; th_wulff = th - offset, th_wulff \in <-pi,pi>
    [th,th_mod, th_wulff] = calc_regularized_Wulff_normal_ang(th,input);
    
    % non-regularized anisotropy function
    soaIE = input.intf.params_incl_dep.soaIE;
    
    offset_ang_included = false;
    [fIEijfun,dfIEijfun,~] = AssignAnisotropyFunction(input.intf.params_incl_dep,offset_ang_included); % un-rotated anisofun
    
%     imagesc(S),set(gca,'dataaspectratio',[1,1,1],'ydir','normal'), colorbar
    
    [XYcentroid, XYcontour, r_contour] = GetContourAndCentroid(S,th);
    cond_r = ~isnan(r_contour);
%     plot(th,r_contour)
     %% dev fitting
%     C = contourc(S,[0,0]);
%     C(:,1) = [];
%     C = C;
%     
%     % Wulff parametrized by interface normal angle => find it in every
%     % point of the contour and you can fit the thing
%     %     inspired by https://itectec.com/matlab/matlab-how-to-calculate-normal-to-a-point-on-random-contour/
%     w = cscvn(C); % arclength parametrized spline through contour
%     wtang = fnder(w); % fnval(wtang,w.breaks) is tangent vector
%     wnorm = [0 -1;1 0]*fnval(wtang,w.breaks); % normal, outward-pointing
% %     visualize
% %         fnplt(w),
% %         hold on
% %         plot(C(1,:),C(2,:),'o')
% %         quiver(C(1,:),C(2,:),wnorm(1,:),wnorm(2,:))
% %         hold off
% %         axis equal
%     % get the interface normal angles
%     th_contour = atan2(wnorm(2,:),wnorm(1,:))';
%     
%     th_m  = border_angles(input.intf.params_incl_dep,false);
%     cond_fa = find_forbidden_angles(th_contour,th_m,input.intf.params_incl_dep);
%     cond_aa = ~cond_fa;
%     
%     % for independent fitting of the Wulff shape coordinates
%     wxfun = @(xx,r,t) xx + r*( fIEijfun(t,soaIE).*cos(t) - dfIEijfun(t,soaIE).*sin(t) );
%     wyfun = @(yy,r,t) yy + r*( fIEijfun(t,soaIE).*sin(t) + dfIEijfun(t,soaIE).*cos(t) );
%     % fitting
%     startfitR = (max(C(1,:))-min(C(1,:)))/(1-soaIE);
%     fo = fitoptions('method','NonlinearLeastSquares','Robust','off','StartPoint', [input.Nx/2 startfitR],'Algorithm','Trust-Region');
%     fitx = fit(th_contour(cond_aa),C(1,cond_aa)',fittype(wxfun,'independent','t'),fo);
%     fo.StartPoint = [input.Ny/2 startfitR];
%     fity = fit(th_contour(cond_aa),C(2,cond_aa)',fittype(wyfun,'independent','t'),fo);
%     % coefficients [xcenter, radius1 ; ycenter radius2]
%     cv(1,:) = coeffvalues(fitx);
%     cv(2,:) = coeffvalues(fity);
%     % compound function to get single value of radius from the data
%     compfun = @(r,t) wxfun(cv(1,1),r,t)-wyfun(cv(2,1),r,t);
%     fo = fitoptions('method','NonlinearLeastSquares','Robust','off','StartPoint', startfitR);
%     fitr = fit(th_contour(cond_aa),C(1,cond_aa)'-C(2,cond_aa)',fittype(compfun,'independent','t'),'StartPoint',startfitR);
%     fit_RW = coeffvalues(fitr);
%     
%     
%     plotc = XYcentroid;
% %     plotc = cv(:,1)';
%     
% % %     plot(C(1,:),C(2,:),'o',wxfun(plotc(1),fit_RW,th_contour(cond_aa)),wyfun(plotc(2),fit_RW,th_contour(cond_aa)),'k-')
% % %     axis equal
% % %     
% % %     ssqd(1,1) = sum((C(1,cond_aa)'-wxfun(plotc(1),fit_RW,th_contour(cond_aa))).^2);
% % %     ssqd(1,2) = sum((C(2,cond_aa)'-wyfun(plotc(2),fit_RW,th_contour(cond_aa))).^2);
% % %     
% % %     subplot(211)
% % %     plot(th_contour,C(1,:),'o',th_contour(cond_aa),wxfun(plotc(1),cv(1,2),th_contour(cond_aa)),'-')
% % %     subplot(212)
% % %     plot(th_contour(cond_aa),wxfun(plotc(1),cv(1,2),th_contour(cond_aa))-C(1,cond_aa)','x')
% % %     
% % %     subplot(211)
% % %     plot(th_contour,C(2,:),'o',th_contour(cond_aa),wyfun(plotc(2),cv(2,2),th_contour(cond_aa)),'-')
% % %     subplot(212)
% % %     plot(th_contour(cond_aa),wyfun(plotc(2),cv(2,2),th_contour(cond_aa))-C(2,cond_aa)','x')
% % %     plot(th_contour,C(1,:)-C(2,:),'o',th_contour(cond_aa),compfun(fit_RW,th_contour(cond_aa)),'x')
%     
%     RW = fit_RW;
%     
%     r_contour = CalcRadiusFromXYcontour(C'-plotc,th);
%     XYcontour = [r_contour.*cos(th) r_contour.*sin(th)];
    
    %%
%     
    % minimal radius
    if soaIE<1e-6 || std(r_contour)/mean(r_contour)<0.003
        RW = mean(r_contour);
    else
        % index of minimal value in radius
        rmin_ind = find(r_contour==min(r_contour));
        rmin_ind = rmin_ind(1); % if more points found, one is sufficient. Either they are near each other or they are separate minima - both ok
        % to locate interval around the minimal point of width 2pi/n/3
        th_min_halfw_ind = floor(numel(th)/input.intf.params_incl_dep.nfold/6);
        th_min_interv_ind = rmin_ind+(-th_min_halfw_ind:th_min_halfw_ind);
        th_min_interv_ind = th_min_interv_ind(th_min_interv_ind>1 & th_min_interv_ind<length(r_contour));
        parabola_min = polyfit(th(th_min_interv_ind),r_contour(th_min_interv_ind),2);
        % interpolate radius to the point -b/2a where f'(x)==0, minimum of the parabola
        rmin = interp1(th(th_min_interv_ind),r_contour(th_min_interv_ind), -parabola_min(2)/parabola_min(1)/2);
        RW = rmin/(1-soaIE);
    end
    r_csc = r_contour/RW;
% %     RW = min(r_contour(cond_r))/(1-soaIE);
    %%
%     if input.intf.params_incl_dep.Omega <=1
% %         rc_mid = 0.5*(min(r_contour(cond_r))+max(r_contour(cond_r))) ;% mean is not convenient ... mean(r_contour(cond_r));
%         rc_mid = min(r_contour(cond_r))/(1-soaIE);
%     else % if >1
% %         rc_mid = max(r_contour(cond_r)) / (fIEijfun(th_m(1,1),soaIE)/cos(th_m(1,1)));
%         rc_mid = min(r_contour(cond_r))/(1-soaIE);
%     end
    
    % xy coord of points in wulff shape
    wx = fIEijfun(th_mod,soaIE).*cos(th_mod) - dfIEijfun(th_mod,soaIE).*sin(th_mod) ;
    wy = fIEijfun(th_mod,soaIE).*sin(th_mod) + dfIEijfun(th_mod,soaIE).*cos(th_mod) ;
    
    % radius under angles th_wulff 
    rw = CalcRadiusFromXYcontour([wx,wy],th_wulff);
        cond_r = cond_r & ~isnan(rw);
    
%         th_mod = th_mod(cond_r);
    wx = rw.*cos(th);
    wy = rw.*sin(th);
    cx = XYcontour(:,1);
    cy = XYcontour(:,2);
        
    norm_sum_squared_diff = sum((rw(cond_r)-r_csc(cond_r)).^2);
    HausdorffD = HausdorffDist([wx,wy],XYcontour/RW);
    
    if plotting
        % plotting 1D
        figure(22)
            subplot(121)
                plot(th,r_contour,'.-',th,rw,'.-')
                grid on
                xlabel('theta')
                ylabel('radius')
                legend('PF contour','Wulff plot')
        
        % plotting in 2D
        figure(22)
            subplot(122)
                plot(cx,cy,'.',wx,wy,'.')
                xlabel('x')
                ylabel('y')
                grid on
                legend('PF contour','Wulff plot')
                axis equal
    end
    
end % func CompareToWulff


%% GetContourAndCentroid
% finds the zero level set => input must be difference of the 2 OPs
function [XYcentroid, XYinterp, r_interp] = GetContourAndCentroid(S,th)
% contourc always chooses one of x or y coordinates on the existing grid
% and finds the other so that (x,y) lays on contour
% check of this is sum(any(mod(c1,1)==0,2)) == length(c1)
    
    c1 = contourc(S,[0 0]);
    c1(:,1)=[];
    c1 = c1';
    
    % !!!!!
%     XYcentroid = mean(c1);
    XYcentroid = GetCentroid(S); % S should be the inner PF
%     XYcentroid = [size(S,2)/2 , size(S,2)/2];

    c1(:,2) = c1(:,2)-XYcentroid(2);
    c1(:,1) = c1(:,1)-XYcentroid(1);
    
    r_interp = CalcRadiusFromXYcontour(c1,th);
%     plot(c1(:,1),c1(:,2),'o'),axis equal
    
%     ang1 = atan2(c1(:,2),c1(:,1));
%     r1 = sqrt(c1(:,1).^2+c1(:,2).^2);
%     [ang1, i1] = unique(ang1);
%     r1 =r1(i1);
% 
%     r_interp(:,1) = interp1(ang1,r1,th,'linear');
    
    XYinterp = [r_interp.*cos(th) r_interp.*sin(th)];
end % func GetContourAndCentroid
%%
function XYcentroid = GetCentroid(S) % S should be the inner PF
    [Ny, Nx] = size(S);
    [X,Y] = meshgrid(1:Nx,1:Ny);
    sympref('HeavisideAtOrigin', 1);
    S = heaviside(S);
    grainarea = trapz(trapz(S));

    XYcentroid(1) = trapz(trapz(X.*S))/grainarea;
    XYcentroid(2) = trapz(trapz(Y.*S))/grainarea;
%     imagesc(Y.*S),axis equal, colorbar
%     XYcentroid(1) = sum(sum(X.*S))/grainarea;
%     XYcentroid(2) = sum(sum(Y.*S))/grainarea;
end % func GetCentroid

%%
function th = AddVorticesAngleCoord(input,th)
%     offset_ang_included = false;
%     [fIEijfun,~,~] = AssignAnisotropyFunction(input,offset_ang_included);
    for k = 1: input.nfold
        la = length(th);
        seg_width = 2*pi/input.nfold;
        th(la+1) = (k-1)*seg_width + input.offset_ang;
%         rw(la+1) = fIEijfun(lim_ang,input.soaIE)/cos(lim_ang);
    end % for
    th = rotate_to_first_inetrval(th,0);
end % func