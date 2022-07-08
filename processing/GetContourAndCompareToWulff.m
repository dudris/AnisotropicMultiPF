%% [XYcentroid, r_contour, norm_sum_squared_diff,HausdorffD] = GetContourAndCompareToWulff(S, th, input,plotting)
% INPUT
%   S      ... single matrix of phase field values where the contour is assumed to lie on value 0
%              ... recommended to take p{1}-p{2} or reverse
%   th      ... column vector of angles in rad to which the PF contour and Wulff will be displayed
%   input   ... input structure 
%   plotting ... bool to supress or allow results plotting
% OUTPUT
%   XYcentroid             ... size(XYcentroid)=[1,2], geometric center of the inner PF (centroid)
%   RW                          ... scalar radius of Wulff shape
%   r_contour               ... size(r_contour)=size(th), the radial coordinate of the PF contour interpolated to th
%   norm_sum_squared_diff ... scalar, mean of squared differences between the scaled contour and unit-radius Wulff shape
%   HausdorffD             ... Hausdorff distance computed for the scaled contour and unit-radius Wulff shape

function [XYcentroid, RW ,r_contour, norm_sum_squared_diff,HausdorffD] = GetContourAndCompareToWulff(S, th, input,plotting)
    
    th=th(:);
    
    % __ workaround to assure backward compatibility with older input structures
    if ~isfield(input.intf.params_incl_dep,'offset_ang') % if offset_ang is not one of fields in params_incl_dep
        input.intf.params_incl_dep.offset_ang = input.misori(1); % implying new input including misorientation
    end
    
    % __ th_mod ... allowed angles to plot Wulff; th_wulff = th - offset, th_wulff \in <-pi,pi>
    [th,th_mod, th_wulff] = calc_regularized_Wulff_normal_ang(th,input);
    
    % __ non-regularized anisotropy function
    soaIE = input.intf.params_incl_dep.soaIE;
    
    offset_ang_included = false;
    [fIEijfun,dfIEijfun,~] = AssignAnisotropyFunction(input.intf.params_incl_dep,offset_ang_included); % un-rotated anisofun
    
    % __ finds centroid, interpolates the radial coordinate of contour points to polar angles th
    [XYcentroid, XYcontour, r_contour] = GetContourAndCentroid(S,th);
    cond_r = ~isnan(r_contour);    
%     
    % __ determination of scalar Wulff radius
    if soaIE<1e-6 || std(r_contour)/mean(r_contour)<0.003
        RW = mean(r_contour);
    else
        % __ index of minimal value in radius
        rmin_ind = find(r_contour==min(r_contour));
        rmin_ind = rmin_ind(1); % if more points found, one is sufficient. Either they are near each other or they are separate minima - both ok
        % __ to locate interval around the minimal point of width 2pi/n/3
        th_min_halfw_ind = floor(numel(th)/input.intf.params_incl_dep.nfold/6);
        th_min_interv_ind = rmin_ind+(-th_min_halfw_ind:th_min_halfw_ind);
        th_min_interv_ind = th_min_interv_ind(th_min_interv_ind>1 & th_min_interv_ind<length(r_contour));
        parabola_min = polyfit(th(th_min_interv_ind),r_contour(th_min_interv_ind),2);
        % __ interpolate radius to the point -b/2a where f'(x)==0, minimum of the parabola
        rmin = interp1(th(th_min_interv_ind),r_contour(th_min_interv_ind), -parabola_min(2)/parabola_min(1)/2);
        RW = rmin/(1-soaIE);
    end
    r_csc = r_contour/RW;
    
    % __ xy coord of points in wulff shape
    wx = fIEijfun(th_mod,soaIE).*cos(th_mod) - dfIEijfun(th_mod,soaIE).*sin(th_mod) ;
    wy = fIEijfun(th_mod,soaIE).*sin(th_mod) + dfIEijfun(th_mod,soaIE).*cos(th_mod) ;
    
    % __ radius under angles th_wulff 
    rw = CalcRadiusFromXYcontour([wx,wy],th_wulff);
        cond_r = cond_r & ~isnan(rw);
    
    wx = rw.*cos(th);
    wy = rw.*sin(th);
    cx = XYcontour(:,1);
    cy = XYcontour(:,2);
        
    norm_sum_squared_diff = sum((rw(cond_r)-r_csc(cond_r)).^2)/sum(cond_r);
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
    
    XYcentroid = GetCentroid(S); % S should be the inner PF

    c1(:,2) = c1(:,2)-XYcentroid(2);
    c1(:,1) = c1(:,1)-XYcentroid(1);
    
    r_interp = CalcRadiusFromXYcontour(c1,th);
    
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
end % func GetCentroid
