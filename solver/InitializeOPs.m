%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Three strategies of parameters assignment
% in multi phase field model of grain growth with anisotorpic grain boundary properties”, Mendeley Data, 
% v1 http://dx.doi.org/10.17632/5wrv3ky9pp.1>, coupled to publication of the same name by 
% Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754)
% 
%% InitializeOPs
% produces several simple geometric initial conditions for phase field simulations.
% INPUT
%   - in ... structure; used fields: ICcode,ICparam,Nx,Ny,nOP
%       - in.ICcode ... string, either of: 'CircleInMatrix', 'Wulff_weak', 'Tjunction' or '2CirclesInMatrix'
%       - in.ICparam ... vector defining the geometry, specific to ICcode
%       - in.Nx, in.Ny ... grid dimensions in x and y direction. Note that
%       in matrix notation the point (x,y) is in the matrix element [y,x]
%       - in.nOP ... number of phase fields. Assigned automatically from ICcode in 'make_input.m'
% OUTPUT
%   - p ... cell, size(p) = [in.nOP,1], every cell is a matrix of size [in.Ny,in.Nx]
% 
function p = InitializeOPs(in,showIC)
    Nx = in.Nx;
    Ny = in.Ny;
    p = cell(in.nOP,1);
    roundlim = eps;
    param = in.ICparam;
    for i = 1:in.nOP
        p{i} = roundlim*ones([Ny,Nx]); % 
    end
    [X,Y] = meshgrid((1:Nx),(1:Ny));
    
    switch in.nOP
        case 2
            switch in.ICcode
                case 'CircleInMatrix'
                    % circular grain in matrix param = [centerx centery radius] 
                    assert(length(param)==3,'wrong number of InitializeOPs defining parameters')
                    center = [ param(1), param(2)];
                    radius = param(3);
                    condINIT= sqrt((Y-center(2)).^2+ (X-center(1)).^2) <= radius;

                case 'Wulff'
                    assert(length(param)==4,'wrong number of InitializeOPs defining parameters')
                    center = [ param(1), param(2)];
                    radius_W = param(3);
                    Omega = param(4);
                    soaIE = Omega/(in.intf.params_incl_dep.nfold^2-1);
                    if ~isfield(in.intf.params_incl_dep,'offset_ang') % if offset_ang is not one of fields in params_incl_dep
                        in.intf.params_incl_dep.offset_ang = in.misori(1); % implying new input including misorientation
                    end
                    ba_ncnvx = border_angles(in.intf.params_incl_dep,false);
%                         ba_ncnvx = border_angles_v2(in.intf.params_incl_dep,'forbb',false);
                    ba_ncnvx = rotate_to_first_inetrval(ba_ncnvx,in.intf.params_incl_dep.offset_ang);
                    [th_allwd_ncnvx, ~] = GetAllowedAngles(ba_ncnvx,100,false);
                    th_mod = th_allwd_ncnvx;
                    offset_ang_included = true;
                    [fIEijfun,dfIEijfun,~] = AssignAnisotropyFunction(in.intf.params_incl_dep,offset_ang_included);
                    wx =  radius_W * ( fIEijfun(th_mod,soaIE).*cos(th_mod) - dfIEijfun(th_mod,soaIE).*sin(th_mod) );
                    wy =  radius_W * (fIEijfun(th_mod,soaIE).*sin(th_mod) + dfIEijfun(th_mod,soaIE).*cos(th_mod) );
%                         figure(33),plot(wx,wy,'o'),axis equal

                    [TH,R] = cart2pol(X-center(1),Y-center(2));
                    TH(isnan(TH)) = pi;
%                         TH_W = rotate_to_first_inetrval(reshape(TH,[numel(X),1]),-in.intf.params_incl_dep.offset_ang);
                    TH_W = rotate_to_first_inetrval(reshape(TH,[numel(X),1]),0);
                    RW = CalcRadiusFromXYcontour([wx,wy],TH_W);
                    RW = reshape(RW,[Ny,Nx]);
                    condINIT = R<=RW;
            end % switch ICcode, nOP =2 

        p{2}(condINIT) = 1-roundlim;% initial condition
        p{1}(~condINIT) = 1-roundlim;% initial condition

        case 3 % switch in.nOP==3
            switch in.ICcode
                case 'Tjunction'
                    % vertical plane ... p{1} liquid, RIGHT ; p{2},p{3} LEFT, param = [posSL , posSS]
                    assert(length(param)==2,'wrong number of InitializeOPs defining parameters')
                    posSL = param(1);
                    posSS = param(2);
                    condINIT_S1=X <= posSS;
                    condINIT_L=Y >= posSL;
                    p{1}(condINIT_L) = 1- roundlim;
                    p{2}(~condINIT_L&condINIT_S1) = 1- roundlim;
                    p{3}(~condINIT_L&~condINIT_S1) = 1- roundlim;

                case '2CirclesInMatrix'
                    % 2 circles in matrix, param = [center1x center1y radius1 ; center2x center2y radius2]
                    assert(numel(param)==6,'wrong number of InitializeOPs defining parameters')
                    centers = [param(1,1) param(1,2) ; param(2,1) param(2,2) ];
                    radius = [param(1,3) , param(2,3)];
                    condINIT_1 = sqrt((Y-centers(1,2)).^2+ (X-centers(1,1)).^2) <= radius(1);
                    condINIT_2 = sqrt((Y-centers(2,2)).^2+ (X-centers(2,1)).^2) <= radius(2);
                    p{1}(~condINIT_1 & ~condINIT_2) = 1- roundlim;
                    p{2}(condINIT_1) = 1- roundlim;
                    p{3}(condINIT_2) = 1- roundlim;
            end
    end % switch nOP
    
    disp([num2str(in.nOP) ' phase fields initialized: ''' in.ICcode ''' with parameters [' num2str(param(:)') ']'])
    
    if showIC
        showOPs = zeros(size(p{1}));
        for k = 1:in.nOP
            showOPs = showOPs + p{k}*k;
        end
    figIS = figure(1);
    imagesc(showOPs),daspect([1 1 1]),
    cb = colorbar;
    cb.Ticks = 1:in.nOP;
    cb.Ticks = 1:in.nOP;
    cb.Limits = [0.9 in.nOP+0.1];
    set(gca,'Ydir','normal')
    title(['\Sigma_{i=1}^{nOP} i*p_i'])
%     pause
%     close(figIS)
    end % if show Initial condition

end % func InitializeOPs

%% CalcRadiusFromXYcontour
function r_interp = CalcRadiusFromXYcontour(XYcontour,th)
% XYcontour ... size(XYcontour) = [ptscount 2] = [ xcoord ycoord ]
    ang1 = atan2(XYcontour(:,2),XYcontour(:,1));
    r1 = sqrt(XYcontour(:,1).^2+XYcontour(:,2).^2);
    [ang1, i1] = unique(ang1);
    r1 =r1(i1);

    r_interp(:,1) = interp1(ang1,r1,th,'linear');
end
