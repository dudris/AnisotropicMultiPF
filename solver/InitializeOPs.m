%% InitializeOPs
% function p = InitializeOPs(in,showIS)
% from 'in' structure it uses fields: ICcode,ICparam,Nx,Ny,nOP
function p = InitializeOPs(in,showIS)
    Nx = in.Nx;
    Ny = in.Ny;
    p = cell(in.nOP,1);
    roundlim = eps;
    param = in.ICparam;
    for i = 1:in.nOP
        p{i} = roundlim*ones([Ny,Nx]); % 
    end
    [X,Y] = meshgrid((1:Nx),(1:Ny));
    
    if ~isfield(in,'iterno') % field counting energ. assessment iterations does not exist
        in.iterno = 1;
    end
    
    % in case of iterative energy assessment
    if in.iterno > 1
        assert(in.nOP == 2,['in.iterno==' num2str(in.iterno) ' and in.nOP~=2. '])
        C = in.ICcontour{in.iterno};
        % calculate 'scale' such that I am at least 1.5IW far from boundary
        IWboundaryoffset = 2;
        scale(1) = ( in.Nx - 2*IWboundaryoffset*(in.IW/in.dx) )/(max(C(:,1))-min(C(:,1)));
        scale(2) = ( in.Ny - IWboundaryoffset*(in.IW/in.dx) )/(max(C(:,2))-min(C(:,2)));
        scale = min(scale); 
%         scale = min([scale,1.95]); % maximal scale is 2
        [condINIT,dimsscaled] = interp_yval_of_wetting_contour(C,[Nx,Ny],scale, true);
        p{1}(~condINIT) = 1-roundlim; 
        p{2}(condINIT) = 1-roundlim;
        
    else % in.iterno == 1
        switch in.nOP
            case 2
                switch in.ICcode
                    case 'CircleInMatrix'
                        % circular grain in matrix param = [centerx centery radius] 
                        assert(length(param)==3,'wrong number of InitializeOPs defining parameters')
                        center = [ param(1), param(2)];
                        radius = param(3);
                        condINIT= sqrt((Y-center(2)).^2+ (X-center(1)).^2) <= radius;
                        
                    case 'EAVarAng' % Energy Assessment with Variable Angle                        
                        % compound of 2 circular arcs with different contact angle with x axis connected by tangent  line
                        % param= [angLdeg angRdeg init_rad_meters ]
                        assert(length(param)==3,'wrong number of InitializeOPs defining parameters')
                        
                        fitfactor = 0.7; % ratio of maximal shape width (and height) to domain width (and height)
                        
                        condINIT = JoinVarangCircSegments(in,param,fitfactor,false);
                        
                    case 'VerticalPlane'
                        % vertical plane ... p{1} LEFT ; P{2} RIGHT, param = [pos]
                        assert(length(param)==1,'wrong number of InitializeOPs defining parameters')
                        pos = param(1);
                        condINIT=X >= pos;
                    case 'TiltedPlane'
                        % straight line intersecting y=1 and y=Ny in points x1 and x2 
                        assert(length(param)==2,'wrong number of InitializeOPs defining parameters')
                        x1 = param(1); 
                        x2 = param(2);
                        ab = [x1 1; x2 1]\[1;Ny];
                        % to make sure OP1 is left all the time
                        if ab(1) <= 0
                            condINIT = Y > ab(1)*X + ab(2);
                        elseif ab(1) > 0
                            condINIT = Y < ab(1)*X + ab(2);
                        end
                    case 'TiltedPlaneAng'
                        % straight line intersecting y=1 and y=Ny in points x1 and x2 
                        assert(length(param)==1,'wrong number of InitializeOPs defining parameters')
                        ang = param(1); 
                        slope = tand(ang);
                        line = slope*(X-Nx/2) + Ny/2;
                        % to make sure OP1 is left all the time
                        if slope <= 0
                            condINIT = Y > line;
                        elseif slope > 0
                            condINIT = Y < line;
                        end
                    case 'VerticalSlab'
                        % slab in a matrix
                        assert(length(param)==2,'wrong number of InitializeOPs defining parameters')
                        pos = param(1);
                        halfwidth = param(2);
                        condINIT=abs(X-pos) <= halfwidth;
                    case 'EllipseInMatrix'
                        % elliptical grain in matrix
                        assert(length(param)==4,'wrong number of InitializeOPs defining parameters')
                        semiaxY = param(4);
                        semiaxX = param(3);
                        center = [param(1) param(2)];
                        condINIT= (Y-center(2)).^2/semiaxY.^2 + (X-center(1)).^2/semiaxX.^2 <= 1;
                    case 'RectangleInMatrix'
                        % square grain in matrix
                        % x-y
                        condINIT = (abs(X-ceil(Nx/2)) <= param(1) ) & ( abs(Y-ceil(Ny/2)) <= param(2) );
                    case 'tiltedsquare' % xy-(-x)y
                        pp1 = -(X-ceil(Nx/2))+ sqrt(2)/4*Nx+Ny/2;
                        pp3 = -(X-ceil(Nx/2))- sqrt(2)/4*Nx+Ny/2;
                        pp2 = (X-ceil(Nx/2))- sqrt(2)/4*Nx+Ny/2;
                        pp4 = (X-ceil(Nx/2))+ sqrt(2)/4*Nx+Ny/2;
                        condINIT = Y<pp1 & Y>pp2 & Y>pp3 & Y<pp4;
                    case 'WulffOmega>1'
                        warning('As of 26/8/2020 the ICcode ''WulffOmega>1'' had some problems. Be sure to check what you are doing.' )
                        assert(length(param)==4,'wrong number of InitializeOPs defining parameters')
                        center = [ param(1), param(2)];
                        mean_radius = param(3);
                        soaIE = param(4);
                        assert(soaIE*(in.nfold^2-1)>1,'Omega <= 1. Modify ICparam(4).')
                        if ~isfield(in.intf.params_incl_dep,'offset_ang') % if offset_ang is not one of fields in params_incl_dep
                            in.intf.params_incl_dep.offset_ang = in.misori(1); % implying new input including misorientation
                        end
                        offset_ang_included = false;
                        [fIEijfun,dfIEijfun,~] = AssignAnisotropyFunction(in,offset_ang_included);
                        th = linspace(-pi,pi,300);
                        th = th';
                        th_m = border_angles(in,false);
                        th_for_forbb = linspace(-pi,pi,15e3)';
                        cond_forb_ang = find_forbidden_angles(th_for_forbb,th_m,in);
                        th_mod = th_for_forbb(~cond_forb_ang);
                        % rotate the wulff plot 
    %                     th_wulff = rotate_to_first_inetrval(th,-results_entry.offset_ang); % not sure why this takes the negative of offset ang
                        wx = mean_radius*( fIEijfun(th_mod,soaIE).*cos(th_mod) - dfIEijfun(th_mod,soaIE).*sin(th_mod) );
                        wy = mean_radius*( fIEijfun(th_mod,soaIE).*sin(th_mod) + dfIEijfun(th_mod,soaIE).*cos(th_mod) );
    %                     TH = reshape(atan2(Y,X),[numel(X),1]);
                        [TH,R] = cart2pol(X-center(1),Y-center(2));
                        TH_W = rotate_to_first_inetrval(reshape(TH,[numel(X),1]),-in.offset_ang);
                        RW = CalcRadiusFromXYcontour([wx,wy],TH_W);
                        RW = reshape(RW,[Ny,Nx]);
                        condINIT = R<=RW;
                    case 'Wulff_weak'
                        assert(length(param)==4,'wrong number of InitializeOPs defining parameters')
                        center = [ param(1), param(2)];
                        mean_radius = param(3);
                        Omega = param(4);
                        soaIE = Omega/(in.intf.params_incl_dep.nfold^2-1);
                        if ~isfield(in.intf.params_incl_dep,'offset_ang') % if offset_ang is not one of fields in params_incl_dep
                            in.intf.params_incl_dep.offset_ang = in.misori(1); % implying new input including misorientation
                        end
                        offset_ang_included = false;
                        [fIEijfun,dfIEijfun,~] = AssignAnisotropyFunction(in.intf.params_incl_dep,offset_ang_included);
                        th = linspace(-pi,pi,300);
                        th = th';
                        % rotate the wulff plot 
    %                     th_wulff = rotate_to_first_inetrval(th,-results_entry.offset_ang); % not sure why this takes the negative of offset ang
                        wx = mean_radius*( fIEijfun(th,soaIE).*cos(th) - dfIEijfun(th,soaIE).*sin(th) );
                        wy = mean_radius*( fIEijfun(th,soaIE).*sin(th) + dfIEijfun(th,soaIE).*cos(th) );
    %                     TH = reshape(atan2(Y,X),[numel(X),1]);
                        [TH,R] = cart2pol(X-center(1),Y-center(2));
                        TH_W = rotate_to_first_inetrval(TH(:),-in.misori(1));
                        RW = CalcRadiusFromXYcontour([wx,wy],TH_W);
                        RW = reshape(RW,[Ny,Nx]);
                        condINIT = R<=RW;
                    case 'IS7'
                        % half + circular
                        center = [ -ceil(Nx/4), ceil(Ny/2)];
                        radius = 3/4*Nx;
                        condINIT= (sqrt((Y-center(2)).^2+ (X-center(1)).^2) <= radius) | X < Nx/4;
                    case 'IS8'
                        % ellipse from both sides
                        semiaxY = Ny/1.5;
                        semiaxX = Nx/3;
                        center1 = [1 ceil(Ny/2)]; center2 = [Nx Ny/2];
                        condINIT= (Y-center1(2)).^2/semiaxY.^2 + (X-center1(1)).^2/semiaxX.^2 >= 1;
                        condINIT= condINIT & (Y-center2(2)).^2/semiaxY.^2 + (X-center2(1)).^2/semiaxX.^2 >= 1;
                        condINIT= (Y-center2(2)).^2/semiaxY.^2 + (X-center2(1)).^2/semiaxX.^2 >= 1;
                    case 'IS9'
                        % dents against
                        fracdent = 1/20;
                        dentlength = 5;
                        x1 = Nx*fracdent; x2 = dentlength*x1;
                        ab1 = [x1 1; x2 1]\[1;Ny];
                        ab2 = [x2 1; x1 1]\[1;Ny];
                        condINIT = (Y > ab1(1)*X + ab1(2)) & (Y < ab2(1)*X + ab2(2));
                        x1 = Nx*(1-fracdent); x2 = Nx*(1 - dentlength*fracdent);
                        ab1 = [x1 1; x2 1]\[1;Ny];
                        ab2 = [x2 1; x1 1]\[1;Ny];
                        condINIT = condINIT |  ( (Y > ab1(1)*X + ab1(2)) & (Y < ab2(1)*X + ab2(2)) );
                end % switch ICcode, nOP =2 

            p{2}(condINIT) = 1-roundlim;% initial condition
            p{1}(~condINIT) = 1-roundlim;% initial condition

            case 3 % switch in.nOP==3
                switch in.ICcode
                    case '3junctions'
                        assert(length(param)==2,'wrong number of InitializeOPs defining parameters')
                        assert(param(1)<=Ny/2,'adjust 3-junction position to be <=Ny/2')
                        trijun_pos = param(1);
                        angle = param(2);
                        condINIT_1= X >= Nx/2;
                        XX = X - Nx/2;
                        A = tan(angle);
                        pp1 = A*XX + (Ny-trijun_pos);
                        pp2 = -A*XX + (Ny-trijun_pos);
                        pp3 = A*XX + trijun_pos;
                        pp4 = -A*XX + trijun_pos;
                        condINIT_2 = ( (Y >= pp1) & (Y >= pp2) ) | ( (Y <= pp3) & (Y <= pp4) );
    %                     imagesc(condINIT_2)
                        p{2}(condINIT_1 & ~condINIT_2) = 1-roundlim;% initial condition
                        p{1}(~condINIT_1 & ~condINIT_2) = 1-roundlim;% initial condition
                        p{3}(condINIT_2) = 1-roundlim;

                    case 'SemiCircleOnPlane'
                        assert(length(param)==2,'wrong number of InitializeOPs defining parameters')
                        plane_pos = param(1);
                        radius = param(2);
                        condINIT_1= Y >= plane_pos;
                        center = [ Nx/2, plane_pos];
                        condINIT_2= (sqrt((Y-center(2)).^2+ (X-center(1)).^2) <= radius) & Y > plane_pos;
    %                     imagesc(condINIT_2)
                        p{1}(condINIT_1 & ~condINIT_2) = 1-roundlim;% initial condition
                        p{2}(~condINIT_1 & ~condINIT_2) = 1-roundlim;% initial condition
                        p{3}(condINIT_2) = 1-roundlim;

                    % !!! to bedeprecated
                    case 'SemiCircleOnPlaneDiffuse' 
                        assert(length(param)==2,'wrong number of InitializeOPs defining parameters')
                        plane_pos = param(1);
                        radius = param(2);
                        condINIT_1= Y >= plane_pos;
                        p{2} = CalcDiffuseIntf(in, condINIT_1)
                        center = [ Nx/2, plane_pos];
                        condINIT_2= (sqrt((Y-center(2)).^2+ (X-center(1)).^2) <= radius) & Y > plane_pos;
    %                     imagesc(condINIT_2)
                        p{1}(condINIT_1 & ~condINIT_2) = 1-roundlim;% initial condition
                        p{2}(~condINIT_1 & ~condINIT_2) = 1-roundlim;% initial condition
                        p{3}(condINIT_2) = 1-roundlim;

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

            case 4 % switch in.nOP==4
                switch in.ICcode
                    case 'SemiCircleOnTjunction'
                        assert(length(param)==2,'wrong number of InitializeOPs defining parameters')
                        plane_pos = param(1);
                        radius = param(2);
                        condINIT_1= Y >= plane_pos;
                        center = [ Nx/2, plane_pos];
                        condINIT_2= (sqrt((Y-center(2)).^2+ (X-center(1)).^2) <= radius) & Y > plane_pos;
                        condINIT_3 = X >= Nx/2;
    %                     imagesc(condINIT_2)
                        p{1}(condINIT_1 & ~condINIT_2) = 1-roundlim;% liquid around nucleus
                        p{2}(~condINIT_1 & ~condINIT_2 & condINIT_3) = 1-roundlim;% underlying grain right
                        p{3}(~condINIT_1 & ~condINIT_2 & ~condINIT_3) = 1-roundlim;% underlying grain left
                        p{4}(condINIT_2) = 1-roundlim; % nucleus (semicircle)

                    case 'TestEnergyCalc_4PFs'
                        IWpx = in.IW/in.dx;
                        assert(Nx>=18*IWpx,'error in funcion InitializeOPs: In.ICcode= ''TestEnergyCalc_4PFs'' but Nx<15*IWpx')
                        assert(Ny>=9*IWpx,'error in funcion InitializeOPs: In.ICcode= ''TestEnergyCalc_4PFs'' but Ny<9*IWpx')
                        plane_pos = [Nx 2*Nx]/3;
                        radius = min([0.1*Nx , 1.5/11*Ny]); % 1.5/5*(Nx/3) = 0.1*Nx
                        condINIT_plane1= X <= plane_pos(1);
                        condINIT_plane2= X > plane_pos(1) & X <= plane_pos(2);
                        p{1}(condINIT_plane1) = 1-roundlim;% left slab
                        p{2}(condINIT_plane2) = 1-roundlim;% center slab
                        p{3}(~condINIT_plane1 & ~condINIT_plane2) = 1-roundlim;% right slab
                        % 4 circles
                        centers = [ Nx/2/3, 2.5/9*Ny ; Nx/2/3, 6.5/9*Ny ; Nx/2, 2.5/9*Ny ; 5*Nx/2/3, 2.5/9*Ny ];
                        condINIT_circ14= (sqrt((Y-centers(1,2)).^2+ (X-centers(1,1)).^2) <= radius);
                        condINIT_circ13= (sqrt((Y-centers(2,2)).^2+ (X-centers(2,1)).^2) <= radius);
                        condINIT_circ24= (sqrt((Y-centers(3,2)).^2+ (X-centers(3,1)).^2) <= radius);
                        condINIT_circ34= (sqrt((Y-centers(4,2)).^2+ (X-centers(4,1)).^2) <= radius);
    %                     imagesc(condINIT_2)
                        p{1}(condINIT_circ14 | condINIT_circ13) = roundlim;% left slab without circles (3 and 4)
                        p{2}(condINIT_circ24) = roundlim;% center slab without circle (4)
                        p{3}(condINIT_circ34) = roundlim;% right slab without circle (4)
                        p{3}(condINIT_circ13) = 1-roundlim;% circle (3) within left slab
                        p{4}(condINIT_circ14 | condINIT_circ24 | condINIT_circ34) = 1-roundlim; % circles

                    case '3CirclesInMatrix'
                        % 3 circles in matrix, param = [center1x center1y radius1 ; center2x center2y radius2 ; center3x center3y radius3]
                        assert(numel(param)==9,'wrong number of InitializeOPs defining parameters')
                        centers = [param(1,1) param(1,2) ; param(2,1) param(2,2) ; param(3,1) param(3,2) ];
                        radius = [param(1,3) , param(2,3) , param(3,3)];
                        condINIT_1 = sqrt((Y-centers(1,2)).^2+ (X-centers(1,1)).^2) <= radius(1);
                        condINIT_2 = sqrt((Y-centers(2,2)).^2+ (X-centers(2,1)).^2) <= radius(2);
                        condINIT_3 = sqrt((Y-centers(3,2)).^2+ (X-centers(3,1)).^2) <= radius(3);
                        p{1}(~condINIT_1 & ~condINIT_2 & ~condINIT_3) = 1- roundlim;
                        p{2}(condINIT_2) = 1- roundlim;
                        p{3}(condINIT_1) = 1- roundlim;
                        p{4}(condINIT_3) = 1- roundlim;
                end % switch ICcode 4 OPs
        end % switch nOP
    end % if in.iterno>1
    
    
    
    if isa(in.ICparam,'double')
        disp([num2str(in.nOP) ' phase fields initialized: ''' in.ICcode ''' with parameters [' num2str(param(:)') ']'])
    elseif iscell(in.ICparam)
        disp([num2str(in.nOP) ' phase fields initialized: ''' in.ICcode ''' with parameters [' param{1} ' ,' num2str([param{2} , param{3}]) ']'])
    end
    
    if showIS
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

%% 
function condINIT = JoinVarangCircSegments(in,param,fitfactor,plotting)
    
    [X,Y] = meshgrid((1:in.Nx)-in.Nx/2,(1:in.Ny)-1); % grid in px centered to the middle of domain
    al = param(1)*pi/180;
    ar = param(2)*pi/180;
    angs = param(1:2)'*pi/180;

    % radius of final circle arc with contact angle th (input in radians, output in multiples of initial radius)
    Rth = @(th) sqrt(pi./(2*th-sin(2*th)));
    rad = param(3)/in.dx;
    R = rad*Rth(angs);
    % to set the shape within the domain
    % half-width from 0 to outermost point, different for different contact angle
    halfwidths = R.*heaviside(angs-pi/2-1e-10) + R.*sin(angs).*heaviside(-angs+pi/2+1e-10);
    width = sum(halfwidths);
    height = max(R.*(1-cos(angs)));
%     fitfactor = 0.7;
    tooclose = ( width > fitfactor*in.Nx ) | ( height > fitfactor*in.Ny );
    while tooclose
    rad = rad*0.98;
    R = rad*Rth(angs);
    halfwidths = R.*heaviside(angs-pi/2-1e-10) + R.*sin(angs).*heaviside(-angs+pi/2+1e-10);
    width = sum(halfwidths);
    height = max(R.*(1-cos(angs)));
    tooclose = ( width > fitfactor*in.Nx ) | ( height > fitfactor*in.Ny );
    end

    grtr = angs > [angs(2);angs(1)];
    sign = [1;-1];
    % set 0 in such a way that the asymmetric thing is centered in domain (add halfwidth to left outermost point)
    X = X -(halfwidths(1) - width/2);
    % centers of both circles on vertical line, to be shifted from x=0 later on
    centers = [0 ,-R(1)*cos(al) ;  0 , -R(2)*cos(ar) ];
    for kk = 1:2
    arc{kk} = sqrt((Y-centers(kk,2)).^2+ (X-centers(kk,1)).^2) <= R(kk);
    end
    %                         imagesc(arc{1} + arc{2} ), set(gca,'YDir','normal'),axis equal

    if abs(param(1) - param(2)) < 1e-1
        condINIT = arc{kk};

    else
        % tangent line parameters
        b = (R(1) - R(2) )/(R(1)*cos(al) - R(2)*cos(ar));
        a = sign(grtr)*sqrt(1-b^2);
        c = R(2)*(cos(ar)*b -1);
        % tangent line
        tline = -(a*X+c)/b>Y;
        % make sure that for small and large contat angle tline does not cover the contact point
        if any(param(1:2)<25) && any(param(1:2)>60)
            x_interp_tline = -c/a;
            x_interp_smallang = sign(grtr)*halfwidths(~grtr);
            if abs(x_interp_tline) > abs(x_interp_smallang) % x intercept of tline farther from 0 than contact point
                cshifted = -0.95*x_interp_smallang*a;
                % shift tline but keep former param c for othercomputations
                tline = -(a*X+cshifted)/b>Y;
            end
        end
        
        %                             imagesc(tline), set(gca,'YDir','normal'),axis equal

        % tangent points
        for kk = 1:2
            A = 1 + a^2/b^2;
            B = 2*a/b*( c/b-R(kk)*cos(angs(kk)) );
            C = (c/b)^2  -2*R(kk)*cos(angs(kk))*c/b -R(kk)^2*sin(angs(kk))^2 ;
            tangptx_tworoots =  roots([A,B,C]); 
            assert( abs(B^2-4*A*C)<1e-4 , 'unexpected solution for tangent points. Check InitializaOPs, ICcode: EAVarAng')
            tangpt(kk,1) = real(tangptx_tworoots(1,:));
        end
        tangpt(:,2) = -(tangpt(:,1)*a+c)/b;                        
        % normal lines passing tangent points
        c_nl = b*tangpt(:,1) - a*tangpt(:,2);

        if al<ar
            nline{1} = ( b*X - c_nl(1) )/a<Y;
            nline{2} =( b*X - c_nl(2) )/a>Y;
        else 
            nline{1} = ( b*X - c_nl(1) )/a>Y;
            nline{2} =( b*X - c_nl(2) )/a<Y;
        end
    %     imagesc((tline &  nline{1} & nline{2})), set(gca,'YDir','normal'),axis equal

        if plotting
            if exist('cshifted','var')
                c_tl = cshifted;
            else
                c_tl = c;
            end
            figure(17)
            for kk = 1:2
                plot(X(1,:),sqrt(R(kk)^2-(X(1,:)-centers(kk,1)).^2 )+centers(kk,2),'-') % arcs
                hold on
                plot(X(1,:),(b*X(1,:) - c_nl(kk))/a,'-'), % normal line in tangent point
            end
            axis equal
            plot(X(1,:),-(a*X(1,:) + c_tl)/b,'-k'),  % tangent line
            plot(tangpt(:,1),tangpt(:,2),'xk','LineWidth',1.5) % tangent points
            axis([X(1,1) X(1,end) Y(1,1) Y(end,1)])
            hold off
        end

        % what is below tangent line and not in the arcs
        tangjoint = (tline &  nline{1} & nline{2}) &  ~(arc{1} | arc{2});
        %                             imagesc(tangjoint), set(gca,'YDir','normal'),axis equal

        % throw away right part of the left-angle arc (what is farther than left-arc tangent point OR 0)
        if ( ar > pi/2 && al<pi/2 ) || ( ( ar > pi/2 && al > pi/2 ) && al<ar ) 
            cutoff(1) = 0;
            arc{1}(X>cutoff(1)) = false;
        elseif ( ar < pi/2 && al < pi/2 ) 
            cutoff(1) = tangpt(2,1);
            arc{1}(X>cutoff(1)) = false;
        end

        % throw away left part of the right-angle arc (what is farther than right-arc tangent point OR 0)
        if ( al > pi/2 && ar<pi/2 ) || ( ( ar > pi/2 && al > pi/2 ) && ar<al ) 
            cutoff(2) = 0;
            arc{2}(X<cutoff(2)) = false;
        elseif ( al < pi/2 && ar < pi/2 ) || ( al > pi/2 && ar > pi/2 ) 
            cutoff(2) = tangpt(1,1);
            arc{2}(X<cutoff(2)) = false;
        end

        % tangjoint true outside the grain => eliminate tangjoint behind contact points x coord
        if ( al > pi/2 && ar>pi/2 )
            % shift of contact point from 90 deg to th with conserved volume
            shift = @(th) (1 - sin(th)./sqrt((2*th-sin(2*th))/pi));
            contptx = rad * ([-1;1] + sign.* shift(angs));
            tangjoint((X<contptx(1)) | (X>contptx(2))) = false;
        %                                 imagesc((X<contptx(1)) | (X>contptx(2))), set(gca,'YDir','normal'),axis equal

        end

        condINIT = arc{1} | arc{2} | tangjoint;
%         imagesc(arc{1} | arc{2} | tangjoint ), set(gca,'YDir','normal'),axis equal
    end% if equal angles
end % func
