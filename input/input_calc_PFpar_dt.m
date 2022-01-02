%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Three strategies of parameters assignment
% in multi phase field model of grain growth with anisotorpic grain boundary properties”, Mendeley Data, 
% v1 http://dx.doi.org/10.17632/5wrv3ky9pp.1>, coupled to publication of the same name by 
% Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% 
%% input_calc_PFpar_dt
% input structure completion (preceded by in = inputfile_GGaniso)
% calculates the PF parameters and time step length
% USAGE:
% mode 1: 
%   in = input_calc_PFpar_dt(in)
% mode 2: 
%   in = input_calc_PFpar_dt(in,Courant_nr)
%       if Courant_nr=[] a default value is used

function in = input_calc_PFpar_dt(varargin)
in = varargin{1};
    
in.nOP = Assign_nOP_from_ICcode(in.ICcode);

[in.IE,in.GBmobility,in.is_locally_aniso_IE,in.is_locally_aniso_L,in.misori,IE_aniso_minmax] = AssignInterfaceProperties(in.nOP,in.intf,in.PFori, in.PFphase);

IE_minmax = get_IEminmax(in.IE,in.is_locally_aniso_IE,IE_aniso_minmax);
in.IEinit = determineIEinit(IE_minmax,in.IE,in.model);

[in.kpp0, in.gam0, in.m, in.Lij, in.IWs, gsq] = get_PF_parameters(in.model, in.IE, in.GBmobility, in.IW,in.IEinit);

if in.is_inclination_dependent_IE
    assert(~in.intf.isStrongAniso,'Only weak anisotropy approximation was validated. Set in.intf.isStrongAniso=false and rerun.')
    in.A = 9/2*gsq; %to simplify anisotropy expression in gamma 
    if strcmp(in.model,'IWc') || strcmp(in.model,'IWvK')
        in.soaKP = in.intf.params_incl_dep.soaIE; % strength of anistropy of kappa equals that of interface energy
    elseif strcmp(in.model,'IWvG')
        in.soaKP = nan;
    end
end % if incl. dep. IE

assert(all(in.gam0)>=0.52 & all(in.gam0)<=40,'some value of gam0 outside interval (0.52,40)')
clear IEresh

if length(varargin) == 2 && ~isempty(varargin{2})
    in.Courant_nr = varargin{2};
else % length(varargin) ~= 2 OR isempty(varargin{2})
    in.Courant_nr = GetCourantNumber(in.model,in.is_inclination_dependent_IE,in.intf);
end

% find maximal value of anisotropy function in L
% only 1 type of anisotorpy assumed
maxLaniso = find_Lmax_analytically(in);

in.dt = in.Courant_nr*in.dx^2/max(maxLaniso*in.Lij.*in.kpp0);
in.Ndt = round(in.simtime/in.dt,-1);
in.anisoNdt =in.Ndt - in.precycle;

if in.ctrcnt>in.Ndt
    in.tsteptsctr = int16(1:in.Ndt);
    warning('in.ctrcnt>in.Ndt, output length set to in.Ndt')
else
    in.tsteptsctr = int16(floor(linspace(1,in.Ndt,in.ctrcnt)));
end
    
if in.outputAtChckpt{1}
    tsteptsChckpt = int16(floor(linspace(1,in.Ndt,(in.outputAtChckpt{2}+1))));
    in.outputAtChckpt{3} = tsteptsChckpt(2:end);
end

%% check some input consistency
assert(size(in.IE,1)==in.nOP*(in.nOP-1)/2,['unexpected number of interfaces declared, length(IE)=' num2str(length(in.IE)) ])

end % func


%% GetCourantNumber
function C = GetCourantNumber(model,is_inclination_dependent_IE,intf)

cond_no_or_weak_aniso = ~is_inclination_dependent_IE || (is_inclination_dependent_IE && intf.params_incl_dep.Omega<=1);
    minmaxratio = min(intf.IE_phases)/max(intf.IE_phases);
    if cond_no_or_weak_aniso
        if minmaxratio<=0.2 && minmaxratio>=0.129
            C = 0.12;
        elseif minmaxratio<=0.3
            C = 0.17;
        else
            C = 0.3;
        end
    else % strong aniso, Omega>1
        
        mods = categorical({'IWc','IWvG','IWvK'});
        condmodel = mods == model;
        
        if intf.params_incl_dep.nfold == 4
            O = intf.params_incl_dep.Omega;
            pars = [1.3800    1.4300    0.9800;     1.9500    1.9100    2.3500 ] ; % 
            C = 1/(pars(1,condmodel)*O+pars(2,condmodel));
            
        else % thorough analysis of stability was only done for 4fold as of 28/10/2021
            % valid for older data, as of 21/10/2021 the stability is much better than that
            pars = [8.27, 2.4704]; 
            dO = intf.params_incl_dep.soaIE*intf.params_incl_dep.Omega;
            C = 1/(pars(1)*dO+pars(2));
            if C>0.3
                C=0.3;
            end
        end
        
        
%         % {3fold ; 4fold  ; 6fold } - fits not applicable for soa>0.4 (4fold)
%         polyparams = {[2.1362   -1.9472    0.5083] ; [1.9074   -1.6546    0.3988] ; [ 14.9730  -19.9730   10.7406   -2.9487    0.3759] };
%         if intf.params_incl_dep.nfold == 3
%              iind = 1;
%              C = polyval(polyparams{iind},intf.params_incl_dep.soaIE);
%         elseif intf.params_incl_dep.nfold == 4
%              iind = 2;
%              C = polyval(polyparams{iind},intf.params_incl_dep.soaIE);
%          elseif intf.params_incl_dep.nfold == 6
%              iind = 3;
%              C = polyval(polyparams{iind},intf.params_incl_dep.soaIE);
%         end % if nfold
    end % if cond_no_or_weak_aniso

end % func


%% 
function IE_minmax = get_IEminmax(IEs,is_locally_aniso_IE,IE_aniso_minmax)

IE_min_all = IEs;
IE_max_all = IEs;

if any(is_locally_aniso_IE)
    IE_min_all(is_locally_aniso_IE) = IEs(is_locally_aniso_IE)*IE_aniso_minmax(1);
    IE_max_all(is_locally_aniso_IE) = IEs(is_locally_aniso_IE)*IE_aniso_minmax(2);
end

IE_minmax(1) = min(IE_min_all);
IE_minmax(2) = max(IE_max_all);

end % func

%% GetPFparameters 
function [kpp0, gam0, m, L,IWout] = GetPFparamsIsotropic(IEres, GBmobility ,IWin)
    gam0 = 1.5;
    m = 6*IEres/IWin;
    kpp0 = (3/4)*IEres*IWin;
    L = (4/3)*GBmobility/IWin;
    IWout = IWin;
%     if is_with_constant_IW
%         [kpp0, gam0,~, m, L] = parameters(IEres, GBmobility ,IWin,IEinit);
%         IWout = IWin;
%     else
%         % kpp0 = kappa ... constant computed from IWin and IEinit
%         [kpp0, gam0, m, L, IWout, ~] = parameters_varIW(IEres,GBmobility,IWin,IEinit);
%     end% if
    
%     % test
%     [gamma,g_function,sqrt_f0_function] = load_gfunc_sqrtf0;
%     ggam = interp1(gamma,g_function,gam0_inp);
%     sqrtf0gam = interp1(gamma,sqrt_f0_function,gam0_inp);
%     rel_diff_from_expected.m = abs(6-1/ggam/sqrtf0gam)/6 ;
%     rel_diff_from_expected.kpp0 = abs(3/4-sqrtf0gam/ggam)/(3/4) ;
%     rel_diff_from_expected.L = abs(4/3-ggam/sqrtf0gam)/(4/3) ;
end % func GetPFparamsIsotropic

function [m, Lij, kpp0, gam0, IWout] = GetPFparamsMisorientationOnly(model,IEres, GBmobility ,IWin,IEinit)
    if strcmp(model,'IWc')
        [kpp0, gam0,~, m, Lij] = parameters(IEres, GBmobility ,IWin,IEinit);
         IWout = IWin;
         assert(all(gam0>=0.52) && all(gam0<=40),['IWc, error in ''GetPFparamsMisorientationOnly'' : IEinit = ' num2str(IEinit) ' and some of gam0 < 0.52 or gam0 > 40. gam0 = [' num2str(gam0) ']'])
    elseif strcmp(model,'IWvG')
        % kpp0 = kappa ... constant computed from IW and IEinit
        [kpp0, gam0, m, Lij, IWout, ~] = parameters_varIW(IEres,GBmobility,IWin,IEinit);
    elseif strcmp(model,'IWvK')
        % gam0 = 1.5
        [kpp0, gam0, m, Lij,IWout] = parameters_IWvK(IEres, GBmobility ,IWin,IEinit);
    end %if
end

function [m, L, kpp0, gam0, soaKP, IWout, A] = GetPFparamsAnisoWeak(model,IEres, GBmobility ,IWin,IEinit, soaIE)
    if strcmp(model,'IWc')
        [kpp0, gam0,gsq, m, L] = parameters(IEres, GBmobility ,IWin,IEinit);
        soaKP = soaIE;
         A = 9/2*gsq; %to simplify anisotropy expression in gamma for small anisotropies
         IWout = IWin;
    elseif strcmp(model,'IWvG')
        % kpp0 = kappa ... constant computed from IW and IEinit
        [kpp0, gam0, m, L, IWout, gsq] = parameters_varIW(IEres,GBmobility,IWin,IEinit);
        soaKP = nan;
        A = 9/2*gsq;
    elseif strcmp(model,'IWvK')
        [kpp0, gam0, m, L,IWout] = parameters_IWvK(IEres, GBmobility ,IWin,IEinit);
        soaKP = soaIE;
        A = nan;
    end %if
end % func GetPFparamsAnisoWeak

function [kpp0, gam0, m, Lij,IWs] = parameters_IWvK(IEresh, GBmobility ,IWin,IEinit)
    gam0 = 1.5;
    m = 6*IEinit/IWin;
    kpp0 = (9/2)*IEresh.^2/m;
    IWs = 6*IEresh/m;
%     IWs = (2/3)*IEresh/m;
    Lij = (4/3)*GBmobility./IWs;
end% func

%% Assign_nOP_from_ICcode
function nOP = Assign_nOP_from_ICcode(ICcode)
    codes2OP = {'CircleInMatrix','VerticalPlane','LoadContour','TiltedPlane','TiltedPlaneAng','VerticalSlab','EllipseInMatrix','RectangleInMatrix','Wulff_weak','EAVarAng'};
    codes3OP = {'3junctions','SemiCircleOnPlane', 'SemiCircleOnPlaneDiffuse' ,'Tjunction','2CirclesInMatrix'};
    codes4OP = {'SemiCircleOnTjunction','SemiCircleOnTjunctionDiffuse','TestEnergyCalc_4PFs','3CirclesInMatrix'};
    if any(strcmp(ICcode,codes2OP))
        nOP = 2;
    elseif any(strcmp(ICcode,codes3OP))
        nOP =3;
    elseif any(strcmp(ICcode,codes4OP))
        nOP =4;
    end
end% func   

%%
function maxLaniso = find_Lmax_analytically(in)
    if in.is_inclination_dependent_L
        % anisotropy funciton of L
        % making sure the anisofunction takes 'th' as input
        buffLfun = @(phi) in.intf.params_incl_dep.Lfun(phi);
        Lf = sym(buffLfun);
%         subplot(211)
%             fplot(Lf,[-pi,pi])
%         subplot(212)
%             fplot(diff(Lf,1),[-pi,pi])
%             grid on

        syms phi
        % get x values of zero derivative (periodic function expected)
        [solx, param, cond] = solve(diff(Lf,1)==0, phi, 'ReturnConditions', true);
        numsol = size(solx,1);
        assume(cond)
        % get values of parameter in periodic solution in (-pi,pi) interval
        interval = [solx>-pi , solx<pi];
        solk = cell(numsol,1);
        for kk = 1:numsol
            solk{kk} = solve(interval(kk,:), param);       
        end
        
        nzsol = cellfun(@(x) ~isempty(x),solk);
        
        % get x values in (-pi,pi) interval where diff(Lf,1)==0 
        valx = subs(solx(nzsol), param, solk(find(nzsol)));
        % valx where 2nd derivative is <0 and hence with local maximum
        maxx = valx(simplify(subs(diff(Lf,2),valx)<0));
        % maximal value of Lf in found points 
        maxLaniso = max(subs(Lf,maxx));
        % in MATLAB double precision
        maxLaniso = double(maxLaniso);
    else
        maxLaniso = 1;
    end
end
