%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), ?Benchmarking of different strategies to include anisotropy in a curvature-driven multi-phase-field model?, Mendeley Data, V2, doi: 10.17632/5wrv3ky9pp.2> 
% coupled to publication of the same name by Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n? 714754).
% 
%% input_calc_PFpar_dt
% - input structure completion (preceded by in = make_input)
% - calculates the PF parameters and time step length
% - USAGE:
%   -mode 1: 
%       -in = input_calc_PFpar_dt(in)
%   -mode 2: 
%       -in = input_calc_PFpar_dt(in,Courant_nr)
%           if Courant_nr=[] a default value is used

function in = input_calc_PFpar_dt(varargin)
in = varargin{1};

% ___ determine number of phase fields automatically
in.nOP = Assign_nOP_from_ICcode(in.ICcode);

% ___ list of pair-wise interface specificators
[in.IE,in.GBmobility,in.is_locally_aniso_IE,in.is_locally_aniso_L,in.misori,IE_aniso_minmax] = AssignInterfaceProperties(in.nOP,in.intf,in.PFori, in.PFphase);

% ___ max and min interface nergy in the system
IE_minmax = get_IEminmax(in.IE,in.is_locally_aniso_IE,IE_aniso_minmax);
% ___ determine initializing interface eenrgy value such that the narrowest
% interface width is in.IW
in.IEinit = determineIEinit(IE_minmax,in.IE,in.model);

% ___ get phase field parameters for every pair-wise interface from the input
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

% ___ assign Courant number
if length(varargin) == 2 && ~isempty(varargin{2})
    in.Courant_nr = varargin{2};
else % length(varargin) ~= 2 OR isempty(varargin{2})
    in.Courant_nr = GetCourantNumber(in.model,in.is_inclination_dependent_IE,in.intf);
end

% find maximal value of anisotropy function in L
% only 1 type of anisotorpy assumed
maxLaniso = find_Lmax_analytically(in);

% ___ determine time step length
in.dt = in.Courant_nr*in.dx^2/max(maxLaniso*in.Lij.*in.kpp0);
in.Ndt = round(in.simtime/in.dt,-1);
in.anisoNdt =in.Ndt - in.precycle;

%__ select appropriate integer type for time steps depending on Ndt
inttypes = {'int16','int32','int64'};
intmaxs = cell2mat(cellfun(@(x) int64(intmax(x)),inttypes,'UniformOutput',false)');
ind_inttype = min(find(in.Ndt<=intmaxs));
if in.ctrcnt>in.Ndt
%     in.tsteptsctr = int16(1:in.Ndt);
    in.tsteptsctr = eval([inttypes{ind_inttype} '(1:in.Ndt)']);
    warning('in.ctrcnt>in.Ndt, output length set to in.Ndt')
else
%     in.tsteptsctr = int16(floor(linspace(1,in.Ndt,in.ctrcnt)));
    in.tsteptsctr = eval([inttypes{ind_inttype} '(floor(linspace(1,in.Ndt,in.ctrcnt)))']);
end
    
% ___ time steps when phase fields are saved
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

%% Assign_nOP_from_ICcode
function nOP = Assign_nOP_from_ICcode(ICcode)
    codes2OP = {'CircleInMatrix','VerticalPlane','LoadContour','TiltedPlane','TiltedPlaneAng','VerticalSlab','EllipseInMatrix','RectangleInMatrix','Wulff_weak','EAVarAng','Wulff'};
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
