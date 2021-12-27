%input_calc_PFpar_dt
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

% [in.IE,in.is_locally_isotropic,in.misori] = AssignInterfaceProperties(in.nOP,in.intf,in.ind_is_solid,in.PFori);
[in.IE,in.GBmobility,in.is_locally_aniso_IE,in.is_locally_aniso_L,in.misori] = AssignInterfaceProperties(in.nOP,in.intf,in.PFori, in.PFphase);

% in.GBmobility = 7.5e-16*ones(size(in.IE));

in.IEinit = determineIEinit(in.intf.IE_phases(:),in.model,in.intf,in);
% decides on how maxmin IE are treated
IEresh = mean(in.IE,2);


if (~in.is_inclination_dependent_IE && ~in.is_misori_dependent ) % isotropic, single interface type
    [in.kpp0, in.gam0, in.m, in.Lij,in.IWs] = GetPFparamsIsotropic(IEresh, in.GBmobility ,in.IW);
    
elseif ~in.is_inclination_dependent_IE && in.is_misori_dependent % pairwise isotropic
    [in.m, in.Lij, in.kpp0, in.gam0, in.IWs] = GetPFparamsMisorientationOnly(in.model,IEresh, in.GBmobility ,in.IW,in.IEinit);
    
elseif in.is_inclination_dependent_IE && in.is_misori_dependent % inclination dependent multiPF
    [in.m, in.Lij, in.kpp0, in.gam0, in.soaKP, in.IWs, in.A] = GetPFparamsAnisoWeak(in.model,IEresh, in.GBmobility ,in.IW,in.IEinit, in.intf.params_incl_dep.soaIE);
    
%     error('Mighty creator asks: is this really implemented? ... in.is_inclination_dependent_IE=true & in.is_misori_dependent=true')
elseif in.is_inclination_dependent_IE && in.intf.isStrongAniso
    [in.m, in.Lij,  in.kpp0, in.gam0, in.soaKP, in.soaGM,in.IWout] = GetPFparamsAnisoStrong(in.model,IEresh, in.GBmobility ,in.IW,in.IEinit);
    in.IE = mean(IE,2); % ??
    
elseif ~in.intf.isStrongAniso || in.is_misori_dependent
    [in.m, in.Lij, in.kpp0, in.gam0, in.soaKP, in.IWs,in.A] = GetPFparamsAnisoWeak(in.model,IEresh, in.GBmobility ,in.IW,in.IEinit, in.intf.params_incl_dep.soaIE);
    
end
assert(all(in.gam0)>=0.52 & all(in.gam0)<=40,'some value of gam0 outside interval (0.52,40)')
clear IEresh

if length(varargin) == 2 && ~isempty(varargin{2})
    in.Courant_nr = varargin{2};
else % length(varargin) ~= 2 OR isempty(varargin{2})
    in.Courant_nr = GetCourantNumber(in.model,in.is_inclination_dependent_IE,in.intf);
%     in.Courant_nr=0.05;
    % Courant_nr=0.05;
    % Courant_nr = Courant_nr*0.5;
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

%% AssignGPdifferencePrefactor
function pref_GP_diff = AssignGPdifferencePrefactor(nOP,PF_to_conserve)
    pref_GP_diff = -ones(nOP,2);
    pref_GP_diff(PF_to_conserve{1},1) =1;
    pref_GP_diff(PF_to_conserve{2},2) =1;
end

%% determineIEinit
function IEinit = determineIEinit(IE_phases,model,intf,in)
IE_phases = unique(IE_phases(~isnan(IE_phases)));
    if strcmp(model,'IWc')
%         if min(IE_phases)/max(IE_phases)<0.2
%             IEinit = 0.5*max(IE_phases);
%         else
%             IEinit = max(IE_phases);
%         end
%         IEinit = min(IE_phases); % only working for cca 0.45<IEratio<5
%         IEinit = max(IE_phases); % working for cca 0.0.03<IEratio<22 
%         IEinit = 0.5;
        IEinit = mean(IE_phases); % works well for interval IEratio cca 0.03 to 45
    
    elseif strcmp(model,'IWvG')
% % %         IEinit = determineIEinit_varIW(IW,intf,plotting);
            
        if any(intf.is_incl_dep_IE) % inclination dependent
            % taking mean IE to compute mean PF parameters
            if strcmp(in.PFpar_compspec,'mean')
                IEinit = determineIEinit_varIW(IE_phases); 
            % taking into account largest&smallest aniso IE to compute mean PF parameters
            elseif strcmp(in.PFpar_compspec,'fullaniso')
                % maximal reached aniso energy is meanIE*(1 + soaIE)
                maxIEaniso = max(intf.IE_phases(intf.is_incl_dep_IE))*(1 + intf.params_incl_dep.soaIE);
                % minimal reached aniso energy is meanIE*(1 - soaIE)
                minIEaniso = min(intf.IE_phases(intf.is_incl_dep_IE))*(1 - intf.params_incl_dep.soaIE);
                % some other isotropic interface can have smaller IE
                maxIE = max([maxIEaniso, intf.IE_phases(~intf.is_incl_dep_IE)]);
                minIE = min([minIEaniso, intf.IE_phases(~intf.is_incl_dep_IE)]);
                IEminmax =  [minIE, maxIE];
                % minimal and maximal IE of the system needed
            end % PF par computation specification
            
        else
            IEminmax =  [min([intf.IE_phases]), max([intf.IE_phases])];
        end % if inclination dep
        IEinit = determineIEinit_varIW(IEminmax);
        
    elseif strcmp(model,'IWvK')
        
        if any(intf.is_incl_dep_IE) % inclination dependent
            % taking mean IE to compute mean PF parameters
            if strcmp(in.PFpar_compspec,'mean')
                IEinit = min(intf.IE_phases);
            elseif strcmp(in.PFpar_compspec,'fullaniso')
    %         to account for the smallest with inclination dep. below
                % minimal reached aniso energy is meanIE*(1 - soaIE)
                minIEaniso = min(intf.IE_phases(intf.is_incl_dep_IE))*(1 - intf.params_incl_dep.soaIE);
                % some other isotropic interface can have smaller IE
                IEinit = min([minIEaniso, intf.IE_phases(~intf.is_incl_dep_IE)]);
            end% PF par computation specification
        else
            IEinit = min(intf.IE_phases);
        end % if inclination dep
        
    end % if model
end
%% determineIEinit_varIW
% find the most convenient IEinit - such that gamma is within reasonable
% bounds irrespective of the ratio of IEmax/IEmin
function SIGMA_INIT = determineIEinit_varIW(IE_phases)
    if length(unique(IE_phases)) == 1
        SIGMA_INIT = unique(IE_phases);

    else
        % G = g(gamma)
        G1 = 0.098546; % gamma = 0.52
        G2 = 0.765691 ; %  gamma=40
        assert(min(IE_phases)/max(IE_phases)>G1/G2,'input>determineIEinit_varIW msg: IEmin/IEmax too small')
        IEinit_lim(1) = max(IE_phases)*sqrt(2/9)/G2;
        IEinit_lim(2) = min(IE_phases)*sqrt(2/9)/G1;
        assert(IEinit_lim(1)<IEinit_lim(2),'input>determineIEinit_varIW msg: something is wrong, IEinit_lim(1)>=IEinit_lim(2)')
    %     SIGMA_INIT = mean(IEinit_lim);
        SIGMA_INIT = IEinit_lim(1)+0.1*(IEinit_lim(2)-IEinit_lim(1));
    end% if isotropic
end% func

% %% get_IEresh
% function IEresh = get_IEresh(in)
%     
%     if in.is_inclination_dependent_IE
%         if strcmp(in.model,'IWc')
%             IEresh = mean(in.IE,2);
%             
%         elseif strcmp(in.model,'IWvG')
%             IEresh = max(in.IE,[],2);
%             
%         elseif strcmp(in.model,'IWvK')
%             IEresh = min(in.IE,[],2);
%         end
%     else
%         IEresh = in.IE(:);
%     end
% end

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

function [m, L,  kpp0, gam0, soaKP, soaGM,IWout] = GetPFparamsAnisoStrong(model,IEres, GBmobility ,IWin,IEinit)
    if any(size(IEres)~=size(GBmobility))
        GBmobility = ones(size(IEres))*GBmobility(1);
    end
    if strcmp(model,'IWc')
        % check reshape for more PFs
        [kpp00, gam00,~, m, L] = parameters(IEres, GBmobility ,IWin,IEinit);
        L = mean(L); % especially here
        % !!!! should not regularized anisotropy function be used here to compute the soa's????
        soaKP = (kpp00(:,1)-kpp00(:,2))/(kpp00(:,1)+kpp00(:,2));
        soaGM = (gam00(:,1)-gam00(:,2))/(gam00(:,1)+gam00(:,2)); % (max-min)/(max+min)
        kpp0 = mean(kpp00,2);
        gam0 = mean(gam00,2);
        IWout = IWin;
    elseif strcmp(model,'IWvG')
        [kpp00, gam00, m, L, IWout, ~] = parameters_varIW(IEres,GBmobility,IWin,IEinit);
        soaKP = nan;
        soaGM = (gam00(:,1)-gam00(:,2))/(gam00(:,1)+gam00(:,2)); % (max-min)/(max+min)
        kpp0 = mean(kpp00,2);
        gam0 = mean(gam00,2);
    elseif strcmp(model,'IWvK')
        error('Meaningless specification: is_StrongAniso=true but model=IWvK. Set false.')
    end
    
    if max(gam0)>40
        error('some gam0 parameter larger than 4 - check it out')
    end
end % func GetPFparamsAnisoStrong


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
