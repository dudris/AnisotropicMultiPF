%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Benchmarking of different strategies to include anisotropy in a curvature-driven multi-phase-field model”, Mendeley Data, V2, doi: 10.17632/5wrv3ky9pp.2> 
% coupled to publication of the same name by Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754).
% 
%% get_PF_parameters
% [kpp0, gam0, m, L, IWout, gsq] = get_PF_parameters(model,IEs, GBmobility ,IWmin, IEinit)
% - the narrowest interface in models with variable interface width to be IWmin
% - the used formulas are described in Supplementary material the Paper
% 
% INPUT
%   - model ... string vector, either 'IWc', 'IWvG' or 'IWvK' (as in Minar, !!!, 2022)
%       'IWc' ... model with CONSTANT interafce width (as in Moelans, Phys.Rev.B, 2008)
%       'IWvG' ... model with VARIABLE interafce width and all anisotropy in parameter GAMMA (as in Ravash, !!!!, 2017)
%       'IWvK' ... model with VARIABLE interafce width and all anisotropy in parameter KAPPA (as in Minar, !!!, 2022)
%   - IEs ... vector, interface energies for which the phase field (PF) parameters are determined
%   - GBmobility ... vector, grain boundary mobilities for which the phase field (PF) parameters are determined
%   - IWmin ... scalar, minimal interface width in meters
%   - IEinit ... scalar, returned by function 'determineIEinit'
% note that it should be size(IEs)==size(GBmobility)
% 
% OUTPUT
%   - kpp0 ... vector, size(kpp0)=size(IEs), values of parameter kappa for
%   each interface (a constant in IWvG)
%   - gam0 ... vector, size(gam0)=size(IEs), values of parameter gamma for
%   each interface (a constant in IWvK, gam0=1.5 )
%   - m ... a constant, parameter m
%   - L ... vector, size(L)=size(IEs), values of parameter L for each interface 
%   - IWout ... size(IWout)=size(IEs), interface width for each interface
%   (a constant in IWc)
%   - gsq ... vector size(gsq)=size(IEs), vaue of (g(gamma)).^2, needed in
%   'weak anisotorpy approximation' used in the paper in calucation of
%   values of gamma at interfaces with inclination-dependent IE (not used in IWvK)

function [kpp0, gam0, m, L, IWout, gsq] = get_PF_parameters(model,IEs, GBmobility ,IWmin, IEinit)

    if strcmp(model,'IWc')
        IWout = IWmin;
        % ___ polynomial sqrt(f0)*g -> 1/gam; fitted for gamma 0.52-40
        prodgsqrtf0c_gaminv_coef = [103.397      -165.393      105.3469     -44.55661       24.7348     -11.25718      1.999642];
        % ___ polynomial 1/gam -> sqrt(f0(1/gam)); fitted for gamma 0.52-40
        sqrtf0_gaminv_coef = [-0.072966     0.35784    -0.68325     0.63578    -0.48566     0.53703];
        
        m = 6*IEinit/IWmin;
        prodgsqrtf0c = IEs/6/IEinit;
        gamma_inverse = polyval(prodgsqrtf0c_gaminv_coef, prodgsqrtf0c);
        gam0 = 1./gamma_inverse;

        sqrtf0c = polyval(sqrtf0_gaminv_coef,1./gam0);
        kpp0 = sqrtf0c.^2.*6*IEinit*IWmin;
        
        L = GBmobility.*IEs./kpp0;
        gsq = IEs.^2./kpp0/m;
        
    elseif strcmp(model,'IWvG')
        % ___ both fitted for 0.52 <= gamma <= 40
        gaminv_gsq_coef = [-5.73008      18.8615     -23.0557      7.47952      8.33568     -8.01224      2.00013];
        sqrtf0_gaminv_coef = [-0.072966     0.35784    -0.68325     0.63578    -0.48566     0.53703];

        % ___ iteration to find such gammas to give correct IE AND IWmin
        for ctr = 0:1
            if ctr == 0
                IWiter = IWmin;
            end
            kpp0 = 3/4 *IEinit*IWiter ;
            m = 6*IEinit/IWiter;

            gsq = IEs.^2/(m*kpp0);

            gamma_inverse =  polyval(gaminv_gsq_coef,gsq);
            sqrtf_interf =  polyval(sqrtf0_gaminv_coef,gamma_inverse);

            IWiter = IWmin*sqrt(8)*max(unique(sqrtf_interf)); % IW_0 = IWmin*sqrt(f_0(gam_max)/f_0(1.5)) = IWmin*sqrt(8*f_0(gam_max))
        end
        
        IWout = sqrt(kpp0./m)./sqrtf_interf;
        gam0 = 1./gamma_inverse;
        L = GBmobility.*IEs/kpp0;
        
    elseif strcmp(model,'IWvK')
        gsq = nan;
        
        gam0 = 1.5;
        m = 6*IEinit/IWmin;
        kpp0 = (9/2)*IEs.^2/m;
        IWout = 6*IEs/m;
        L = (4/3)*GBmobility./IWout;
        
    end

end % func


