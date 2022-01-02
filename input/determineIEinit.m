%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Three strategies of parameters assignment
% in multi phase field model of grain growth with anisotorpic grain boundary properties”, Mendeley Data, 
% v1 http://dx.doi.org/10.17632/5wrv3ky9pp.1>, coupled to publication of the same name by 
% Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% 
%% determineIEinit
% IEinit = determineIEinit(IEminmax,IEs,model)
% - the function checks that the anisotropy is not too strong to compute
% the parameters
% 
% INPUT
%   - IEminmax ... minimal and maximal interface energy in the system (may not be trivial with inclination dependence)
%   - IEs ... mean interface energy of every pair-wise interface
%   - model ... either of 'IWc', 'IWvG' or 'IWvK', see Paper for more details
% OUTPUT
%   - IEinit ... scalar
%       - if all interfaces have equal IE, then IEinit=unique(IEs)
%       irrespective of the model

function IEinit = determineIEinit(IEminmax,IEs,model)

    if strcmp(model,'IWc')
        g_function	= [0.098546		0.765691]'; % g(gamma) for gamma=0.52 and 40
        sqrt_f0 = [	0.07014984		0.5309548]'; % sqrt(f0(gamma)) for gamma=0.52 and 40
        Gf = g_function([1,end]).*sqrt_f0([1,end]);
        IEinit_min = max(IEminmax)./6/Gf(2); 
        IEinit_max = min(IEminmax)./6/Gf(1); 
        assert(min(IEminmax)/max(IEminmax)>(Gf(1)/Gf(2)),['input_calc_PFpar_dt>determineIEinit error msg: IEmin/IEmax too small (' num2str(min(IEminmax)/max(IEminmax)) '). Must be >' num2str(Gf(1)/Gf(2))])
        IEinit = mean(unique(IEs)); % this usually yields convenient values of gamma not far from 1.5
        % for cases when the above IEinit does mot meet the inequalities
        if (IEinit > IEinit_max) || (IEinit < IEinit_min)
            % to be replaced by value in between the limits
            IEinit = mean([IEinit_min,IEinit_max]);
        end

    elseif strcmp(model,'IWvG')
        if length(unique(IEminmax)) == 1
            IEinit = unique(IEminmax);

        else
            % G = g(gamma)
            g_function	= [0.098546		0.765691]'; % g(gamma) for gamma=0.52 and 40
            assert(min(IEminmax)/max(IEminmax)>(g_function(1)/g_function(2)),['input_calc_PFpar_dt>determineIEinit error msg: IEmin/IEmax too small (' num2str(min(IEminmax)/max(IEminmax)) '). Must be >' num2str(g_function(1)/g_function(2))])
            IEinit_lim_min = max(IEminmax)*sqrt(2/9)/g_function(2);
            IEinit_lim_max = min(IEminmax)*sqrt(2/9)/g_function(1);
        %     IEinit = mean(IEinit_lim);
            IEinit = IEinit_lim_min+0.1*(IEinit_lim_max-IEinit_lim_min);
        end% if isotropic
        
    elseif strcmp(model,'IWvK')
        IEinit = min(IEminmax);
        
    end % if model
end