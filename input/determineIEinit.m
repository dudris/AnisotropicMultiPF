%% determineIEinit
% INPUT
%   IEs ... interface energy of every pair-wise interface in J/m^2
%       - if any(intf.is_incl_dep_IE)
%           - size(IEs) = [npairs,2]
%           - k-th interface has inclination-dependent IE, then
%           IEs(k,1) = minimal reached IE, IEs(k,2) = maximal reached IE
%       else 
%           size(IEs) = [npairs,1]
%       - if all interfaces have equal IE, then IEinit=unique(IEs)
%       irrespective of the model

function IEinit = determineIEinit(IEminmax,IEs,model)

% IEminmax(1) = min(min(IEs));
% IEminmax(2) = max(max(IEs));

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