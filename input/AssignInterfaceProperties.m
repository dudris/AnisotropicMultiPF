%% AssignInterfaceProperties
% [IE,is_locally_isotropic] = AssignInterfaceProperties(nOP,intf,ind_is_solid)
% - to each pair-wise interface assigns mean IE, misorientation angles and a bool containing information whether that interface is inclination dependent
% - for nOP there are N = (nOP-1)nOP/2 interfaces
% IE                         ... vector of mean interface energies size(IE) = [N,1],
% is_locally_aniso_IE ... vector of bools,  size(is_locally_aniso_IE) = [N,1]
% is_locally_aniso_L ... vector of bools,  size(is_locally_aniso_L) = [N,1]
% misori                 ... vector of misorientation angles,  size(misori) = [N,1]

function [IE,is_locally_aniso_IE,is_locally_aniso_L,misori] = AssignInterfaceProperties(nOP,intf,ind_is_solid,PFori)
    npairs = nOP*(nOP-1)/2;
    indpairs = combnk(1:nOP,2);
%     is_locally_isotropic = false(npairs,1);
    is_locally_aniso_IE = false(npairs,1);
    is_locally_aniso_L = false(npairs,1);
    misori = zeros(npairs,1);
    if any(intf.is_incl_dep_IE) % incl. dep IE
        IE = ones(npairs,2);
    else
        IE = ones(npairs,1);
    end
    
    single_phase = false;
    if length(ind_is_solid)==nOP || isempty(ind_is_solid)
        single_phase = true;
    end
    
    for k = 1:npairs
        i = indpairs(k,1);
        j = indpairs(k,2);
        is_solid(1) = any(i == ind_is_solid);
        is_solid(2) = any(j == ind_is_solid);
        assert(any(is_solid),'Check the simulation input: two PFs of liquid phase encountered.')
        if all(is_solid) && ~single_phase
            intf_type = 2;
        elseif any(is_solid) 
            intf_type = 1;
        end
        
        is_locally_aniso_IE(k) = intf.is_incl_dep_IE(intf_type);
        is_locally_aniso_L(k) = intf.is_incl_dep_L(intf_type);
%         is_locally_isotropic(k) = is_locally_isotorpic_IE && ~intf.is_incl_dep_L(intf_type) ;
        misori(k) = PFori(j)- PFori(i); % CAREFUL - the order ij or ji could be relevant
        IE_mean(k) = intf.IE_phases(intf_type);
            
        if is_locally_aniso_IE(k) % locally inclination-dep. IE
            IE(k,[1,2]) = IE_mean(k)*([1 1] + intf.params_incl_dep.soaIE*[1 -1]); % [IE_max IE_min] for each interface , i.e. (npairs x 2)
        elseif ~is_locally_aniso_IE(k) && any(intf.is_incl_dep_IE) % locally isotropic intf in a system with some incl. dep IE
            IE(k,1) = AssignIE_locally_iso(intf,misori(k),IE_mean(k)); 
            IE(k,2) = IE(k,1);
        else % is locally IE-isotropic and no 
            IE(k,1) = AssignIE_locally_iso(intf,misori(k),IE_mean(k)); 
        end % if locally isotropic IE
    end % for npairs
    
    
end% func

%% test

% nOP = 5 ;
% ind_is_solid = 2:5;
% intf.IE_phases = [0.3 0.03]; % [SL , GB]
% intf.is_incl_dep = [false false]; % [SL , GB]
% intf.is_incl_dep = [true false]; % [SL , GB]
% intf.isStrongAniso = false;
% [IE,is_locally_isotropic] = AssignInterfaceProperties(nOP,intf,ind_is_solid)
% indpairs = combnk(1:nOP,2);
% T = table(indpairs,IE, is_locally_isotropic)

%%
function IEk = AssignIE_locally_iso(intf,misorik,IEmeank)
    switch intf.params_misor_dep.code
        case 'none'
            IEk = IEmeank;
        case 'readshockley'
            ratio = abs(misorik)/intf.params_misor_dep.width;
            if ratio>1
                IEk = IEmeank;
            else
                IEk = IEmeank*ratio*(1-log(ratio));
            end
    end
end

%% test single phase
% nOP = 50
% ind_is_solid = 1:50;
% PFori = linspace(pi/4,-pi/4,nOP); % uniform orientations distribution
% PFori = rand(nOP,1); % uniform orientations distribution
% intf.is_incl_dep = [false];
% intf.IE_phases = 0.3;
% intf.params_misor_dep.code = 'readshockley';
% intf.params_misor_dep.width = 15*pi/180;

% [IE,is_locally_isotropic,misori] =AssignInterfaceProperties(nOP,intf,ind_is_solid,PFori);
% w = intf.params_misor_dep.width;
% [misori,indsort] = sort(misori);
% indcusp = abs(misori)<=w;
% plot(misori,IE(indsort),'o',misori(indcusp),0.3*abs(misori(indcusp))/w.*(1-log(abs(misori(indcusp))/w)),'-k',misori(~indcusp),0.3*ones(sum(~indcusp),1),'-k'), grid on