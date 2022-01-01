%% AssignInterfaceProperties
% [IE,is_locally_isotropic] = AssignInterfaceProperties(nOP,intf,ind_is_solid)
% - to each pair-wise interface assigns mean IE, misorientation angles and a bool containing information whether that interface is inclination dependent
% - for nOP there are npairs = (nOP-1)nOP/2 interfaces
% IE                         ... vector of mean interface energies size(IE) = [npairs,1],
% is_locally_aniso_IE ... vector of bools,  size(is_locally_aniso_IE) = [npairs,1]
% is_locally_aniso_L ... vector of bools,  size(is_locally_aniso_L) = [npairs,1]
% misori                 ... vector of misorientation angles,  size(misori) = [npairs,1]

function [IE,GBmobility,is_locally_aniso_IE,is_locally_aniso_L,misori,IE_aniso_minmax] = AssignInterfaceProperties(nOP,intf,PFori, PFphase)
    npairs = nOP*(nOP-1)/2;
    indpairs = combnk(1:nOP,2); % all possible 2-element combinatinos os numbers 1:nOP
%     is_locally_isotropic = false(npairs,1);
    is_locally_aniso_IE = false(npairs,1);
    is_locally_aniso_L = false(npairs,1);
    misori = zeros(npairs,1);
    GBmobility = zeros(npairs,1);
%     if any(intf.is_incl_dep_IE) % incl. dep IE
%         IE = ones(npairs,2); % motivated by 'strong anisotropy
%         approximation', not needed in weak
%     else
%         IE = ones(npairs,1);
%     end
    IE = ones(npairs,1);
    
    for k = 1:npairs
        i = indpairs(k,1);
        j = indpairs(k,2);
        
        i_phase = PFphase(indpairs(k,1));
        j_phase = PFphase(indpairs(k,2));
        is_locally_aniso_IE(k) = intf.is_incl_dep_IE(i_phase,j_phase);
        is_locally_aniso_L(k) = intf.is_incl_dep_L(i_phase,j_phase);
        misori(k) = PFori(j)- PFori(i); % CAREFUL - the order ij or ji could be relevant
        IE_mean(k) = intf.IE_phases(i_phase,j_phase);
        GBmobility(k) = intf.mob_phases(i_phase,j_phase);
            
        if is_locally_aniso_IE(k) % locally inclination-dep. IE
            IE(k,1) = IE_mean(k);
%             IE(k,[1,2]) = IE_mean(k)*([1 1] + intf.params_incl_dep.soaIE*[-1 1]); % [IE_min IE_max] for each interface , i.e. (npairs x 2)
        elseif ~is_locally_aniso_IE(k) && any(intf.is_incl_dep_IE(:)) % locally isotropic intf in a system with some incl. dep IE
            IE(k,1) = AssignIE_locally_iso(intf,misori(k),IE_mean(k)); 
%             IE(k,2) = IE(k,1);
        else % is locally IE-isotropic and no 
            IE(k,1) = AssignIE_locally_iso(intf,misori(k),IE_mean(k)); 
        end % if locally isotropic IE
    end % for npairs
    
    if any(any(intf.is_incl_dep_IE))
        % min of anisotorpy function
        IE_aniso_minmax(1,1) = 1-intf.params_incl_dep.soaIE;
        % max of REGULARIZED anisotorpy function (equal to that of non-regularized for weak anisotropy)
        th_m = border_angles(intf.params_incl_dep,false);
        [fIE,~,~] = AssignAnisotropyFunction(intf.params_incl_dep,false);
        IE_aniso_minmax(1,2) = fIE(th_m(1),intf.params_incl_dep.soaIE)/cos(th_m(1));
    else
        IE_aniso_minmax = [];
    end
    
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