%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Benchmarking of different strategies to include anisotropy in a curvature-driven multi-phase-field model”, Mendeley Data, V2, doi: 10.17632/5wrv3ky9pp.2> 
% coupled to publication of the same name by Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754).
% 
%% [fIEijfun,dfIEijfun,ddfIEijfun] = AssignAnisotropyFunction(input,offset_ang_included)
function [fIEijfun,dfIEijfun,ddfIEijfun] = AssignAnisotropyFunction(input, offset_ang_included)
% [fIEijfun,dfIEijfun,ddfIEijfun] = AssignAnisotropyFunction(codeIEaniso,is_isotropic,nfold,offset_ang)
    soaIE = input.soaIE;
    nfold = input.nfold;
    codeIEaniso = input.codeIEaniso;

    switch codeIEaniso
        case 'IEanisofun_1'
            if offset_ang_included
                offset_ang = -input.offset_ang;
                fIEijfun = @(phi,soa) ones(size(phi)) + soa*cos(nfold*(phi+offset_ang));
                dfIEijfun = @(phi,soa) -nfold*soa*sin(nfold*(phi+offset_ang));
                ddfIEijfun = @(phi,soa) -nfold*nfold*soa*cos(nfold*(phi+offset_ang));
            else % offset_ang_included = false
                % offset removed because anisofun is calculated in a more complicated way and I need the un-rotated anisofun
                fIEijfun = @(phi,soa) ones(size(phi)) + soa*cos(nfold*phi);
                dfIEijfun = @(phi,soa) -nfold*soa*sin(nfold*phi);
                ddfIEijfun = @(phi,soa) -nfold*nfold*soa*cos(nfold*phi);
            end
    end % switch codeIEaniso
end % function AssignAnisotropyFunction

