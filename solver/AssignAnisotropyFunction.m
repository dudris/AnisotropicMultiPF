% [fIEijfun,dfIEijfun,ddfIEijfun] = AssignAnisotropyFunction(input)
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

