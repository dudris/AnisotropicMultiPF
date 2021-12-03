%% calcDifferentialOperators
% [gradX,gradY,gradXX,gradYY,gradYX,lap,simsize] = calcDifferentialOperators(in)
% except for 'simsize' all are function handles containing matrix multiplication
% simsize is 1x2 vector with system dmensions (either [NxNy,1] or [Ny,Nx] depending on 'solvermethod' switch)
function [gradX,gradY,gradXX,gradYY,gradYX,lap,simsize] = calcDifferentialOperators(in)
    solvermethod = in.solvermethod;
    laplacianmethod = in.laplacianmethod;
    
    if strcmp(in.BCs,'Nspecial')
        BCs = 'N';
    else
        BCs = in.BCs;
    end
    Nx = in.Nx;
    Ny = in.Ny;
    dx = in.dx;
    dy = in.dy;
    if strcmp(solvermethod,'1D')
        % laplacean matrix <=> div grad, chcked against built-in function 'laplacian'
        % laplacian(S)= lapSx + lapSy = S*lapX + lapY*S
        if strcmp(BCs,'P')
            BCcell = {'P' []};
        elseif strcmp(BCs,'N')
            BCcell = {'N' 'N'};
        end
        gradXXmatrix = laplacian1D(Nx,dx, BCcell,'x'); % S*lapX = lapSx
        gradYYmatrix = laplacian1D(Ny,dy, BCcell,'y'); % lapY*S = lapSy
        gradXmatrix = create_cdiff(Nx,dx,'x'); % S*gradXmatrix = gradSx 
        gradYmatrix = create_cdiff(Ny,dy,'y'); % gradYmatrix*S = gradSy
        gradX = @(S) S*gradXmatrix;
        gradY = @(S) gradYmatrix*S;
        gradYX = @(S) (gradYmatrix*S)*gradXmatrix;
        gradXX = @(S) S*gradXXmatrix;
        gradYY = @(S) gradYYmatrix*S ;
        lap = @(S) S*gradXXmatrix +  gradYYmatrix*S ;
        simsize = [Ny, Nx];
    elseif strcmp(solvermethod,'2Dlin')
        lapmatrix = laplacian_2Dlin([Nx,Ny],[dx dy],laplacianmethod,BCs);
        gradXmatrix = create_cdiff_lin([Nx,Ny],[dx dy],'x',BCs);
        gradYmatrix = create_cdiff_lin([Nx,Ny],[dx dy],'y',BCs);
        gradXYmatrix = create_cdiff_lin([Nx,Ny],[dx dy],'xy',BCs);
        gradXXmatrix = create_cdiff_lin([Nx,Ny],[dx dy],'xx',BCs);
        gradYYmatrix = create_cdiff_lin([Nx,Ny],[dx dy],'yy',BCs);
        gradX = @(S) gradXmatrix*S;
        gradY = @(S) gradYmatrix*S;
        gradYX = @(S) gradXYmatrix*S;
        gradXX = @(S) gradXXmatrix*S;
        gradYY = @(S) gradYYmatrix*S ;
        lap = @(S) lapmatrix*S;
        simsize = [Nx*Ny, 1];
    end % if solvermethod
end % func calcDifferentialOperators