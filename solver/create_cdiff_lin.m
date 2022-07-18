%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Benchmarking of different strategies to include anisotropy in a curvature-driven multi-phase-field model”, Mendeley Data, V2, doi: 10.17632/5wrv3ky9pp.2> 
% coupled to publication of the same name by Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754).
% 
%% cdiff = create_cdiff_lin(systemsize,spacing,direction,BCs)
% ... creates finite differences matrix cdiff, size(cdiff) = [Nx*Ny Nx*Ny]
% ... cdiff stands for partial derivative operator in x, y direction and also
% for 2nd order derivatives, i.e. xx, xy and yy
% ... periodic and Neuman and mixed BCs possible (mixed... P @ LR sides and N @ TB sides)
% - systemsize = [Nx Ny], where
%         Nx is number of gridpoints in x-dir (i.e. columns in system matrix)
%         Ny is number of gridpoints in y-dir (i.e. rows in system matrix)
% - spacing = [dx dy]
% - direction switch options:  'x', 'y', 'xx', 'yy', 'xy'
% - BCs switch options: 'P', 'N', 'mix'

function cdiff = create_cdiff_lin(systemsize,spacing,direction,BCs)

% systemsize = [4,4];
% spacing = [1 1];
% BCs = 'N';
% direction = 'xy';
Ny = systemsize(2);
Nx = systemsize(1);
dx = spacing(1);
dy = spacing(2);

switch BCs
    case 'P' %periodic

    if strcmp(direction,'x')
        matrixGrid = sparse(toeplitz([0, -1, zeros(1,Nx-3),1],[0,1, zeros(1,Nx-3),-1]));
        M1 = speye(Ny);
        cdiff = kron(matrixGrid,M1)/(2*dx);

    elseif strcmp(direction,'y')
        matrixGrid = speye(Nx);
        M1 = sparse(toeplitz([0, -1, zeros(1,Ny-3),1],[0,1, zeros(1,Ny-3),-1]));
        cdiff = kron(matrixGrid,M1)/(2*dy);
        
    elseif strcmp(direction,'xx')
        matrixGrid = sparse(toeplitz([-2, 1, zeros(1,Nx-3),1]));
        M1 = speye(Ny);
        cdiff = kron(matrixGrid,M1)/(dx*dx);
%         imagesc(cdiff),grid on,colorbar
    
    elseif strcmp(direction,'yy') 
        matrixGrid = speye(Nx);
        M1 = sparse(toeplitz([-2, 1, zeros(1,Ny-3),1]));
        cdiff = kron(matrixGrid,M1)/(dy*dy);
%         imagesc(cdiff),grid on,colorbar
    
    elseif strcmp(direction,'xy') 
        matrixGrid = sparse(toeplitz([0, -1, zeros(1,Nx-3),1],[0,1, zeros(1,Nx-3),-1]));
        M1 = speye(Ny);
        cdiff1 = kron(matrixGrid,M1)/(2*dx);
        %
        matrixGrid = speye(Nx);
        M1 = sparse(toeplitz([0, -1, zeros(1,Ny-3),1],[0,1, zeros(1,Ny-3),-1]));
        cdiff2 = kron(matrixGrid,M1)/(2*dy);
        %
        cdiff = cdiff1 +cdiff2;
%         imagesc(cdiff),grid on,colorbar
    else
        error('Invalid switch ''direction'' used. Must be either ''x'',''y'',''xx'',''yy'' or ''xy''')
    end % if direction
    
    case 'N' %neumann 
        if strcmp(direction,'x')
            matrixGrid = sparse(toeplitz([0, -1, zeros(1,Nx-3),1],[0,1, zeros(1,Nx-3),-1]));
            matrixGrid([1,Nx],:) = 0;
            M1 = eye(Ny);
            cdiff = kron(matrixGrid,M1)/(2*dx);
        
        elseif strcmp(direction,'y')
            matrixGrid = eye(Nx);
            M1 = sparse(toeplitz([0, -1, zeros(1,Ny-3),1],[0,1, zeros(1,Ny-3),-1]));
            M1([1,Ny],:) = 0;
            cdiff = kron(matrixGrid,M1)/(2*dy);
%             imagesc(cdiff),grid on,daspect([1 1 1]),colorbar
%             spy(cdiff),grid on,daspect([1 1 1]),colorbar
    
        elseif strcmp(direction,'xx')
            matrixGrid = sparse(toeplitz([-2, 1, zeros(1,Nx-2)]));
            matrixGrid(1,2) = 2;
            matrixGrid(end,end-1) = 2;
            M1 = speye(Ny);
            cdiff = kron(matrixGrid,M1)/(dx*dx);
%             imagesc(cdiff),grid on,colorbar
    
        elseif strcmp(direction,'yy')
            matrixGrid = speye(Nx);
            M1 = sparse(toeplitz([-2, 1, zeros(1,Ny-2)]));
            M1(1,2) = 2; 
            M1(end,end-1) = 2;
            cdiff = kron(matrixGrid,M1)/(dy*dy);
%         imagesc(cdiff),grid on,colorbar
        
        elseif strcmp(direction,'xy') 
            matrixGrid = sparse(toeplitz([0, -1, zeros(1,Nx-3),1],[0,1, zeros(1,Nx-3),-1]));
            matrixGrid([1,Nx],:) = 0;
            M1 = eye(Ny);
            cdiff1 = kron(matrixGrid,M1)/(2*dx);
            %
            matrixGrid = eye(Nx);
            M1 = sparse(toeplitz([0, -1, zeros(1,Ny-3),1],[0,1, zeros(1,Ny-3),-1]));
            M1([1,Ny],:) = 0;
            cdiff2 = kron(matrixGrid,M1)/(2*dy);
            %
            cdiff = cdiff1 +cdiff2;
%             imagesc(cdiff),grid on,colorbar
        else
            error('Invalid switch ''direction'' used. Must be either ''x'',''y'',''xx'',''yy'' or ''xy''')
        end % if direction
        
    case 'mix' % mixed NBC and PBC
        % matrixGrid from PBC and M1 from NBC
        if strcmp(direction,'x')
            matrixGrid = sparse(toeplitz([0, -1, zeros(1,Nx-3),1],[0,1, zeros(1,Nx-3),-1]));
            M1 = eye(Ny);
            cdiff = kron(matrixGrid,M1)/(2*dx);
        
        elseif strcmp(direction,'y')
            matrixGrid = speye(Nx);
            M1 = sparse(toeplitz([0, -1, zeros(1,Ny-3),1],[0,1, zeros(1,Ny-3),-1]));
            M1([1,Ny],:) = 0;
            cdiff = kron(matrixGrid,M1)/(2*dy);
% %             imagesc(cdiff),grid on,daspect([1 1 1]),colorbar
% %             spy(cdiff),grid on,daspect([1 1 1]),colorbar
%     
        elseif strcmp(direction,'xx')
            matrixGrid = sparse(toeplitz([-2, 1, zeros(1,Nx-3),1]));
            M1 = speye(Ny);
            cdiff = kron(matrixGrid,M1)/(dx*dx);
% %             imagesc(cdiff),grid on,colorbar
%     
        elseif strcmp(direction,'yy')
            matrixGrid = speye(Nx);
            M1 = sparse(toeplitz([-2, 1, zeros(1,Ny-2)]));
            M1(1,2) = 2; 
            M1(end,end-1) = 2;
            cdiff = kron(matrixGrid,M1)/(dy*dy);
% %         imagesc(cdiff),grid on,colorbar
%         
        elseif strcmp(direction,'xy') 
            matrixGrid = sparse(toeplitz([0, -1, zeros(1,Nx-3),1],[0,1, zeros(1,Nx-3),-1]));
            M1 = eye(Ny);
            cdiff1 = kron(matrixGrid,M1)/(2*dx);
            %
            matrixGrid = speye(Nx);
            M1 = sparse(toeplitz([0, -1, zeros(1,Ny-3),1],[0,1, zeros(1,Ny-3),-1]));
            M1([1,Ny],:) = 0;
            cdiff2 = kron(matrixGrid,M1)/(2*dy);
            %
            cdiff = cdiff1 +cdiff2;
% %             imagesc(cdiff),grid on,colorbar
        else
            error('Invalid switch ''direction'' used. Must be either ''x'',''y'',''xx'',''yy'' or ''xy''')
        end % if direction
        
end% switch BCs
        
end % function