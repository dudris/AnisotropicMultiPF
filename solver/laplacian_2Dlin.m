%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Three strategies of parameters assignment
% in multi phase field model of grain growth with anisotorpic grain boundary properties”, Mendeley Data, 
% v1 http://dx.doi.org/10.17632/5wrv3ky9pp.1>, coupled to publication of the same name by 
% Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754)
% 
%% lap = laplacian_2Dlin(systemsize,spacing,StencilType,BCs)
% - systemsize = [Nx Ny], where
%         Nx is number of gridpoints in x-dir (i.e. columns in system matrix)
%         Ny is number of gridpoints in y-dir (i.e. rows in system matrix)
% - spacing = [dx dy]
% - StencilType switch options:  '5pt', 9pt8, '9pt20'
% - BCs switch options: 'P', 'N'


function lap = laplacian_2Dlin(systemsize,spacing,StencilType,BCs)
% error('must check the prefactors in laplacians')
Ny = systemsize(2);
Nx = systemsize(1);
h =  spacing(1);

assert(spacing(1)==spacing(2),'laplacian_2Dlin error: Spacing dx~=dy')

switch StencilType
    case '5pt'
        stencil_coef = h^2;
        switch BCs
            case 'P'
                M1 = sparse(toeplitz([-4, 1, zeros(1,Ny-3),1]));
                M2 = sparse(eye(Ny));
                M3 = sparse(zeros(Ny));
                M = {M1 , M2 , M3};
        %         imagesc(M1),colorbar
                matrixGrid = toeplitz([1,2,3*ones(1,Nx-3),2]);
                lap = cell2mat(M(matrixGrid))/stencil_coef;
%                 imagesc(lap),colorbar

            case 'N'
                M1 = sparse(toeplitz([-4, 1, zeros(1,Ny-2)]));
                M1(1,[1 2]) = [-4 2] ;
                M1(end,[Ny-1 Ny]) = [2 -4] ;
%                 imagesc(M1),colorbar
                M2 = sparse(eye(Ny));
                M2_2 = 2*M2;
                M3 = sparse(zeros(Ny));
                M = {M1 , M2, M2_2 , M3};
                
                matrixGrid = toeplitz([1,2,4*ones(1,Nx-2)]);
                matrixGrid(1,2) = 3 ;
                matrixGrid(end,Nx-1) = 3 ;
                lap = cell2mat(M(matrixGrid))/stencil_coef;
%                 imagesc(lap),colorbar, daspect([1 1 1])
%                 sum(cdiff,2)

            case 'mix' % periodic in x, Neumann in y
                M1_N = sparse(toeplitz([-4, 1, zeros(1,Ny-2)]));
                M1_N(1,[1 2]) = [-4 2] ;
                M1_N(end,[Ny-1 Ny]) = [2 -4] ;
                M2_P = sparse(eye(Ny));
                M3_P = sparse(zeros(Ny));
                M = {M1_N , M2_P , M3_P};
        %         imagesc(M1),colorbar
                matrixGrid = toeplitz([1,2,3*ones(1,Nx-3),2]);
                lap = cell2mat(M(matrixGrid))/stencil_coef;
%                 imagesc(lap),colorbar, daspect([1 1 1])
%                 sum(lap,2)
                
        end% switch BCs
        
    case '9pt8'
        stencil_coef = 3*h^2;
        switch BCs
            case 'P'
                M1 = sparse(toeplitz([-8, 1, zeros(1,Ny-3),1]));
        %         imagesc(M1),colorbar
                M2 = sparse(toeplitz([1, 1, zeros(1,Ny-3),1]));
        %         imagesc(M2),colorbar     
                M3 = sparse(zeros(Ny));
                M = {M1 , M2 , M3};
        %         imagesc(M1),colorbar
                matrixGrid = toeplitz([1,2,3*ones(1,Nx-3),2]);
                lap = cell2mat(M(matrixGrid))/stencil_coef; % 
    %                 imagesc(M1),colorbar,daspect([1 1 1])
%                     imagesc(lap),colorbar,daspect([1 1 1])
            case 'N'
                 M1 = sparse(toeplitz([-8, 1, zeros(1,Ny-2)]));
                M1(1,[1 2]) = [-8 2] ;
                M1(end,[Ny-1 Ny]) = [2 -8] ;
%                 imagesc(M1),colorbar, daspect([1 1 1])
                M2 = sparse(toeplitz([1, 1, zeros(1,Ny-2)]));
                M2(1,[1 2]) = [1 2] ;
                M2(end,[Ny-1 Ny]) = [2 1] ;
%                 imagesc(M2),colorbar, daspect([1 1 1])
                M2_2 = 2*M2;
                M3 = sparse(zeros(Ny));
                M = {M1 , M2, M2_2 , M3};
                
                matrixGrid = toeplitz([1,2,4*ones(1,Nx-2)]);
                matrixGrid(1,2) = 3 ;
                matrixGrid(end,Nx-1) = 3 ;
                lap = cell2mat(M(matrixGrid))/stencil_coef;
%                 imagesc(lap*3),colorbar, daspect([1 1 1])
%                 sum(lap*3,2)

            case 'mix' % periodic in x, Neumann in y
                M1_N = sparse(toeplitz([-8, 1, zeros(1,Ny-2)]));
                M1_N(1,[1 2]) = [-8 2] ;
                M1_N(end,[Ny-1 Ny]) = [2 -8] ;
        %         imagesc(M1_N),colorbar
                M2 = sparse(toeplitz([1, 1, zeros(1,Ny-2)]));
                M2(1,[1 2]) = [1 2] ;
                M2(end,[Ny-1 Ny]) = [2 1] ;
        %         imagesc(M2_P),colorbar     
                M3_P = sparse(zeros(Ny));
                M = {M1_N , M2 , M3_P};
        %         imagesc(M1),colorbar
                matrixGrid = toeplitz([1,2,3*ones(1,Nx-3),2]);
                lap = cell2mat(M(matrixGrid))/stencil_coef; % 
%                 imagesc(lap*3),colorbar, daspect([1 1 1])
%                 sum(lap*3,2)

        end % switch BCs

    case '9pt20'
        stencil_coef = 6*h^2;
        
        switch BCs
            case 'P'
                M1 = sparse(toeplitz([-20, 4, zeros(1,Ny-3),4]));
        %         imagesc(M1),colorbar
                M2 = sparse(toeplitz([4, 1, zeros(1,Ny-3),1]));
        %         imagesc(M2),colorbar     
                M3 = sparse(zeros(Ny));
                M = {M1 , M2 , M3};
        %         imagesc(M1),colorbar
                matrixGrid = toeplitz([1,2,3*ones(1,Nx-3),2]);
                lap = cell2mat(M(matrixGrid))/stencil_coef;
%                     imagesc(lap*6),colorbar,daspect([1 1 1])
%                     sum(lap*6,2)
            case 'N'
                M1 = sparse(toeplitz([-20, 4, zeros(1,Ny-2)]));
                M1(1,[1 2]) = [-20 8] ;
                M1(end,[Ny-1 Ny]) = [8 -20] ;
%                 imagesc(M1),colorbar, daspect([1 1 1])
                M2 = toeplitz([4, 1, zeros(1,Ny-2)]);
                M2(1,[1 2]) = [4 2] ;
                M2(end,[Ny-1 Ny]) = [2 4] ;
%                 imagesc(M2),colorbar, daspect([1 1 1])
                M2_2 = 2*M2;
                M3 = sparse(zeros(Ny));
                M = {M1 , M2, M2_2 , M3};
                
                matrixGrid = sparse(toeplitz([1,2,4*ones(1,Nx-2)]));
                matrixGrid(1,2) = 3 ;
                matrixGrid(end,Nx-1) = 3 ;
                lap = cell2mat(M(matrixGrid))/stencil_coef;
%                 imagesc(lap*6),colorbar,daspect([1 1 1])
%                 sum(lap*6,2)
            case 'mix'
                M1_N = sparse(toeplitz([-20, 4, zeros(1,Ny-2)]));
                M1_N(1,[1 2]) = [-20 8] ;
                M1_N(end,[Ny-1 Ny]) = [8 -20] ;
        %         imagesc(M1_N),colorbar
                M3_P = sparse(zeros(Ny));
                
%                 M2_P = sparse(toeplitz([4, 1, zeros(1,Ny-3),1]));
%                 M = {M1_N , M2_P , M3_P};
                
                M2 = sparse(toeplitz([4, 1, zeros(1,Ny-2)]));
                M2(1,[1 2]) = [4, 2] ;
                M2(end,[Ny-1 Ny]) = [2 , 4] ;
        %         imagesc(M2_P),colorbar     
                M = {M1_N , M2 , M3_P};
                
        %         imagesc(M1),colorbar
                matrixGrid = toeplitz([1,2,3*ones(1,Nx-3),2]);
                lap = cell2mat(M(matrixGrid))/stencil_coef;
%                 imagesc(lap*6),colorbar,daspect([1 1 1]),grid on
%                 sum(lap*6,2)
        end % switch BCs
        
end % switch StencilType
        

end % func


% test mixed BCs

