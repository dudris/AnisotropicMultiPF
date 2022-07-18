%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Benchmarking of different strategies to include anisotropy in a curvature-driven multi-phase-field model”, Mendeley Data, V2, doi: 10.17632/5wrv3ky9pp.2> 
% coupled to publication of the same name by Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754).
% 
%% laplacian1D
% lap = laplacean(Nx,dx, BCmethod)
% calculates matrix for numerical 1D finite-difference expression of laplace operator
% equidistant grid with Nx pts
% Neumann BC not tested as of 7/4/2020
% BCmethod... cell 1x2 with 'N' or 'D' to specify type of BC on the top {1}
% or bottom {2}
% Nx  ... number of elements
% dx  ... istance between grid points

function lap = laplacian1D(Nx,dx, BCmethod,direction)

r = zeros(Nx,1);
r([1,2]) = [-2,1];
lap = sparse(toeplitz(r)); % central difference in -y (vertical direction) including periodic boundaries

if BCmethod{1} == 'D' % DIRICHLET ; BCmethod{1} is UPPER BC
    lap(1,[1 2])= [0 0];
end
if BCmethod{2} == 'D' % DIRICHLET ; BCmethod{2} is LOWER BC
    lap(Nx,[Nx-1 Nx])= [0 0];
end
if BCmethod{1} == 'N'% NEUMANN ; BCmethod{1} is UPPER BC
    lap(1,[1 2])= [-2 2];
%     lap(2,[1 2 3])= [0 0 0];
end
if BCmethod{2} == 'N' % NEUMANN ; BCmethod{2} is LOWER BC
%     lap(Nx-1,[Nx-2 Nx-1 Nx])= [0 0 0];
    lap(Nx,[Nx-1 Nx])= [2 -2];
end
if BCmethod{1} == 'P'% periodic
    lap(1,end)= 1;
    lap(end,1)= 1;
%     lap(2,[1 2 3])= [0 0 0];
end

% if BCmethod{1} == 'D' % DIRICHLET ; BCmethod{1} is UPPER BC
%     lap(1,[1 2])= [1 0];
% end
% if BCmethod{2} == 'D' % DIRICHLET ; BCmethod{2} is LOWER BC
%     lap(Nx,[Nx-1 Nx])= [0 1];
% end
% if BCmethod{1} == 'N'% NEUMANN ; BCmethod{1} is UPPER BC
%     lap(1,[1 2])= [1 0];
%     lap(2,[1 2 3])= [1 0 0];
% end
% if BCmethod{2} == 'N' % NEUMANN ; BCmethod{2} is LOWER BC
%     lap(Nx-1,[Nx-2 Nx-1 Nx])= [0 0 1];
%     lap(Nx,[Nx-1 Nx])= [0 1];
% end

if strcmp(direction,'y')
    lap = lap/(dx*dx); 
elseif strcmp(direction,'x')
    lap = lap'/(dx*dx); % transposed
end
    
%
% imagesc(lap), colorbar
  %  spy(lap, '.'), grid on
end%func

%% test 2D
% counts nonzero elements in each row of lap <=> no of neighbours
% every line has 5 elements now 

%clear nonzero
%nonzero = 0;
%q = 0;
%for j = 1:size(lap2,1)
%  if nnz(lap2(j,:)) ~= 5
%    q = q+1;
%    nonzero(q,:) = [nnz(lap2(j,:)), j];
%  end
%end

%% test 1D
%Nx = 100;
%dx = 1;
%[~,~,A] = laplacian(Nx,{'P'}); % function from specfun package
%lap = laplacean_(Nx,dx, 1); % my function
%isequal(lap,-A); % WHY NEGATIVE OF THE OTHER?
