%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Three strategies of parameters assignment
% in multi phase field model of grain growth with anisotorpic grain boundary properties”, Mendeley Data, 
% v1 http://dx.doi.org/10.17632/5wrv3ky9pp.1>, coupled to publication of the same name by 
% Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% 
%% cdiff = create_cdiff(systemsize,spacing,direction) 
% for x direction derivative systemsize is the number of columns and the
% matrix is applied as S*cdiff
% for y direction derivative systemsize is the number of rows and the
% matrix is applied as cdiff*S

% function cdiff = create_cdiff(systemsize,spacing,direction,order)
function cdiff = create_cdiff(systemsize,spacing,direction)
M = systemsize;

% switch order % order of derivative
%     case 1
        diag_plus1 = -ones(M,1);
        diag_minus1= ones(M,1);
        cdiff = spdiags([diag_plus1,diag_minus1],[1 -1],M,M)/(2*spacing);

% NOT PERIODIC - BACKWARD and FORWARD difference at the boundary
%         cdiff(1,1) = -1/spacing;
%         cdiff(2,1) = 1/spacing;
%         cdiff(M-1,M) = -1/spacing;
%         cdiff(M,M) = 1/spacing;

% PERIODIC - BACKWARD and FORWARD diffrence at boundary, 
%         cdiff(1,1) = -1/spacing;
%         cdiff(1,M) = 1/spacing;
%         cdiff(M,1) = 1/spacing;
%         cdiff(M,M) = -1/spacing;
%         cdiff(2,1) = 0;
%         cdiff(M-1,M) = 0;
        
% PERIODIC - CENTERED at boundary 
        cdiff(1,end) = 1/(2*spacing);
        cdiff(end,1) = -1/(2*spacing);
        
        if direction == 'y'
            cdiff = cdiff';
        end
%     case 2 % 2nd order derivative in  given direction, for laplacian
%         %second-order central difference, see wiki https://en.wikipedia.org/wiki/Finite_difference#Higher-order_differences
%         diag_plus1 = ones(M,1);
%         diag_minus1 = ones(M,1);
%         diag_  = -2*ones(M,1);
% end

end
    
% test
% dd = 0.1;
% [X,Y] = meshgrid(-5:dd:49*dd); % 
% X([1:9],:) = [];
% Y([1:9],:) = [];
% cdiffX = create_cdiff(size(X,2),X(1,2)-X(1,1),'x');
% cdiffY = create_cdiff(size(Y,1),Y(2,1)-Y(1,1),'y');
% Z1 = exp(-(X.^2)/10);
% Z1 = sin(-X*2*pi/(X(1,end)-X(1,1)));
% Z1 = sin(-X*2*pi/10);
% Z1g = Z1*cdiffX;
% pcolor(Z1), colorbar,daspect([1 1 1]),shading flat
% plot(X(1,:),Z1(1,:),'o',X(1,:),Z1g(1,:),'o'),grid on
% pcolor(Z1*cdiffX), colorbar,daspect([1 1 1])
% pcolor(Y*cdiffX), colorbar,daspect([1 1 1])
% pcolor(cdiffY*Y), colorbar,daspect([1 1 1])
% pcolor(cdiffY*X), colorbar,daspect([1 1 1])
% pcolor(X), colorbar,daspect([1 1 1]),shading flat
% pcolor(Y), colorbar,daspect([1 1 1])
% imagesc(cdiffX)
% imagesc(cdiffY)

% compare to centraldiff => matrix multiplication is 2.4 times faster
% [X,Y] = meshgrid(1:50);
% X([1:20],:) = [];
% Y([1:20],:) = [];
% % cdiffX = create_cdiff(size(X,2),1,'x');
% % cdiffY = create_cdiff(size(Y,1),1,'y');
% for i = 1:5000
%     tic 
%     A1 = X*cdiffX;
%     t1(i) = toc;
%     tic
%     A2 = centraldiff(X,[1 1],1);
%     t2(i) = toc;
%     [t1 t2];
% end
% plot(1:1000,t1,'b.',1:1000,t2,'r.'), legend('matrix multiplication','for loop')
% plot(t1./t2, '.'),title('time matrix multip/time for loop')
% mean(t2./t1) 
% pcolor(AA), colorbar,daspect([1 1 1]),shading interp
% all(all(A1==A2));

% compare to built-in gradient function - that is 3times faster
% dx = 0.02;
% dy = 0.05;
% [X,Y] = meshgrid(-2:dx:2,-3:dy:3);
% % % % cdiffX = create_cdiff(size(X,2),dx,'x');
% % % % cdiffY = create_cdiff(size(Y,1),dy,'y');
% Z = peaks(X,Y);
% pcolor(X,Y,Z), colorbar,daspect([1 1 1]),shading interp
% for i = 1:5000
%     tic 
%     dZdx = Z*cdiffX;
%     dZdy = cdiffY*Z;
%     t1(i) = toc;
%     tic
%     [dZdx2 dZdy2] = gradient(Z,dx,dy);
%     t2(i) = toc;
% end
% [t1 t2]
% pcolor(X,Y,Z), colorbar,daspect([1 1 1]),shading interp
% NN = size(dZdx);
% sum(abs(dZdx-dZdx2)<1e-13)/NN(1)
% 
% mean(t2./t1)