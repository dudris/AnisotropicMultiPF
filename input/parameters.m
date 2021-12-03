% author: Nele Moelans, 2007,
% modifications: Martin Minar
function [kappa, gam,gsq, m, L] = parameters(sigma,mu,thickness,sigma_0)

% kappa, gamma, m, L : model parameters as defined in
% Phys. Rev. B, 78 (2), 024113, July 2008.ttt
% sigma : vector with grb energies, e.g. [0.1 0.2 0.25]
% m : vector with grb mobilities; e.g. [1 1 1]
% thickness : numerical grb thickness, e.g. 1e-7
% sigma_0 :  start_value for the iterative calculation, gamma will equal
% 1.5 for sigma = sigma_0 (it is best to choose sigma_0 so that gamma is between 0.7 and 4 for all elements of sigma)

% example :  [kappa, gamma, m, L] = parameters([0.1 0.2 0.3],[1 1 1],1e-8,0.25) 
% gives kappa =  1.0e-007 * [  0.0739    0.1497    0.1875 ]
% gamma = [0.7412    1.1629    1.4999]
% m =   1.50e+007
% L =  1.0e+007 * [ 1.3538    1.3356    1.3334 ]

%% former Nele's
%  %%for gamma < 4,fitted for gamma(3:67)
% p1 =    -5.288;  
% p2 =    -0.09364;  
% p3 =     9.965;  
% p4 =    -8.183;  
% p5 =     2.007;  
% 
% %for calculation of g(gamma)^2
% gsq_range = 0.05:0.01:0.45;
% gamma_inversed = polyval([p1 p2 p3 p4 p5],gsq_range); % 1/gamma
% 
% pp1 =    0.05676;  
% pp2 =   -0.2924;  
% pp3 =    0.6367;  
% pp4 =   -0.7749;  
% pp5 =    0.6107;  
% pp6 =   -0.4324;  
% pp7 =    0.2792;  
%   
% %calculate double well height 
% % holds only for gamma_init=1.5, because gfuncinit*sqrt(f0int)=(4/3)*f0interf = (4/3)*(1/8) = 1/6
% % but for other value g(gamma_init) and f_0interf(gamma_init)) would have
% % to be calculated
% m = 6*sigma_0/(thickness); 
% 
% g_0 = sqrt(2)/3;
% sqrt_f0 = sqrt(1/8);
% a_0 = sqrt_f0/g_0;
%  for k = 1 : length(sigma)
%      a_star = a_0;
%      a_0 = 0;
%      while abs(a_0 - a_star) > 1e-9
%          a_0 = a_star;
%          kappa_star = sigma(k)*thickness*a_0;% calculate kappa
%          %Numerically calculated g
%          g = sigma(k)/sqrt(kappa_star*m);
%          x = g^2;
%          y = p1*x^4 + p2*x^3 + p3*x^2 + p4*x + p5;
%          %why?
%          gamma_star = 1/y;
%          f_interf =   pp1*y^6 + pp2*y^5 + pp3*y^4 + pp4*y^3 + pp5*y^2 + pp6*y^1 + pp7;
%          a_star = sqrt(f_interf)/g;
%      end % while
%      kappa(k) = kappa_star;
%      gam(k) = gamma_star;  
%      L(k) = (mu(k)/(thickness))*(1/a_star);
%      
%      gsq(k) = interp1(1./gamma_inversed,gsq_range,gam(k),'spline');
%  end% for


%% mine
% both fitted for 0.52 <= gamma <= 40
gaminv_gsq_coef = [-5.73008      18.8615     -23.0557      7.47952      8.33568     -8.01224      2.00013];
sqrtf0_gaminv_coef = [-0.072966     0.35784    -0.68325     0.63578    -0.48566     0.53703];

%for calculation of g(gamma)^2
gsq_range = linspace(.098546^2,0.765691^2,100);
gamma_inversed = polyval(gaminv_gsq_coef,gsq_range); % 1/gamma
% gamma_inversed = polyval([p1 p2 p3 p4 p5],gsq_range); % 1/gamma
  
%calculate double well height 
% holds only for gamma_init=1.5, because gfuncinit*sqrt(f0int)=(4/3)*f0interf = (4/3)*(1/8) = 1/6
% but for other value g(gamma_init) and f_0interf(gamma_init)) would have
% to be calculated
m = 6*sigma_0/(thickness); 

g_0 = sqrt(2)/3;
sqrt_f0 = sqrt(1/8);
a_0 = sqrt_f0/g_0;
 for k = 1 : length(sigma)
     a_star = a_0;
     a_0 = 0;
     ctr = 0;
     while abs(a_0 - a_star) > 1e-9 
         a_0 = a_star;
         kappa_star = sigma(k)*thickness*a_0;% calculate kappa
         %Numerically calculated g
         g = sigma(k)/sqrt(kappa_star*m);
         x = g^2;
         y = polyval(gaminv_gsq_coef,x);
         %why?
         gamma_star = 1/y;
         sqrtf_interf =   polyval(sqrtf0_gaminv_coef,y);
         % in Nele's original the polynomial is for f_0i, I have it for sqrt(f_0i)
%          a_star = sqrt(f_interf)/g;
         a_star = sqrtf_interf/g;
         ctr = ctr+1;
         if ctr==60
             error('Parameters determination algorithm not converging.')
         end
     end % while
     kappa(k) = kappa_star;
     gam(k) = gamma_star;  
     L(k) = (mu(k)/(thickness))*(1/a_star);
     
 end% for
 
 gsq = interp1(1./gamma_inversed,gsq_range,gam);


 end % func
 
 %% Nele's original
 
% %%%%%%%fitting for gamma < 4 (gamma(3:67))
% % %%% g_function^2 -- 1/gamma
% % 1/gamma =  p1*(g_function^2)^3 + p2*(g_function^2)^2 + p3*(g_function^2)^1 + p4;
% %   p1 = -4.8457;
% %   p2 = 11.3;
% %   p3 = -8.2937;
% %   p4 = 2.006;
% 
% 
% % %% 1/gamma -- f_interf 
% % f_interf = pp1*(1/gamma)^6 + pp2*(1/gamma)^5 + pp3*(1/gamma)^4 + pp4*(1/gamma)^3 + pp5*(1/gamma)^2 + pp6*(1/gamma)^1 + pp7;
% %  pp1 = 0.028948
% %   pp2 = -0.19627
% %   pp3 = 0.53653
% %   pp4 = -0.7693
% %   pp5 = 0.65764
% %   pp6 = -0.45665
% %   pp7 = 0.28256
 
%  %%for gamma < 4,fitted for gamma(3:67)
% p1 =    -5.288;  
% p2 =    -0.09364;  
% p3 =     9.965;  
% p4 =    -8.183;  
% p5 =     2.007;  
% 
% %for calculation of g(gamma)^2
% gsq_range = 0.05:0.01:0.45;
% gamma_inversed = polyval([p1 p2 p3 p4 p5],gsq_range); % 1/gamma
% 
% pp1 =    0.05676;  
% pp2 =   -0.2924;  
% pp3 =    0.6367;  
% pp4 =   -0.7749;  
% pp5 =    0.6107;  
% pp6 =   -0.4324;  
% pp7 =    0.2792;  
%   
% %calculate double well height 
% % holds only for gamma_init=1.5, because gfuncinit*sqrt(f0int)=(4/3)*f0interf = (4/3)*(1/8) = 1/6
% % but for other value g(gamma_init) and f_0interf(gamma_init)) would have
% % to be calculated
% m = 6*sigma_0/(thickness); 
% 
% g_0 = sqrt(2)/3;
% sqrt_f0 = sqrt(1/8);
% a_0 = sqrt_f0/g_0;
%  for k = 1 : length(sigma)
%      a_star = a_0;
%      a_0 = 0;
%      while abs(a_0 - a_star) > 1e-9
%          a_0 = a_star;
%          kappa_star = sigma(k)*thickness*a_0;% calculate kappa
%          %Numerically calculated g
%          g = sigma(k)/sqrt(kappa_star*m);
%          x = g^2;
%          y = p1*x^4 + p2*x^3 + p3*x^2 + p4*x + p5;
%          %why?
%          gamma_star = 1/y;
%          f_interf =   pp1*y^6 + pp2*y^5 + pp3*y^4 + pp4*y^3 + pp5*y^2 + pp6*y^1 + pp7;
%          a_star = sqrt(f_interf)/g;
%      end % while
%      kappa(k) = kappa_star;
%      gam(k) = gamma_star;  
%      L(k) = (mu(k)/(thickness))*(1/a_star);
%      
%      gsq(k) = interp1(1./gamma_inversed,gsq_range,gam(k),'spline');
%  end% for