% author: Nele Moelans
% modifications: Martin Minar, 2020
function [kappa, gamma, m, L, thickness, gsq] = parameters_varIW(sigma,mu,thickness_0,sigma_0)

% kappa, gamma, m, L : model parameters as defined in
% Phys. Rev. B, 78 (2), 024113, July 2008.
%But we assume a variable grain boundary width and keep kappa and m fixed
%for all grain boundaries

%input
% sigma : vector with the grb energies, e.g. [0.25 0.14] (the sigma_min and sigma_max for ER=1.78 and ER=0.566)
% mu : vector with grb mobilities; e.g. [2.25e-12 2.25e-12]
% thickness_0 : the width of the boundary with energy sigma_0; for the
% paper: thickness_0=3e-7
% sigma_0 :  grain boundary energy of the boundary for which gamma is taken
% equal to 1.5 ; for the paper: sigma_0=0.25

%%outcome
% kappa : value of the energy gradient coefficient
% gamma : vector with different gamma values for the different boundaries
% with energies as specified in the input sigma
% m : coefficient in the homogeneous part of the interfacial energy
% L : vector with the kinetic coefficients, so that the grain boundary
% mobilities of the different boundaries are as specified in the vector mu
% in the input
% thickness : vector with the thickness of the different boundaries

% % probably fitted for 0.5 < gamma < 4 but works ok for gamma < 40
% gaminv_gsq_coef = [-5.288 , -0.09364 , 9.965 , -8.183 ,2.007 ];
% % below: giving wrong results for gamma < 0.75, see visualization
% f0_gaminv_coef = [0.05676 , -0.2924 , 0.6367  ,-0.7749  ,  0.6107  , -0.4324 ,    0.2792]; 

% both fitted for 0.52 <= gamma <= 40
gaminv_gsq_coef = [-5.73008      18.8615     -23.0557      7.47952      8.33568     -8.01224      2.00013];
sqrtf0_gaminv_coef = [-0.072966     0.35784    -0.68325     0.63578    -0.48566     0.53703];

IWmin = thickness_0;
IWiter = 0;
% iteration to find such gammas to give correct IE AND IWmin==thickness0
for ctr = 0:1
    if ctr == 0
        IWiter = IWmin;
    end
    kappa = 3/4 *sigma_0*IWiter ;
    m = 6*sigma_0/IWiter;

    gsq = sigma.^2/(m*kappa);

    gamma_inverse =  polyval(gaminv_gsq_coef,gsq);
    sqrtf_interf =  polyval(sqrtf0_gaminv_coef,gamma_inverse);
    thickness = sqrt(kappa./m)./sqrtf_interf;
    
    IWiter = IWmin*sqrt(8)*max(unique(sqrtf_interf)); % IW_0 = IWmin*sqrt(f_0(gam_max)/f_0(1.5)) = IWmin*sqrt(8*f_0(gam_max))
end

gamma = 1./gamma_inverse;
L = mu.*sigma/kappa;

% f_interf =  polyval(f0_gaminv_coef,gamma_inverse);
% thickness = sqrt(kappa./(m*f_interf));
 
end % func


%% Nele's version compared to mine - passed
% testpar = inputfile_GGaniso;
% testpar_old = inputfile_GGaniso;
% .
% [testpar.kpp0 == testpar_old.kpp0 ;
% testpar.gam0 == testpar_old.gam0 ;
% testpar.IWs == testpar_old.IWs ;
% testpar.m == testpar_old.m ;
% testpar.Lij == testpar_old.Lij ;
% testpar.soaKP == testpar_old.soaKP]


%% fits visualization, compare g(\gamma) to fig 8 in PRB2008
% clear
% gg_dlim = 0.05;
% gg_uplim = 0.7385; % choose between 0.4 - 0.75
% % gg MUST NOT be above 0.75, ggsq MUST NOT  be above 0.5625 
% q = [-5.288 , -0.09364 , 9.965 , -8.183 ,2.007 ];
% ggsq = linspace(gg_dlim,gg_uplim,200).^2;
% % ggsq = linspace(sqrt(0.05),0.57565,100);
% % plot(ggsq,polyval(q,ggsq),'.'), grid on
% plotstyle = {'-','linewidth',1.5,'markersize',6}
% figure(66)
% plot(1./polyval(q,ggsq),sqrt(ggsq),plotstyle{:})
% xlabel('\gamma'),
% hold on
% 
% qq = [0.05676 , -0.2924 , 0.6367  ,-0.7749  ,  0.6107  , -0.4324 ,    0.2792];
% gamma = 1./polyval(q,ggsq);
% plot(gamma,4/3*sqrt(polyval(qq,1./gamma)),plotstyle{:})
% hold off
% legend('g(\gamma)','4/3( f_{0,i}(\gamma) )^{1/2}','location','best')
% set(gca,'fontsize',13)
% xticks([0:.5:2 ,3:10])
% % ylim([0.1, 0.8])
% xlim([0, 10])
% grid on
% % 
% dest = 'c:\Users\marti\PhD\MySimMatlab\MultiPAniso\figs\'
% save_figure_fig_png(66,[dest 'param_determ_polynonials_g_and_f0'])

%% Nele's way
%%%%%%%fitting for gamma < 4 (gamma(3:67))
% %%% g_function^2 -- 1/gamma
% 1/gamma =  p1*(g_function^2)^3 + p2*(g_function^2)^2 + p3*(g_function^2)^1 + p4;

% %% 1/gamma -- f_interf 
% f_interf = pp1*(1/gamma)^6 + pp2*(1/gamma)^5 + pp3*(1/gamma)^4 + pp4*(1/gamma)^3 + pp5*(1/gamma)^2 + pp6*(1/gamma)^1 + pp7;

  %%for gamma < 4,fitted for gamma(3:67)
%   p1 =      -5.288 ; 
%   p2 =    -0.09364;  
%   p3 =       9.965  ;
%   p4 =      -8.183  ;
%   p5 =       2.007  ;
%        
%        pp1 =     0.05676;  
%        pp2 =     -0.2924;  
%        pp3 =      0.6367  ;
%        pp4 =     -0.7749  ;
%        pp5 =      0.6107  ;
%        pp6 =     -0.4324 ; 
%        pp7 =      0.2792  ;

% kappa = 3/4 *sigma_0*thickness_0 ;
% m = 6*sigma_0/(thickness_0);
% % sigma is a vector
% gsq = sigma.^2/(m*kappa);

%Nele's way
%  for k = 1 : length(sigma)
%     gamma_inverse =  p1*gsq(k)^4 + p2*gsq(k)^3 + p3*gsq(k)^2 + p4*gsq(k) + p5;
%     gamma(k) = 1/gamma_inverse;
%     L(k) = (mu(k)*sigma(k))/kappa;
%     f_interf =   pp1*gamma_inverse^6 + pp2*gamma_inverse^5 + pp3*gamma_inverse^4 ...
%      + pp4*gamma_inverse^3 + pp5*gamma_inverse^2 + pp6*gamma_inverse^1 + pp7;
%     thickness(k) = sqrt(kappa/(m*f_interf));  
%  end

%% mine, works
% kappa = 3/4 *sigma_0*thickness_0 ;
% m = 6*sigma_0/thickness_0;
% 
% gsq = sigma.^2/(m*kappa);
% 
% % % probably fitted for 0.5 < gamma < 4 but works ok for gamma < 40
% % gaminv_gsq_coef = [-5.288 , -0.09364 , 9.965 , -8.183 ,2.007 ];
% % % below: giving wrong results for gamma < 0.75, see visualization
% % f0_gaminv_coef = [0.05676 , -0.2924 , 0.6367  ,-0.7749  ,  0.6107  , -0.4324 ,    0.2792]; 
% 
% % both fitted for 0.52 <= gamma <= 40
% gaminv_gsq_coef = [-5.73008      18.8615     -23.0557      7.47952      8.33568     -8.01224      2.00013];
% sqrtf0_gaminv_coef = [-0.072966     0.35784    -0.68325     0.63578    -0.48566     0.53703];
% 
% gamma_inverse =  polyval(gaminv_gsq_coef,gsq);
% gamma = 1./gamma_inverse;
% L = mu.*sigma/kappa;
% sqrtf_interf =  polyval(sqrtf0_gaminv_coef,gamma_inverse);
% thickness = sqrt(kappa./m)./sqrtf_interf;
% % f_interf =  polyval(f0_gaminv_coef,gamma_inverse);
% % thickness = sqrt(kappa./(m*f_interf));
