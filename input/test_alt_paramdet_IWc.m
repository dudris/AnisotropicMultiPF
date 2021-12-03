%% testing alternative IWc parameters determination strategy
% modify 'slope' between 0 and cca 24 to see the effect of IEinit  on gamma
% any other function valued between IEinit_min and IEinit_max could be used
% instead of the straight line
clear 

g_function	= [0.098546		0.765691]';
sqrt_f0 = [	0.07014984		0.5309548]';
Gf = g_function([1,end]).*sqrt_f0([1,end]);


% polynomial sqrt(f0)*g -> 1/gam; fitted for gamma 0.52-40
prodgsqrtf0c_gaminv_coef = [103.397      -165.393      105.3469     -44.55661       24.7348     -11.25718      1.999642];
% polynomial 1/gam -> sqrt(f0(1/gam)); fitted for gamma 0.52-40
sqrtf0_gaminv_coef = [-0.072966     0.35784    -0.68325     0.63578    -0.48566     0.53703];


ratios = linspace(0.017,1,300)';
IW = 1e-6;
IE_0 = 1;
IEs = IE_0*[ones(size(ratios)) ratios];
IEinit_min = max(IEs,[],2)./6/Gf(2); % 1st column
IEinit_max = min(IEs,[],2)./6/Gf(1); % 2nd column
% IEinit = mean(IEs,2);
% IEinit = IEinit_min;
% IEinit = IEinit_max*0.99;

crossval = (0.017./6/Gf(1));
slope = 0.5; % between 0-24, 
IEinit = crossval+slope*(ratios-0.017);
% IEinit = mean(IEs,2);

for kk=1:length(ratios)
    mus = 7.5e-16*ones(1,2);

    [kappa(kk,:), gam(kk,:),gsq(kk,:), m(kk,1), L(kk,:)] = parameters(IEs(kk,:),mus,IW,IEinit(kk,1)) ;
    
    % alternative approach
    m_alt(kk,1) = 6*IEinit(kk,1)/IW;
    prodgsqrtf0c = IEs(kk,:)/6/IEinit(kk,1);
    invgam_alt = polyval(prodgsqrtf0c_gaminv_coef, prodgsqrtf0c);
    gam_alt(kk,:) = 1./invgam_alt;
    
    sqrtf0c_alt = polyval(sqrtf0_gaminv_coef,1./gam_alt(kk,:));
    kappa_alt(kk,:) = sqrtf0c_alt.^2.*6*IEinit(kk,1)*IW;
    
end

sqrtf0c = polyval(sqrtf0_gaminv_coef,1./gam);
IWs_check = sqrt(kappa./m)./sqrtf0c;
sqrtf0c_alt = polyval(sqrtf0_gaminv_coef,1./gam_alt);
IWs_check_alt = sqrt(kappa_alt./m_alt)./sqrtf0c_alt;

figure(3)
plot(ratios,(IWs_check-IW)/IW,'bo')
hold on
plot(ratios,(IWs_check_alt-IW)/IW,'ro')
hold off

figure(2)
plot(ratios,IEinit_min ,'b-' ,ratios,IEinit_max,'r-',ratios,IEinit,'k-','LineWidth',1.5)
legend('bottom limit','upper limit','IE init')


gamcheck = [min(gam(:)) min(gam_alt(:))]<0.52;
if any(gamcheck)
    warning(['smallest gamma <0.52 in [original , alternative] model: [' num2str([min(gam(:)) min(gam_alt(:))]) ']'])
end

figure(1)
    subplot(221)
        plot(ratios,kappa,'o')
        hold on
        plot(ratios,kappa_alt,'k-')
        hold off
        title('kappa')
        grid on
    subplot(222)
        plot(ratios,gam,'o')
        hold on
        plot(ratios,gam_alt,'k-')
        plot(ratios([1,end]),0.52*[1,1],'k--')
        hold off
        grid on
%         ylim([0.5,3])
        xticks(0:.1:1)
        title('gamma')
    subplot(223)
%         plot(ratios,kappa-kappa_alt,'o')
        dk= (kappa-kappa_alt)./kappa_alt;
        plot(ratios,dk,'o')
        hold on
        plot(ratios([1,end]),[0,0],'k-')
        hold off
        grid on
        title('diff(kappa)')
    subplot(224)
%         plot(ratios,gam-gam_alt,'o')
        dg = (gam-gam_alt)./gam_alt;
        plot(ratios,dg,'o')
        hold on
        plot(ratios([1,end]),[0,0],'k-')
        hold off
        grid on
        title('diff(gamma)')
        
disp(['mean difference relative to the alternative method'])
disp(['kappa: ' num2str(mean(mean(dk)))])
disp(['gamma: ' num2str(mean(mean(dg)))])
disp(['max abs difference relative to the alternative method'])
disp(['kappa: ' num2str(max(max(abs(dk))))])
disp(['gamma: ' num2str(max(max(abs(dg))))])

% % polynomial 1/gam -> sqrt(f0(1/gam)); fitted for gamma 0.52-40
% sqrtf0_gaminv_coef = [-0.072966     0.35784    -0.68325     0.63578    -0.48566     0.53703];
% sqrtf0c = polyval(sqrtf0_gaminv_coef,1./gam);
% IWs_check = sqrt(kappa./m)./sqrtf0c;
% 
% gamma_val =	[0.52	0.53	0.54	0.55	0.6	0.65	0.7	0.75	0.8	0.85	0.9	0.95	1	1.05	1.1	1.15	1.2	1.25	1.3	1.35	1.4	1.45	1.5	1.55	1.6	1.65	1.7	1.75	1.8	1.85	1.9	1.95	2	2.05	2.1	2.15	2.2	2.25	2.3	2.35	2.4	2.45	2.5	2.55	2.6	2.65	2.7	2.75	2.8	2.85	2.9	2.95	3	3.05	3.1	3.15	3.2	3.25	3.3	3.35	3.4	3.45	3.5	3.55	3.6	3.65	3.7	3.75	3.8	3.85	3.9	3.95	4	4.05	4.1	4.15	4.2	4.25	4.3	4.35	4.4	4.45	4.5	4.55	4.6	4.65	4.7	4.75	4.8	4.85	4.9	4.95	5	5.05	5.1	5.15	5.2	5.25	5.3	5.35	5.4	5.45	5.5	5.55	5.6	5.65	5.7	5.75	5.8	5.85	5.9	5.95	6	6.05	6.1	6.15	6.2	6.25	6.3	6.35	6.4	6.45	6.5	6.55	6.6	6.65	6.7	6.75	6.8	6.85	6.9	6.95	7	7.05	7.1	7.15	7.2	7.25	7.3	7.35	7.4	7.45	7.5	7.55	7.6	7.65	7.7	7.75	7.8	7.85	7.9	7.95	8	8.2	8.3	8.4	8.5	8.6	8.7	8.8	8.9	9	9.1	9.2	9.3	9.4	9.5	9.6	9.7	9.8	9.9	10	11	12	13	14	15	16	20	24	32	40]';
% g_function	= [0.098546	0.119839	0.137421	0.1526	0.208989	0.248543	0.279304	0.304484	0.325758	0.344131	0.360261	0.374603	0.387487	0.399161	0.409813	0.419594	0.428622	0.436994	0.444791	0.452079	0.458913	0.465341	0.471405	0.477137	0.48257	0.48773	0.49264	0.49732	0.501788	0.506062	0.510155	0.51408	0.517849	0.521473	0.524961	0.528321	0.531563	0.534692	0.537716	0.540641	0.543472	0.546215	0.548873	0.551452	0.553956	0.556388	0.558752	0.561051	0.563289	0.565468	0.56759	0.569659	0.571676	0.573644	0.575565	0.577441	0.579274	0.581064	0.582815	0.584528	0.586203	0.587843	0.589449	0.591022	0.592563	0.594073	0.595553	0.597005	0.598429	0.599826	0.601197	0.602543	0.603865	0.605164	0.606439	0.607693	0.608925	0.610136	0.611327	0.612498	0.613651	0.614785	0.615901	0.616999	0.618081	0.619145	0.620194	0.621227	0.622245	0.623248	0.624236	0.625211	0.626171	0.627118	0.628052	0.628973	0.629882	0.630778	0.631663	0.632535	0.633397	0.634247	0.635087	0.635916	0.636734	0.637543	0.638341	0.63913	0.639909	0.640679	0.64144	0.642192	0.642935	0.643669	0.644395	0.645113	0.645823	0.646525	0.647219	0.647906	0.648585	0.649257	0.649922	0.650579	0.65123	0.651874	0.652511	0.653142	0.653766	0.654384	0.654996	0.655602	0.656201	0.656795	0.657383	0.657966	0.658542	0.659114	0.65968	0.66024	0.660796	0.661346	0.661891	0.662432	0.662967	0.663498	0.664024	0.664545	0.665061	0.665573	0.666081	0.666584	0.667083	0.669037	0.66999	0.670927	0.671849	0.672757	0.67365	0.674529	0.675395	0.676247	0.677087	0.677915	0.67873	0.679533	0.680325	0.681106	0.681876	0.682635	0.683384	0.684122	0.690998	0.697079	0.702515	0.707419	0.711877	0.715955	0.729414	0.73972	0.75484	0.765691]';
% IE_in = interp1(gamma_val,g_function,gam).*sqrt(kappa*m);
% 
% figure(1)
% plot(0.52:.01:40,polyval(sqrtf0_gaminv_coef,1./(0.52:.01:40)))
% hold on
% plot(gam,sqrtf0c,'o')
% hold off
% 
% % alternative way
% prodgsqrtf0c_gaminv_coef = [103.397      -165.393      105.3469     -44.55661       24.7348     -11.25718      1.999642];
% 
% IEinit_min = max(IEs)./6/0.4065;
% IEinit_max = min(IEs)./6/0.0069;
% assert(IEinit>IEinit_min&IEinit<IEinit_max)
% % [kappa, gam,gsq, m, L] = parameters([0.1 0.2 0.3],[1 1 1],1e-8,0.25) 
% 
% m_alt = 6*IEinit/IW;
% prodgsqrtf0c = IEs/6/IEinit;
% invgam_alt = polyval(prodgsqrtf0c_gaminv_coef, prodgsqrtf0c);
% gam_alt = 1./invgam_alt;
% 
% sqrtf0c_alt = polyval(sqrtf0_gaminv_coef,1./gam_alt);
% kappa_alt = sqrtf0c_alt.^2.*6*IEinit*IW;