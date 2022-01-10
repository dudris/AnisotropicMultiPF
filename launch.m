%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Three strategies of parameters assignment
% in multi phase field model of grain growth with anisotorpic grain boundary properties”, Mendeley Data, 
% v1 http://dx.doi.org/10.17632/5wrv3ky9pp.1>, coupled to publication of the same name by 
% Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.

%% EXAMPLE - run a custom simulation defined in make_input.m
clear 
addpath('input\','solver\')

in = make_input;
% compute model parameters and time step
in = input_calc_PFpar_dt(in);

[p,A,F,in] = run_simulation(in);

figure(1)
imagesc(p{2}), colorbar, axis equal
% validate the energy output
%1...total IE, 2 ... mean total IE

[ttt,~] = get_sim_timeline(in);
xxlim = ttt(end)*[0.05,1];
figure(11)
area = A(:,2)*in.dx^2*in.Nx*in.Ny;

subplot(121)
    plot(ttt,F(:,1),'b.')
    hold on
    CW = 1-in.intf.params_incl_dep.soaIE*in.intf.params_incl_dep.Omega/2;
    totIE_anal = 2*in.intf.IE_phases*sqrt(area*pi*CW);
    plot(ttt,totIE_anal,'b--')
    hold off
    xlim(xxlim)
subplot(122)
    plot(ttt,F(:,2),'r.')
    hold on 
    meanIE_anal = in.intf.IE_phases*CW;
    plot(ttt,meanIE_anal*ones(size(ttt)),'r--.')
    hold off
    xlim(xxlim)

    sumpsq = p{1}.^2+p{2}.^2;
    S_ctr = cellfun(@(x) mean(mean(x.^2./sumpsq)),p)';
    
%% EXAMPLES - sample simulations from the Paper
% in 'input\examples\' are ready-made input files with sample simulations
% of the four experiments from the Paper
clear
addpath('input\examples\','solver\','input\')

% uncomment the line corresponding to the simulation to be run
% in = input_shrinking_circles; % IWvK, IEtop = 0.3IEbot, IEbot = 0.3 J/m^2
% in = input_trijunction; % IWvG, ratio 1/0.6=1.667
% in = input_wulff; % IWvK, Omega=5, fourfold
% in = input_compensated_aniso; % IWc, Omega=0.6, fourfold

% compute model parameters and time step
in = input_calc_PFpar_dt(in);
% run the simulation
[p,A,F,in] = run_simulation(in);

%% EXAMPLE - only determination of parameters
% see help comments of the functions below for more details
clear 
IEminmax = [0.1, 0.3]; % minimal and maximal interface energy in the system
IEs = [0.3 0.3 0.2 0.1]; % interface energies of the pair-wise interfaces. 
model = 'IWvG'; % either of 'IWc', 'IWvG' or 'IWvK'

GBmobility = 7.5*1e-16*ones(size(IEs)); % interface mobility in m^4/Js
IWmin = 1e-9; % in meters; no interface will be narrower than IWmin

IEinit = determineIEinit(IEminmax,IEs,model);
[kpp0, gam0, m, L, IWout, gsq] = get_PF_parameters(model,IEs, GBmobility ,IWmin, IEinit);

%% EXAMPLE - visualization of time evolution of area
% usage of get_sim_timeline function

clear
addpath('input\examples\','solver\','input\')

in = input_shrinking_circles; % IWvK, IEtop = 0.3IEbot, IEbot = 0.3 J/m^2
in.plotcond = false; % turn off plotting during simulation
in.ctrcnt = 40; % define number of output points 
in = input_calc_PFpar_dt(in);
[p,A,F,in] = run_simulation(in);

[t_ctr, t_ctr_p]= get_sim_timeline(in);
figure(1)
plot(t_ctr,A,'o-')
xlabel('time (s)')
ylabel('area fraction')
legend('matrix','bottom circle','top circle','Location','east')



