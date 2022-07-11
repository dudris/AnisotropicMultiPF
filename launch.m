%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Three strategies of parameters assignment
% in multi phase field model of grain growth with anisotorpic grain boundary properties”, Mendeley Data, 
% v1 http://dx.doi.org/10.17632/5wrv3ky9pp.1>, coupled to publication of the same name by 
% Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754)

%% EXAMPLE - run a custom simulation defined in make_input.m
clear 
addpath('input\','solver\','processing\')

% ___ creates input file 'in', type 'help make_input' for details
in = make_input;
% ___ adds model parameters and time step to 'in', type 'help input_calc_PFpar_dt' for details
in = input_calc_PFpar_dt(in);

% ___ run the simulation, type 'help run_simulation' for details
[pctr,S,F,in] = run_simulation(in);

% ___ visualize output
figure(1)
imagesc(pctr{2}), colorbar
    set(gca,'YDir','normal','DataAspectRatio',[1,1,1])
    title('phase field 2')
    xlabel('x (grid pts)')
    ylabel('y (grid pts)')
    
%% EXAMPLES - sample simulations from the Paper
% ___ in 'input\examples\' are ready-made input files with sample simulations
% of the four experiments from the Paper
clear
addpath('input\examples\','solver\','input\','processing\')

% ___ uncomment the line corresponding to the simulation to be run
% in = input_shrinking_circles; % IWvK, IEtop = 0.3IEbot, IEbot = 0.3 J/m^2
% in = input_trijunction; % IWvG, ratio 1/0.6=1.667
% in = input_wulff; % IWvK, Omega=5, fourfold
% in = input_compensated_aniso; % IWc, Omega=0.6, fourfold

% ___ compute model parameters and time step
in = input_calc_PFpar_dt(in);
% ___ run the simulation
[pctr,S,F,in] = run_simulation(in);

%% EXAMPLE - only determination of parameters
% ___ see help comments of the functions below for more details
clear 
IEminmax = [0.1, 0.3]; % (J/m^2), minimal and maximal interface energy in the system
IEs = [0.3 0.3 0.2 0.1]; % (J/m^2), interface energies of the pair-wise interfaces. 
model = 'IWvG'; % either of 'IWc', 'IWvG' or 'IWvK', as described in the paper

GBmobility = 7.5*1e-16*ones(size(IEs)); % (m^4/Js), interface mobility
IWmin = 1e-9; % (m), no interface will be narrower than IWmin

% ___ computation of phase field parameters and time step requires
% appropriate value of initializing interface energy IEinit. More details
% in Supplementary material of the Paper
IEinit = determineIEinit(IEminmax,IEs,model);
% ___ get phase field parameters, type help get_PF_parameters' for more 
[kpp0, gam0, m, L, IWout, gsq] = get_PF_parameters(model,IEs, GBmobility ,IWmin, IEinit);

%% EXAMPLE - visualization of time evolution of area in shrinking circles simulation 
% ___ usage of supplementary get_sim_timeline function

clear
addpath('input\examples\','solver\','input\','processing\')

in = input_shrinking_circles; % IWvK, IEtop = 0.3IEbot, IEbot = 0.3 J/m^2
in.plotcond = false; % turn off plotting during simulation
in.ctrcnt = 40; % define number of output points 
in = input_calc_PFpar_dt(in); % compute PF parameters and time step
[pctr,S,F,in] = run_simulation(in); % run the simulation

% ___ extract the simulation time from the input associated with output
% checkpoints
[t_ctr, t_ctr_p]= get_sim_timeline(in);

% ___ visualize time evolution of area fraction of individual phase fields
figure(1)
    plot(t_ctr,S,'o-')
    xlabel('time (s)')
    ylabel('area fraction')
    legend('matrix','bottom circle','top circle','Location','east')

%% EXAMPLE - time evolution of match to Wulff shape and the shrinkage rate
clear
addpath('input\examples\','solver\','input\','processing\')

% __ generate input structure
in = input_wulff; % IWvK
% __ modify input to be returned phase fields in time-equidistant checkpoints
in.outputAtChckpt{1} = true;
% __ set number of the checkpoints
in.outputAtChckpt{2} = 10;

% __ to turn the on-the-fly plotting off
in.plotcond = false;

% __ set the desired normalized strength of anisotropy Omega to value 0-7.5
% __ note that with increased Omega, the initial shape has smaller area for
% constant radius, to be set under in.ICparam(3)
Omega = 5; 
in.intf.params_incl_dep.Omega = Omega; % anisotropy of interface energy
in.intf.params_incl_dep.soaIE = Omega/(4^2-1); % denoted delta in the paper
in.ICparam(4) = Omega; % set the strength of anisotropy of the initial shape

% ___ compute model parameters and time step
in = input_calc_PFpar_dt(in);
% ___ run the simulation
[pctr,S,F,in] = run_simulation(in);

% __ set the polar angles (rad) to have the contour interpolated into (column vector)
phi = linspace(-pi,pi,500);
phi = phi';
% __ process the contours in the checkpoints
% __ type 'help GetContourAndCompareToWulff'
for k = 1:length(pctr)
    [XYcentroid(k,:), RW(k,:) ,r_contour(:,k), norm_sum_squared_diff(k,1),HausdorffD(k,1)] = GetContourAndCompareToWulff(pctr{k}{2}-pctr{k}{1}, phi, in, false);
end

% ___ extract the simulation time from the input associated with output
% checkpoints
[t_ctr, t_ctr_p]= get_sim_timeline(in);

% __ visualize the results in the simulation time
figure(2)
    subplot(121)
    plot(t_ctr_p,HausdorffD,'o-')
    xlabel('time (s)')
    ylabel('Hausdorff Distance (-)')
%     legend('matrix','bottom circle','top circle','Location','east')
    
    subplot(122)
    domarea = in.Nx^2*in.dx^2;
    dAdt = domarea*diff(S(:,2))./diff(t_ctr);
    plot(t_ctr(2:end),dAdt,'-')
    CW = polyval([-0.00032,0.00639,-0.04219,0.00034,1.00000],Omega)/(1-in.intf.params_incl_dep.soaIE);
    dAdt_anal = -2*pi*in.GBmobility*in.IE*CW;
    hold on
    plot(t_ctr([1,end]),ones(1,2)*dAdt_anal,'k--')
    hold off
    ylim([-inf,-5e-16])
    legend('shrinkage rate','analytic')
