clear 
addpath input\
addpath solver\

in = make_input;
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
