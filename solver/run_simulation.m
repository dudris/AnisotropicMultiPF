function [pctr,S,F,in] = run_simulation(in)

%% PrintSimulationSpecification
PrintSimulationSpecification(in)

%% initialization & preallocation 

if in.nOP == 2
    L = in.Lij;
end

global gradX gradY gradXX gradYY gradYX lap size2D simsize
size2D = [in.Ny in.Nx];
[gradX,gradY,gradXX,gradYY,gradYX,lap,simsize] = calcDifferentialOperators(in);

basic_plotparams.size2D = size2D;
basic_plotparams.spacing = [in.dx,in.dy];

if in.Ndtstart_is_1
    p = InitializeOPs(in,false); 
else
    p = in.p;
end
x = in.dx*(1:in.Nx);
y = in.dy*(1:in.Ny);

for i = 1:in.nOP
    p{i} = reshape(p{i},simsize);
end

% assign limiting forbidden angle - zero for weak anisotropy
lim_forbb_ang = Assign_lim_forbb_ang(in);

ndim = 2;
gp = cell(in.nOP,1);
df0dp = cell(in.nOP,1);
df0dGp = cell(in.nOP,1);
DFvG = cell(in.nOP,1);
DFdivG = cell(in.nOP,1);
DFdivK = cell(in.nOP,1);
DFvK = cell(in.nOP,1);
DFgK = cell(in.nOP,1);
DFlapp = cell(in.nOP,1);
% assumes each OP has the SL interface
dhLdp = cell(in.nOP,1);
dhGdp = cell(in.nOP,1);

npairs = in.nOP*(in.nOP-1)/2;
indpairs = combnk(1:in.nOP,2); % combinations of OPs , matrix size (npairs,2) 
nn = cell(npairs,1); % normal at the interface, phi_ij in Neles paper
nnl = cell(npairs,1); % length of nn in a given point
nnlsq = cell(npairs,1);
ori = cell(npairs,1); % normal at the interface, phi_ij in Neles paper

sumpsq = zeros(simsize);
% Fdens = cell(2,1);
Fdens = zeros(simsize);

dkppdGp = cell(in.nOP,1);
dgmmdGp = cell(in.nOP,1);
for i = 1:in.nOP
    dkppdGp_sumpsqpsq{i} = cell(1,ndim);
    dgmmdGp{i} = cell(1,ndim);
    for ddim = 1:ndim
        dkppdGp_sumpsqpsq{i}{ddim} = zeros(simsize);
        dgmmijdGp{i}{ddim} = zeros(simsize);
    end
end

% model specification - switch structure to 
models = categorical({'IWc','IWvG','IWvK'})';
specs = logical([1 1 ; 0 1 ; 1 0 ]);
modspec.kpp = specs(models == in.model,1);
modspec.gam = specs(models == in.model,2);
gkL = struct;
aniso = struct;

% check values
pvalUnder = cell(in.nOP,1);
pvalAbove = cell(in.nOP,1);


% output
if in.outputAtAllCtr{1} || in.outputAtChckpt{1}
    
    if in.outputAtAllCtr{1}
        tctrout_p = int16([ 1 , in.outputAtAllCtr{2}:in.outputAtAllCtr{2}:in.Ndt in.Ndt]);
    else
        tctrout_p = in.Ndt;
    end
    if in.outputAtChckpt{1}
        tchck_p = in.outputAtChckpt{3};
    else
        tchck_p = in.Ndt ;
    end
    
    tctrout_p_all = unique([tctrout_p, tchck_p]);
    numctrout_p = length(tctrout_p_all);
    pctr = cell(numctrout_p,1);
    
    tctrout_all = unique([ in.tsteptsctr , tctrout_p_all]);
    numctrout_all = length(tctrout_all);
else
    tctrout_all = in.tsteptsctr;
    numctrout_all = in.ctrcnt;
end

F = nan(numctrout_all, 2);
ind_singlePF_IEcalc = get_ind_singlePF_calcIE(in.spec_singlePF_IEcalc,in.nOP);
S = nan(numctrout_all , in.nOP);

ctrerr = 0;
ctr_out = 1; % counter for output
ctr = 1;

%% simulation
tstep = in.Ndtstart;
simulation_proceeds = true;
tic
while simulation_proceeds
% for tstep = Ndtstart:Ndt 
    % ploting condition
    PlotAfter = tstep > in.PlotAftertstep;
    plotting = ( any(tstep==tctrout_all ) || PlotAfter) && in.plotcond;

    basic_plotparams.tstep = tstep;
    basic_plotparams.plotting = plotting;
    
    is_locaniso_maincalc_IE = in.is_locally_aniso_IE & tstep >= in.precycle;
    is_incl_dep_maincalc = any([in.is_inclination_dependent_IE;in.is_inclination_dependent_L]) & tstep >= in.precycle;
    is_misor_dep_or_precycle = ( ~in.is_locally_aniso_IE & ~in.is_locally_aniso_L ) | tstep < in.precycle;
    
    
    sumpsq = zeros(simsize); 
    sumgradpsq = zeros(simsize);
    sumpsqpsq = zeros(simsize);
    
    for i = 1:in.nOP
        sumpsq = sumpsq + p{i}.^2;
        lapp{i} = lap(p{i});
        
        if is_incl_dep_maincalc
            gp{i}{1} = gradX(p{i});
            gp{i}{2} = gradY(p{i});
            
            sumpsqpsq = Calc_sumpsqpsq(i,p,sumpsqpsq);
            
            if in.is_inclination_dependent_IE
            sumgradpsq = sumgradpsq + ( gp{i}{1}.^2 + gp{i}{2}.^2 );
            
                if in.all_expr_anal
                    [gradXXp{i}, gradYXp{i}, gradYYp{i},gradsumgradpsq] = Calc_p_derived_fields_indiv_analytically(p{i} , gp{i} , gradsumgradpsq);
                end

                if modspec.kpp
                    % null before next step cumulative calculation
                    DFin.dkppdGp_sumpsqpsq{i}{1} = zeros(simsize);
                    DFin.dkppdGp_sumpsqpsq{i}{2} = zeros(simsize);
                end
            end % if is_incl_dep_maincalc_IE
            
        elseif in.is_misori_dependent
            sumpsqpsq = Calc_sumpsqpsq(i,p,sumpsqpsq);
        end % if is_anisotropic_maincalc
    end % for nOP    
    
    if ~in.all_expr_anal && any(is_locaniso_maincalc_IE) && modspec.kpp
        if in.nOP == 2
            DFin.sumgradpsq_mod = sumgradpsq;
            DFin.grad_sumgrad_over_sumpsqpsq{1} = gradX(sumgradpsq);
            DFin.grad_sumgrad_over_sumpsqpsq{2} = gradY(sumgradpsq);
        else
            DFin.sumgradpsq_mod = sumgradpsq./sumpsqpsq;
            DFin.grad_sumgrad_over_sumpsqpsq{1} = gradX(sumgradpsq./sumpsqpsq);
            DFin.grad_sumgrad_over_sumpsqpsq{2} = gradY(sumgradpsq./sumpsqpsq);
        end % if nOP
    end % if incl. dep. IE with aniso kappa
    %% 
    kappa = zeros(simsize);
    gam_sumpsqpsq = zeros(simsize);
    out_interf_all = false(simsize);
    inn_interf_all = false(simsize);
    
    if in.is_misori_dependent || in.is_inclination_dependent_L
        L = zeros(simsize);
    end
    
    % kappa, gamma, L
    for k = 1:npairs
        i = indpairs(k,1);
        j = indpairs(k,2);
        
        if is_misor_dep_or_precycle(k)
            
            gkL = CalcKapGam_iso(in,k,gkL);

        else  % larger timestep nr AND in.is_locally_isotropic(k)=false
            
            [nnlsq{k} , nnl{k}] = CalcMagnitudeInterfaceNormal_ij(gp{i},gp{j});
            %
            [out_interf_ij , inn_interf_ij] = DetermineInnerOuterInterfaceRegion_ij(is_misor_dep_or_precycle(k), nnl{k}, in.intf.limval,basic_plotparams);
            
            for ddim = 1:ndim
                nn{k}{ddim}= ( gp{i}{ddim}-gp{j}{ddim})./nnl{k};
            end

            ori{k} = nan(simsize);
            ori{k}(out_interf_ij) = atan2(nn{k}{2}(out_interf_ij),nn{k}{1}(out_interf_ij));
            
            [gkL, aniso] = CalcKapGam_aniso(in, modspec, lim_forbb_ang, k, ori{k}, gkL, aniso);
        
            out_interf_all = out_interf_all | out_interf_ij;
            inn_interf_all = inn_interf_all | inn_interf_ij;
        end % if isotropic vs aniso
        
        [kappa, gam_sumpsqpsq, L] = CalcKappaGammaL_Add_kth_Term(kappa, gam_sumpsqpsq, L, gkL, in, p, i, j, k);
        
    end % for npairs
    
    if in.nOP > 2 
        if modspec.kpp
            kappa = kappa./sumpsqpsq;
        end
        L = L./sumpsqpsq;
    end
        
    if any(is_locaniso_maincalc_IE) && modspec.kpp
        DFin.gradkpp{1} = gradX(kappa);
        DFin.gradkpp{2} = gradY(kappa);
    end
    
    %% evolution equation
    
    for i = 1:in.nOP
        
        if any(is_locaniso_maincalc_IE) 
            
            if modspec.kpp
                DFin.divdkppdGp_mod{i} = zeros(simsize);
            end
            
            % for explanation see see Gderivative_fields_summing_algorithm_example5nOP.png
            % allows to do npairs calculation in total, not nOP(nOP-1)
            ind_rows_with_i_and_j_larger = FindRowsWith_i_And_j_Larger(indpairs,i); 
            
            condBeforeDF = out_interf_all;
            
            for l = 1:length(ind_rows_with_i_and_j_larger)
                k = ind_rows_with_i_and_j_larger(l);        

                if in.is_locally_aniso_IE(k) % l-th term added to dkppdGp{i} is nonzero only when that k-th interface is not locally isotropic in IE
                    Gdermultips_ij = CalcGderPrefactors(in,k,aniso);
                    [dkppijdGp_ij,dgmmijdGp{k}] = CalcGderivativeOfKappaGamma_ij(condBeforeDF,nn{k},nnl{k},Gdermultips_ij,modspec);
                    
                    if modspec.kpp
                        j = GetIndexOfTheOther(k,indpairs,i);
                        DFin.dkppdGp_sumpsqpsq = CalcGderivativeOfKappa_i_And_j_AddTerm(DFin.dkppdGp_sumpsqpsq,dkppijdGp_ij, p, in.nOP, i, j);
                    end
                end % if NOT is_locally_isotropic
                
            end %for l = 1:length(ind_rows_with_i_and_j_larger)
            
            if modspec.kpp
                for ddim = 1:ndim
                    DFin.dkppdGp_sumpsqpsq{i}{ddim}(isnan(DFin.dkppdGp_sumpsqpsq{i}{ddim})) = 0;
                end % for ddim
            end
                
            if modspec.kpp 
                DFin.divdkppdGp_mod{i}  = gradX(DFin.dkppdGp_sumpsqpsq{i}{1}) + gradY(DFin.dkppdGp_sumpsqpsq{i}{2});
                DFin.divdkppdGp_mod{i}(isnan(DFin.divdkppdGp_mod{i})) = 0;
            end
            
        end %if is_incl_dep_maincalc
        
        if modspec.gam
            sumgam_ij_psq_i = zeros(simsize);
            if any(is_locaniso_maincalc_IE)
                DFin.sumdgamdGp_ij_div_i = zeros(simsize);
                DFin.sumdgamdGp_ij_vec_i{1} = zeros(simsize);
                DFin.sumdgamdGp_ij_vec_i{2} = zeros(simsize);
                DFin.sumdgamdGp_ij_vecdotpr_i = zeros(simsize);
            end % if is_inclination_dependent
        
            rows_w_i = find(any(indpairs==i,2));
            SIGN = [1 , -1];
            for l = 1:length(rows_w_i) % pairwise interfaces l ~= i
                k = rows_w_i(l);
                j = GetIndexOfTheOther(k,indpairs,i);
                sign = SIGN(indpairs(k,:)==i); % when i is at the 1st place +1, when at the 2nd -1
                pjsq = p{j}.^2;
                if is_locaniso_maincalc_IE(k) % l-th term added to dgmmijdGp{i} is nonzero only when that k-th interface is not locally isotropic
                    DFin.sumdgamdGp_ij_vec_i{1} = DFin.sumdgamdGp_ij_vec_i{1} + sign*dgmmijdGp{k}{1}.*pjsq;
                    DFin.sumdgamdGp_ij_vec_i{2} = DFin.sumdgamdGp_ij_vec_i{2} + sign*dgmmijdGp{k}{2}.*pjsq;
                    DFin.sumdgamdGp_ij_vecdotpr_i = DFin.sumdgamdGp_ij_vecdotpr_i + sign*p{j}.*( gp{j}{1}.*dgmmijdGp{k}{1} + gp{j}{2}.*dgmmijdGp{k}{2} );
                    DFin.sumdgamdGp_ij_div_i = DFin.sumdgamdGp_ij_div_i + sign*pjsq.*( gradX(dgmmijdGp{k}{1}) + gradY(dgmmijdGp{k}{2}) );
                end

                sumgam_ij_psq_i = sumgam_ij_psq_i + gkL.g_ij{k}.*pjsq;

            end % for
        else
            
            sumgam_ij_psq_i = 1.5*(sumpsq - p{i}.^2);
        end % if modspec.gam
        
    %% new timestep value calculation
        [df0dp{i} , DFlapp{i}] = CalcDrivingForceTerms_Isotropic(p{i}, in.m, sumgam_ij_psq_i, kappa,lapp{i});

        df0dGp{i} = DFlapp{i};
        if any(is_locaniso_maincalc_IE)
            condIntoDF = inn_interf_all;
            normgp_i = sqrt(gp{i}{1}.^2+gp{i}{2}.^2); % to strictly limit DF to interfaces of p{i} only
            condIntoDF = inn_interf_all & (normgp_i >= in.intf.limval.aniso.inner/(2*in.dx));
            
            [DFdivG{i} , DFvG{i} , DFdivK{i} ,  DFvK{i} , DFgK{i}] = CalcDrivingForceTerms_Anisotropic(condIntoDF,in.m, p{i}, gp{i}, DFin,modspec,i);
            df0dGp{i} = df0dGp{i} + DFvG{i} + DFdivG{i} + DFdivK{i}+ DFvK{i} + DFgK{i}; % DFlapp{i} is already in df0dGp{i} 
        end % if is_isotropic
         
        p_new{i} = p{i} - in.dt*L.*(df0dp{i} - df0dGp{i});
        if tstep<=in.precycle % area conserved using Lagrange-multipliers approach
            dinterpf_i = dinterpf(p,i,sumpsq);
            redistr_antiforce(i) = - sum(df0dp{i} - df0dGp{i})/sum(dinterpf_i);
            p_new{i} = p_new{i} - in.dt.*L.*redistr_antiforce(i).*dinterpf_i;
        end % if tstep<=precycle 
        
        % add ith term to energy calculation at checkpoint
        if any(tstep==tctrout_all)
            if i == 1
                Fdens = zeros(simsize);
            end
            
            if any(is_locaniso_maincalc_IE)
                % sumgradpsq calculated in each tstep for driving force
                Fdens = Fdens + in.m*(0.25*p_new{i}.^4 - 0.5*p_new{i}.^2);
            else
                % sumgradpsq not calculated yet
                % in case of specBC and need for longer precycle modify
                sumgradpsq = sumgradpsq + ( gradX(p_new{i}).^2 + gradY(p_new{i}).^2 );
                Fdens = Fdens + in.m*(0.25*p_new{i}.^4 - 0.5*p_new{i}.^2);
            end
        end
        
    end % for nOP
    
    % assign new values to old ones and check the PF values are within (0,1)
    for i = 1:in.nOP
        p{i} = p_new{i};
        % feedback on the evolution
        checkNaNp(i) = all(all(~isnan(p{i}))); % true when there is no NaN
        % tolerance on values of PFs to go out of the interval (0,1)
        pValTol = 1e-2;
        [pvalUnder{i},pvalAbove{i},is_pvalUnder(i),is_pvalAbove(i)] = CheckPValuesInInterval(p{i},pValTol);
    end
    
    %% checkpoint calculations - energy, area, arc length ; print progress
    if any(tstep==tctrout_all)
        
%         if is_incl_dep_maincalc || in.is_misori_dependent
%             if is_incl_dep_maincalc 
%                 normgp_i = sqrt(gp{ind_singlePF_IEcalc}{1}.^2+gp{ind_singlePF_IEcalc}{2}.^2);
%             else
%                 normgp_i = sqrt(gradX(p{ind_singlePF_IEcalc}).^2+gradY(p{ind_singlePF_IEcalc}).^2);
%             end
%             cond_singlePF_IEcalc = normgp_i >= in.intf.limval.aniso.inner/(2*in.dx);
% %             [F(ctr,:), S(ctr,:)] = CalcEnergyAreas(Fdens,p,gam_sumpsqpsq,sumgradpsq,sumpsqpsq,sumpsq,in,tstep,cond_singlePF_IEcalc,kappa);
%     
%         else % precycle or isotropic
%             sumpsqpsq = zeros(simsize);
%             for i = 1:(in.nOP-1) % sumpsqpsq was not calculated in isotropic cycle
%                 for j = (i+1):in.nOP
%                     sumpsqpsq = sumpsqpsq + p{i}.^2.*p{j}.^2;
%                 end
%             end
%             normgp_i = sqrt(gradX(p{ind_singlePF_IEcalc}).^2+gradY(p{ind_singlePF_IEcalc}).^2);
%             cond_singlePF_IEcalc = normgp_i >= in.intf.limval.aniso.inner/(2*in.dx);
%         end
%         [F(ctr,:), S(ctr,:)] = CalcEnergyAreas(Fdens,p,gam_sumpsqpsq,sumgradpsq,sumpsqpsq,sumpsq,in,tstep,cond_singlePF_IEcalc,kappa);
        
        [F(ctr,:), S(ctr,:)] = CalcEnergyAreas(Fdens,p,gam_sumpsqpsq,sumgradpsq,sumpsqpsq,sumpsq,in,tstep,kappa);
        Fdens = zeros(simsize);
        
        if ctr >1
            condSinvariant = S(ctr,2)-S(ctr-1,2)~=0;
            if ~in.is_conserved && ~in.is_conc_conserved
                if ~condSinvariant
                    warning(['Simulation terminated at simtime ' num2str(in.dt*tstep) '  Area of p{2} invariant last ' num2str() ' tsteps'])
                    in.Ndt = tstep;
                    in.simtime = tstep*in.dt;
                    break
                end
%                 assert(condSinvariant, ['Simulation terminated at simtime ' num2str(in.dt*tstep) '  Area of p{2} invariant last ' num2str(in.ctrplot) ' tsteps'])
            end
        end
        
        if in.is_cond_termd.bool && tstep>in.precycle% to be conditionally terminted
            obs_PF = in.is_cond_termd.PFnum;
            if strcmp(in.is_cond_termd.code,'area_ratio')
                area_to_stop = S(1,obs_PF)*in.is_cond_termd.ratio;
                cond_sim_ctd = S(ctr,obs_PF) > area_to_stop;
            elseif strcmp(in.is_cond_termd.code,'area_abs')
                area_to_stop = in.is_cond_termd.area;
                cond_sim_ctd = S(ctr,obs_PF)*(in.Ny*in.Nx*in.dx*in.dx) > area_to_stop; 
            elseif strcmp(in.is_cond_termd.code,'1g_totIE') || strcmp(in.is_cond_termd.code,'1g_meanIE')
%                 energy_chng_rate = sum(F(ctr-[2,1,0],1).*[1/2, -2 , 3/2])/in.dt; % looking for minimum - deprecated
                ptsno_formean = ceil(in.is_cond_termd.meantime/in.ctrplot/in.dt); 
                ptsno_formean = min([ptsno_formean,ctr]); % in the first ctrplots compute mean from available only
                ind_formean = (ctr-(ptsno_formean-1)):ctr;
                if strcmp(in.is_cond_termd.code,'1g_totIE')
                    energcol = 1;
                elseif strcmp(in.is_cond_termd.code,'1g_meanIE')
                    energcol = 2;
                end
                mean_energ_somesteps = mean(F(ind_formean,energcol));
                % when value close enough to mean from the last 3 ctrsteps, terminate
                if (tstep*in.dt) > in.is_cond_termd.meantime
                    cond_sim_ctd = abs(mean_energ_somesteps-F(ctr,energcol))/F(ctr,energcol) >= in.is_cond_termd.rellim;
                else
                    cond_sim_ctd = true;
                end
            else % ctr <= 2 
                cond_sim_ctd = true;
            end
            
            if ~cond_sim_ctd % simulation to be terminated
                in.Ndt = tstep; % adjust the total number of time steps
                in.simtime = tstep*in.dt;
            end
        else % no conditional termination
            cond_sim_ctd = true;
        end % if terminate conditionally
        
        ctr=ctr+1;
        t_temp = toc;
        disp(['time step ' num2str(tstep) '/' num2str(in.Ndt) ' => ' num2str(tstep/in.Ndt*100) '%   (' num2str(t_temp) ')'])
        if tstep~=in.Ndt
            tic
        end
%         if in.is_conc_conserved && tstep > precycle
%             disp(['    sum(c)dx^2 = ' num2str(sum(c)*dx^2)])
%         end
    end % if checkpoint calculation
    
    %% plotting
    if plotting % || tstep > in.precycle
%         % figure(1)
        if in.nOP >2 && modspec.gam
            gam = gam_sumpsqpsq./sumpsqpsq;
        else
            gam = gam_sumpsqpsq;
        end
        
        if modspec.kpp
            partoplot = kappa.*ones(simsize);
            partoplot_nm = '\kappa';
        else
            partoplot = gam.*ones(simsize);
            partoplot_nm = '\gamma';
        end
        
        if in.is_conc_conserved && in.is_inclination_dependent_IE && tstep > in.precycle
            ind_ori_exists = find(~in.is_locally_isotropic);
            ind_ori_exists = ind_ori_exists(1);
            plot2D_from_lin_2x3(size2D,{sumpsq,p_new{1},p_new{2},ori{ind_ori_exists},partoplot,c},{['tstep = ' num2str(tstep) ',\Sigma \xi_i^2\xi_j^2' ],'\xi_1','\xi_2','ori_1',partoplot_nm,'c'})
            
        elseif ~in.is_conc_conserved && in.is_inclination_dependent_IE && tstep > in.precycle
            ind_ori_exists = find(in.is_locally_aniso_IE);
            ind_ori_exists = ind_ori_exists(1);
%             plot_basic_info(sumpsq,p_new,ori{ind_ori_exists},kappa,nnl{ind_ori_exists},size2D,tstep,zeros(simsize))
            plot2D_from_lin_2xn(1, {sumpsq,partoplot, p{1},p{2}},{'sumpsqpsq',partoplot_nm,'xi_1','\xi_2'} , size2D, 2,false)
            
        elseif (~in.is_conc_conserved && ~in.is_inclination_dependent_IE ) || tstep <=in.precycle
            if in.nOP == 3
                plot2D_from_lin_2xn(1, {sumpsq,partoplot,ones(simsize), p{1},p{2},p{3}},{'sumpsqpsq',partoplot_nm,'blank space','xi_1','\xi_2','\xi_3'} , size2D, 3,false)
            elseif in.nOP==2
                
                plot2D_from_lin_2xn(1, {sumpsq,partoplot, p{1},p{2}},{'sumpsqpsq',partoplot_nm,'xi_1','\xi_2'} , size2D, 2,false)
            end
            
        elseif in.is_conc_conserved && ~in.is_inclination_dependent_IE
            if in.nOP == 3
                plot2D_from_lin_2x3(size2D,{sumpsq,partoplot,c,p{1},p{2},p{3}},{['tstep = ' num2str(tstep) ', \Sigma \xi_i^2\xi_j^2' ], partoplot_nm, 'c', '\xi_1','\xi_2','\xi_3'})
            elseif in.nOP==2
                plot2D_from_lin_2xn(1, {sumpsq,partoplot, p{1},p{2}},{'sumpsqpsq',partoplot_nm,'xi_1','\xi_2'} , size2D, 2,false)
            end
        end
        
%         plot2D_from_lin_2xn(1, {df0dGp{1},df0dGp{2}, df0dp{1},df0dp{2}},{'df0dGp 1','\df0dGp 2','df0dp 1','\df0dp 2'} , size2D, 2,false)
%         folder= 'figs/stronganiso_strange_4fold_OMG4-5';
%         savepic(folder,tstep)
        
        %%
        if any(is_locaniso_maincalc_IE) && in.plotDF % plotting driving forces
            i = 2;
            plot2D_from_lin_2xn(2, {in.dt*L.*DFlapp{i}, in.dt*L.*DFdivG{i}.*ones(simsize) , in.dt*L.*DFvG{i}.*ones(simsize) ,...
                in.dt*L.*DFdivK{i}.*ones(simsize) ,  in.dt*L.*DFvK{i}.*ones(simsize) , in.dt*L.*DFgK{i}.*ones(simsize)}, ...
                {['in.dt*L.*DFlapp ' num2str(i)], ['in.dt*L.*DFdivG ' num2str(i)] , ['in.dt*L.*DFvG ' num2str(i)] , ['in.dt*L.*DFdivK ' num2str(i)] , ['in.dt*L.*DFvK ' num2str(i)] , ['in.dt*L.*DFgK ' num2str(i)]} , size2D, 3,false)
                
%             i = 4;
%             plot2D_from_lin_2xn(3, {in.dt*L.*DFlapp{i}, in.dt*L.*DFdivG{i} , in.dt*L.*DFvG{i} , in.dt*L.*DFdivK{i} ,  in.dt*L.*DFvK{i} , in.dt*L.*DFgK{i}}, {['in.dt*L.*DFlapp ' num2str(i)], ['in.dt*L.*DFdivG ' num2str(i)] , ['in.dt*L.*DFvG ' num2str(i)] , ['in.dt*L.*DFdivK ' num2str(i)] , ['in.dt*L.*DFvK ' num2str(i)] , ['in.dt*L.*DFgK ' num2str(i)]} , size2D, 3,false)

%             filterCond2D = reshape(condIntoDF,size2D);
%             iplot = 1;
%             [DF_i,DFname] = DFtoMatrix_i(in.dt,L(),size2D,is_incl_dep_maincalc,df0dp{iplot},DFlapp{iplot},DFdivG{iplot},DFdivK{iplot},DFvG{iplot},DFvK{iplot},DFgK{iplot});
%         
%       figure(9)
%         plot_RHS_overall(in,iplot,nn{1},p{2},tstep,DF_i,DFname,in.dt,L,df0dp,df0dGp)
% %       figure(55) and figure(56)
%             plot_anisotropic_DF(DF_i,DFname,filterCond2D,iplot)
        
%             figure(88)
% % %             scalar terms
%             DFdivG_names = {'\Sigma_{i,j}\xi_i^2\xi_j^2',['\nabla\cdot \partial\gamma/\partial(\nabla\xi) ' num2str(iplot)], ['DFdivG' num2str(iplot)]};
%                 plot2D_from_lin_scalar(size2D,iplot,{sumpsqpsq,divdgmmdGp,DFdivG},DFdivG_names)
%             DFdivK_names = {'\Sigma_i |\nabla\xi_i|^2',['\nabla\cdot \partial\kappa/\partial(\nabla\xi) ' num2str(iplot)], ['DFdivK' num2str(iplot)]};
%                 plot2D_from_lin_scalar(size2D,iplot,{sumgradpsq,divdkppdGp,DFdivK},DFdivK_names)
%             DFlapp_names = {['\nabla^2\xi ' num2str(iplot)],'\kappa ', ['DFlapp' num2str(iplot)]};
%                 plot2D_from_lin_scalar(size2D,iplot,{lapp{iplot},kappa,DFlapp{iplot}},DFlapp_names)
% % %             dot product
%             plot2D_from_lin_dotpr(size2D,gradsumpsqpsq,'gradsumpsqpsq',dgmmdGp{iplot},['\partial\gamma/\partial(\nabla\xi) ' num2str(iplot)],DFvG{iplot},['DFvG ' num2str(iplot)])
%             plot2D_from_lin_dotpr(size2D,gradsumgradpsq,'\nabla(\Sigma_i|\nabla\xi_i|^2)',dkppdGp{iplot},['\partial\kappa/\partial(\nabla\xi) ' num2str(iplot)],DFvK{iplot},['DFvK ' num2str(iplot)])
%             plot2D_from_lin_dotpr(size2D,gp{iplot},['\nabla\xi ' num2str(iplot)],gradkpp,'\nabla\kappa',DFgK{iplot},['DFgK ' num2str(iplot)])

        end % if ~is_isotropic
        
        if in.PauseAfterPlotting 
            pause
        end
    end% if plotting

    %% check values of p's
%         check for NaNs
        assert(all(checkNaNp),['NaN found in p(s): [' num2str(find(checkNaNp)) '] , tstep = ' num2str(tstep)]) % error if there was a NaN found in any p
%         check of values between (0,1)

    if any(is_pvalUnder | is_pvalAbove) % if any OP value gets from (0,1)
        % plotting the DF
        if modspec.kpp
            partoplot = kappa.*ones(simsize);
            partoplot_nm = '\kappa';
            plot2D_from_lin_2xn(1, {sumpsq,partoplot, p{1},p{2}},{'sumpsqpsq','\kappa','xi_1','\xi_2'} , size2D, 2,false)
        else
            if in.nOP >2 && modspec.gam
                gam = gam_sumpsqpsq./sumpsqpsq;
            else
                gam = gam_sumpsqpsq;
            end
            partoplot = gam.*ones(simsize);
            partoplot_nm = '\gamma';
            plot2D_from_lin_2xn(1, {sumpsq,partoplot, p{1},p{2}},{'sumpsqpsq','\kappa','xi_1','\xi_2'} , size2D, 2,false)
        end
        
%         figure(15)
%         imagesc(reshape(L,size2D)), colorbar, axis equal
        
        for kk = 1:in.nOP
            RHS{kk} = -in.dt*L.*(df0dp{kk} - df0dGp{kk});
            pout{kk} = p{kk};
            pout{kk}(p{kk}>-pValTol & p{kk}<1+pValTol) = NaN;
            titles{kk} = ['PF' num2str(kk) ' outside (0,1)'] ;
            titles{kk+in.nOP} = ['-dt.L.RHS ' num2str(kk)] ;
        end
        
        nan_alpha = true;
        plot2D_from_lin_2xn(99, [pout , RHS],titles , size2D, in.nOP,nan_alpha)
        
        for kk = 1: in.nOP
            RHS_{kk} = -in.dt*L.*df0dp{kk};
            RHS_{in.nOP+kk} = -in.dt*L.*(-df0dGp{kk});
            titles{kk} = ['-dt.L. df0dp ' num2str(kk)] ;
            titles{kk+in.nOP} = ['-dt.L. (-df0dGp)' num2str(kk)] ;
        end
        plot2D_from_lin_2xn(98, RHS_,titles , size2D, in.nOP,false)
                
         % Probing the DF behind outliers         
%         for i = 1:in.nOP
%             disp(['OP ' num2str(i)])
%             if is_pvalUnder(i)
%                 DFunder = GetDFInOutliers(pvalUnder{i},'UNDER',is_incl_dep_maincalc,in.dt,L,df0dp{i},DFlapp{i},DFdivG{i},DFvG{i},DFdivK{i},DFvK{i},DFgK{i});
% %                 [xcoorU,ycoorU] = GetXYCoordOfOutliers(pvalUnder{i},Ny);
%                 disp(['DFunder sum ' num2str(sum(DFunder))])
%             end % if errUNDER(i)
%             
%             if is_pvalAbove(i)
%                 DFabove = GetDFInOutliers(pvalAbove{i},'ABOVE',is_incl_dep_maincalc,in.dt,L,df0dp{i},DFlapp{i},DFdivG{i},DFvG{i},DFdivK{i},DFvK{i},DFgK{i});
% %                 [xcoorA,ycoorA] = GetXYCoordOfOutliers(pvalAbove{i},Ny);
%                 disp(['DFabove sum ' num2str(sum(DFabove))])
%             end % if errABOVE
%                 clear DFunder DFabove
%         end%for nOP
               
       % throw standard error
       ctrerr = ctrerr + 1;
       if ctrerr <= 0
           warning(['values p < 0 found in phase fields: [' num2str(find(is_pvalUnder)) '] \n%s'], ['values p > 1 found in phase fields: [' num2str(find(is_pvalAbove)) ']'])
       else
           disp(['terminated after time step ' num2str(tstep)])
%                error(['terminated after time step ' num2str(tstep) ' \n%s'], ['values p < 0 found in phase fields: [' num2str(find(errUNDER)) '] \n%s'], ['values p > 1 found in phase fields: [' num2str(find(errABOVE)) ']'])
%                 error(['terminated after time step ' num2str(tstep) ' , errUNDER = [' num2str(find(errUNDER)) '], errABOVE = [' num2str(find(errABOVE)) ']'])
            error(['tsteps = ' num2str(tstep) ', sim.time = ' num2str(tstep*in.dt) ' , is_pvalUnder = [' num2str(find(is_pvalUnder)) '], is_pvalAbove = [' num2str(find(is_pvalAbove)) '], Srel = ' num2str(mean(p{2})/S(1,2)) ','])
       end%if ctrerr
    end % if any OP value gets from (0,1)
    %%
    
    simulation_proceeds = (tstep<in.Ndt) & cond_sim_ctd ; 
    
    if strcmp(in.solvermethod,'2Dlin')
            if in.outputAtAllCtr{1} || in.outputAtChckpt{1}
                if any(tstep==tctrout_p_all) || ~simulation_proceeds 
                    for i = 1:in.nOP
                        assert(ctr_out<=length(pctr))
%                         outctr = floor(tstep/in.outputAtAllCtr{2}) + 1;
                        pctr{ctr_out}{i} = reshape(p{i},size2D);
                    end % for
                    ctr_out = ctr_out+1;
                end% if to print out
            end
    end
    
    DispMsgPrecycleTerminated(in.precycle,tstep,any([in.is_inclination_dependent_L,in.is_inclination_dependent_IE]))
    tstep = tstep + 1;
end % while simulation_proceeds
% toc

if ~in.outputAtAllCtr{1} && ~in.outputAtChckpt{1} % not to save at each in.ctrplot
    if strcmp(in.solvermethod,'2Dlin')
        for i = 1:in.nOP
            pctr{i} = reshape(p{i},size2D);
        end
    else % for 1D in.solvermethod
        pctr = p;
        c = reshape(c,size2D);
    end % if in.solvermethod
else % if outputAtAllCtr of outputAtChckpt
%     ind_to_keep = zeros(0,1);
%     for k = 1:length(pctr)
        % when simulation ends conditionally
%         if ~isempty(pctr{k})
%             ind_to_keep = [ind_to_keep, k];
%         end
%     end
%     pctr = pctr(ind_to_keep);
    cond_keep_ctr = cellfun(@(x) ~isempty(x), pctr);
    pctr = pctr(cond_keep_ctr);
end

if  ~in.is_conc_conserved
    c = nan;
end

nan_timesteps =  any(isnan(S),2);
if any(nan_timesteps)
    S(nan_timesteps,:) = [];    
    F(nan_timesteps,:) = [];
end

if cond_sim_ctd 
    disp(['Simulation terminated @tstep ' num2str(in.Ndt) ': max number of tsteps reached.'])
else
    disp(['Simulation terminated @tstep ' num2str(in.Ndt) ': condition for termination accomplished - ' in.is_cond_termd.code '.'])
end

end % func
% 
% 
%% savepic
function savepic(folder,tstep)
    if exist(folder)~=7
        mkdir(folder)
    end
    filename = sprintf([folder '/%04d.png'],tstep);
    print('-dpng',filename)
end

%% PrintSimulationSpecification
function PrintSimulationSpecification(in)

disp('*****************************')
disp('*****************************')
disp('Simulation started')
disp(['Model: ' in.model])

    if in.is_conserved
        disp('Volume conserving model: Lagrange multipliers')
    elseif in.is_conc_conserved
        disp(['Volume conserving model: concentration field.'])
        disp(['Number of conc. fields: ' num2str(in.conserve_2PFs+1) '. Field 1 conserves PF(s): ' num2str(in.PF_to_conserve{1}) ', field 2 conserves PF(s): ' num2str(in.PF_to_conserve{2})])
    elseif (~in.is_conserved && ~in.is_conc_conserved)
        disp('Model not conserving volume.')
    end % if volume conservation
    
    anisobool = [in.is_inclination_dependent_IE , in.is_inclination_dependent_L, in.is_misori_dependent];
    if any(anisobool)
        disp([ 'ANISOTROPY ON: [ IE incl.dep. , L incl.dep. , misorientation ] = ' num2str(anisobool)])
    else
        disp('ANISOTROPY OFF: isotropic model.')
    end
    
    disp(['Courant nr = ' num2str(in.Courant_nr) ', dt = ' num2str(in.dt)])
    disp(['Simtime = ' num2str(in.simtime) ', Ndt = ' num2str(in.Ndt)])
    if in.is_cond_termd.bool
        disp(['Conditional termination code: ' in.is_cond_termd.code])
    else
        disp('Termination when Ndt reached.')
    end
    disp(['BC: ' in.BCs ', solver method: ' in.solvermethod ' , laplacian: ' in.laplacianmethod ])
end

%% get_ind_singlePF_calcIE
function ind_singlePF_IEcalc = get_ind_singlePF_calcIE(spec_singlePF_IEcalc,nOP)
    if ischar(spec_singlePF_IEcalc) % string input
        if strcmp(spec_singlePF_IEcalc,'default')
            ind_singlePF_IEcalc = nOP;
        end
    elseif isnumeric(spec_singlePF_IEcalc)
        ind_singlePF_IEcalc = spec_singlePF_IEcalc;
    end
end % func get_ind_singlePF_calcIE

%% calcDifferentialOperators_scpecBC
% only bottom and top to be Nspecial => only gradYspBC needed
function [gradXspBC, gradYspBC] = calcDifferentialOperators_specBC(in)

m = in.Nx*in.Ny;

if strcmp(in.BCs,'Nspecial')
    assert(strcmp(in.solvermethod,'2Dlin'),'Incompatible input. BCs=Nspecial but ''solvermethod''~=2Dlin.')
    
%     S = sparse(i,j,s,m,n,nzmax) uses vectors i, j, and s to generate an
%     m-by-n sparse matrix such that S(i(k),j(k)) = s(k), with space
%     allocated for nzmax nonzeros.
    m = in.Ny;
    Yder_onecolumn = sparse([1,1,1,m,m,m],[1,2,3,m-2,m-1,m],[-3/2,2,-1/2,1/2,-2,3/2],m,m,6);
    gradYspBC = kron(speye(in.Nx),Yder_onecolumn)/in.dx;
%     imagesc(gradYspBC), axis equal, colorbar
    
    gradXspBC = kron(Yder_onecolumn,speye(in.Nx))/in.dx;
%     imagesc(gradXspBC), axis equal,colorbar
    
else % when BCs is not 'Nspecial' I don't need anything else
    m = in.Nx*in.Ny;
    gradXspBC = sparse([],[],[],m,m,0);
    gradYspBC = sparse([],[],[],m,m,0);
end % if BCs Nspecial
    
end % func

%% AssignTh_SpecBC
% angles in BCs_specs specify angle norm1 between NORMAL to the interface (pointing into PF_1) and NORMAL to x axis 
%   norm1 + norm2 = pi
% th_specBC must be polar angle taken between NORMAL to the interface and x axis
%   pol1 - pol2 = pi
%   pol1 = pi/2 + norm1
%   pol2 = norm1 - pi/2
% norm1 is in.BCs_specs (in degrees)

function th_specBC = AssignTh_SpecBC(in)

    % if Nspecial
    if strcmp(in.BCs,'Nspecial')
        % bottom
        if ~isempty(in.BCs_specs.bottom)
            th_specBC.bot = zeros(in.Nx,1);
            for kk = 1:size(in.BCs_specs.bottom,1) % for number of rows
                indth = in.BCs_specs.bottom(kk,1):in.BCs_specs.bottom(kk,2);
                th_specBC.bot(indth) = in.BCs_specs.bottom(kk,3)/180*pi;
%                 th_specBC.bot(indth) = pi/2 + in.BCs_specs.bottom(kk,3)/180*pi;
            end
            assert(in.nOP ==2,'error AssignTh_SpecBC: currently works for 2 PFs only.')
            th_specBC.bot(:,2) = pi - th_specBC.bot(:,1); % the second PF to be antiparallel
%             th_specBC.bot(:,2) = th_specBC.bot(:,1) - pi; % the second PF to be antiparallel
        else
            th_specBC.bot = nan(1,2);
        end % if bottom angle assigned
        
        % top
        if ~isempty(in.BCs_specs.top)
            th_specBC.top = zeros(in.Nx,1);
            for kk = 1:size(in.BCs_specs.top,1) % for number of rows
                indth = in.BCs_specs.top(kk,1):in.BCs_specs.top(kk,2);
                th_specBC.top(indth) = in.BCs_specs.top(kk,3)/180*pi;
%                 th_specBC.top(indth) = pi/2 + in.BCs_specs.top(kk,3)/180*pi;
            end
            assert(in.nOP ==2,'error AssignTh_SpecBC: currently works for 2 PFs only.')
            th_specBC.top(:,2) = pi - th_specBC.top(:,1); % the second PF to be antiparallel
%             th_specBC.top(:,2) = th_specBC.top(:,1) - pi; % the second PF to be antiparallel
        else
            th_specBC.top = nan(1,2);
        end % if top angle assigned
        
        % left
        if ~isempty(in.BCs_specs.left)
            th_specBC.left = zeros(in.Ny,1);
            for kk = 1:size(in.BCs_specs.left,1) % for number of rows
                indth = in.BCs_specs.left(kk,1):in.BCs_specs.left(kk,2);
                th_specBC.left(indth) = in.BCs_specs.left(kk,3)/180*pi;
%                 th_specBC.top(indth) = pi/2 + in.BCs_specs.top(kk,3)/180*pi;
            end
            assert(in.nOP ==2,'error AssignTh_SpecBC: currently works for 2 PFs only.')
            th_specBC.left(:,2) = pi - th_specBC.left(:,1); % the second PF to be antiparallel
%             th_specBC.top(:,2) = th_specBC.top(:,1) - pi; % the second PF to be antiparallel
        else
            th_specBC.left = nan(1,2);
        end % if left boundary angle assigned
        
        % left
        if ~isempty(in.BCs_specs.right)
            th_specBC.right = zeros(in.Ny,1);
            for kk = 1:size(in.BCs_specs.right,1) % for number of rows
                indth = in.BCs_specs.right(kk,1):in.BCs_specs.right(kk,2);
                th_specBC.right(indth) = in.BCs_specs.right(kk,3)/180*pi;
%                 th_specBC.top(indth) = pi/2 + in.BCs_specs.top(kk,3)/180*pi;
            end
            assert(in.nOP ==2,'error AssignTh_SpecBC: currently works for 2 PFs only.')
            th_specBC.right(:,2) = pi - th_specBC.right(:,1); % the second PF to be antiparallel
%             th_specBC.top(:,2) = th_specBC.top(:,1) - pi; % the second PF to be antiparallel
        else
            th_specBC.right = nan(1,2);
        end % if right boundary angle assigned
        
    else % if BCs = Nspecial
        th_specBC = [];
        
    end % if Nspecial
end

%% PrepSpecialBCs
function bndry = PrepSpecialBCs(in)

    if strcmp(in.BCs,'Nspecial')
        bndry.top = logical( kron(ones(in.Nx,1),[zeros(in.Ny-1,1);1]) ); 
        bndry.bot = logical( kron(ones(in.Nx,1),[1; zeros(in.Ny-1,1)])); 
        bndry.left = logical( [ones(in.Ny,1) ; zeros(in.Ny*(in.Nx-1),1) ]); 
        bndry.right = logical( [zeros(in.Ny*(in.Nx-1),1) ; ones(in.Ny,1)]); 
    %     imagesc(reshape(bndry.bot,[in.Ny,in.Nx])), set(gca,'YDir','normal')
    %     imagesc(reshape(bndry.top,[in.Ny,in.Nx])), set(gca,'YDir','normal') 
    %     imagesc(reshape(bndry.left,[in.Ny,in.Nx])), set(gca,'YDir','normal') 
    %     imagesc(reshape(bndry.right,[in.Ny,in.Nx])), set(gca,'YDir','normal') 
        
        assert(strcmp(in.laplacianmethod,'9pt20')|strcmp(in.laplacianmethod,'5pt'), 'BCs = Nspecial but laplacianmethod ~= 9pt20 or 5pt.')
        if strcmp(in.laplacianmethod,'9pt20')
            % periodic BCs @ L and R
            bndry.matrixBT = sparse(toeplitz([4,1,zeros(1,in.Nx-3),1]));
            % Neumann BC @ L and R
            bndry.matrixBT(1,2) = 2; 
            bndry.matrixBT(in.Nx,in.Nx-1) = 2; 
            bndry.matrixBT(in.Nx,1) = 0; 
            bndry.matrixBT(1,in.Nx) = 0;
%             imagesc(bndry.matrixBT), colorbar
            
            % periodic BCs @ B and T
            bndry.matrixLR = sparse(toeplitz([4,1,zeros(1,in.Ny-3),1]));
            % Neumann BC @ B and T
            bndry.matrixLR(1,2) = 2; 
            bndry.matrixLR(in.Ny,in.Ny-1) = 2; 
            bndry.matrixLR(in.Ny,1) = 0; 
            bndry.matrixLR(1,in.Ny) = 0;
            % general Neumann BC @ B and T - to not compute laplacian in corners twice
            bndry.matrixLR(1,1) = 0; 
            bndry.matrixLR(1,2) = 0; 
            bndry.matrixLR(in.Ny,in.Ny) = 0;
            bndry.matrixLR(in.Ny,in.Ny-1) = 0; 
            bndry.matrixLR(in.Ny,1) = 0; 
            bndry.matrixLR(1,in.Ny) = 0;
%             imagesc(bndry.matrixLR), colorbar
            
        end
        
    else
        bndry = [];
    end % if Nspecial
end
%% AssignNormgrad_SpecBC
% calculates norm of the gradient at the top or bottom boundary based on
% the specification in the input
% assumes \gamma=1.5 => norm of the grad from values of PFs
function normgrad_i_specBC = AssignNormgrad_SpecBC(in,p_i,kappa,bndry,normgrad_i_specBC,i)
    if strcmp(in.BCs,'Nspecial')
        if ~isempty(in.BCs_specs.bottom)
            if all(size(kappa)==[1,1]) % if kappa only scalar
                normgrad_i_specBC.bot(:,i) = sqrt(2*in.m./kappa).*p_i(bndry.bot).*(1-p_i(bndry.bot));
            else
                normgrad_i_specBC.bot(:,i) = sqrt(2*in.m./kappa(bndry.bot)).*p_i(bndry.bot).*(1-p_i(bndry.bot));
            end % if kappa scalar

        else
            normgrad_i_specBC.bot = [];
        end % if bottom
    
        if ~isempty(in.BCs_specs.top)
            if all(size(kappa)==[1,1]) % if kappa only scalar
                normgrad_i_specBC.top(:,i) = sqrt(2*in.m./kappa).*p_i(bndry.top).*(1-p_i(bndry.top));
            else
                normgrad_i_specBC.top(:,i) = sqrt(2*in.m./kappa(bndry.top)).*p_i(bndry.top).*(1-p_i(bndry.top));
            end % if kappa scalar

        else
            normgrad_i_specBC.top = [];
        end% if top
        
        % if left
        if ~isempty(in.BCs_specs.left)
            if all(size(kappa)==[1,1]) % if kappa only scalar
                normgrad_i_specBC.left(:,i) = sqrt(2*in.m./kappa).*p_i(bndry.left).*(1-p_i(bndry.left));
            else
                normgrad_i_specBC.left(:,i) = sqrt(2*in.m./kappa(bndry.left)).*p_i(bndry.left).*(1-p_i(bndry.left));
            end % if kappa scalar

        else
            normgrad_i_specBC.left = [];
        end % if left
        
        % if right
        if ~isempty(in.BCs_specs.right)
            if all(size(kappa)==[1,1]) % if kappa only scalar
                normgrad_i_specBC.right(:,i) = sqrt(2*in.m./kappa).*p_i(bndry.right).*(1-p_i(bndry.right));
            else
                normgrad_i_specBC.right(:,i) = sqrt(2*in.m./kappa(bndry.right)).*p_i(bndry.right).*(1-p_i(bndry.right));
            end % if kappa scalar

        else
            normgrad_i_specBC.right = [];
        end % if right
        
    end % if BCs Nspecial
end % func

%% ApplySpecialBCs_lap
% normgrad_specBC_i ... column vector size = [Nx,1]
% stencil 5pt
%     stencil_coef = in.dx^2 ... divides ...  in.dy*cos(th_specBC_i).*(bndry.matrix*normgrad_specBC_i)
% stencil 9pt20
%    bndry.matrix*normgrad_specBC_i sums contribution from corresponding ghost nodes as by the stencil
%    stencil_coef = 6*in.dx^2 ... divides ...  in.dy*cos(th_specBC_i).*(bndry.matrix*normgrad_specBC_i)
function lapp_i = ApplySpecialBCs_lap(in,lapp_i,normgrad_specBC,bndry,th_specBC,i)
    
    % bottom
    if ~isempty(in.BCs_specs.bottom)
        normgrad_specBC_iloc = normgrad_specBC.bot(:,i);
        th_specBC_i = th_specBC.bot(:,i);
        
        switch in.laplacianmethod
            case '5pt'
                lapp_i(bndry.bot) = lapp_i(bndry.bot) - (1/in.dx)*2*cos(th_specBC_i).*normgrad_specBC_iloc;

            case '9pt20'
                lapp_i(bndry.bot) = lapp_i(bndry.bot) - (1/6/in.dx)*2*cos(th_specBC_i).*(bndry.matrixBT*normgrad_specBC_iloc);
%                 plot(lapp_i(bndry.bot),'.-'), hold on, plot(- (1/6/in.dy)*cos(th_specBC_i).*(bndry.matrix*normgrad_specBC_i),'.-'),hold off
        end
    end
    
    % top
    if ~isempty(in.BCs_specs.top)
        normgrad_specBC_iloc = normgrad_specBC.top(:,i);
        th_specBC_i = th_specBC.top(:,i);
        
        switch in.laplacianmethod
            case '5pt'
                lapp_i(bndry.top) = lapp_i(bndry.top) + (1/in.dx)*2*cos(th_specBC_i).*normgrad_specBC_iloc;

            case '9pt20'
                % bndry.matrix*normgrad_specBC_i sums contribution from corresponding ghost nodes as by the stencil
                lapp_i(bndry.top) = lapp_i(bndry.top) + (1/6/in.dx)*2*cos(th_specBC_i).*(bndry.matrixBT*normgrad_specBC_iloc);
        end
    end % if top
    
    % left
    if ~isempty(in.BCs_specs.left)
        normgrad_specBC_iloc = normgrad_specBC.left(:,i);
        th_specBC_i = th_specBC.left(:,i);
        
        switch in.laplacianmethod
            case '5pt'
                lapp_i(bndry.left) = lapp_i(bndry.left) - (1/in.dy)*2*cos(th_specBC_i).*normgrad_specBC_iloc;

            case '9pt20'
                % bndry.matrix*normgrad_specBC_i sums contribution from corresponding ghost nodes as by the stencil
                lapp_i(bndry.left) = lapp_i(bndry.left) - (1/6/in.dy)*2*cos(th_specBC_i).*(bndry.matrixLR*normgrad_specBC_iloc);
        end
    end % if left
    
    % right
    if ~isempty(in.BCs_specs.right)
        normgrad_specBC_iloc = normgrad_specBC.right(:,i);
        th_specBC_i = th_specBC.right(:,i);
        
        switch in.laplacianmethod
            case '5pt'
                lapp_i(bndry.right) = lapp_i(bndry.right) + (1/in.dy)*2*cos(th_specBC_i).*normgrad_specBC_iloc;

            case '9pt20'
                % bndry.matrix*normgrad_specBC_i sums contribution from corresponding ghost nodes as by the stencil
                lapp_i(bndry.right) = lapp_i(bndry.right) + (1/6/in.dy)*2*cos(th_specBC_i).*(bndry.matrixLR*normgrad_specBC_iloc);
        end
    end % if left
    
end% func

%% ApplySpecialBCs_grad
function gp_i = ApplySpecialBCs_grad(in,gp_i,normgrad_specBC,bndry,th_specBC,i)
%     xdirsign = [1,1];
    
    if ~isempty(in.BCs_specs.bottom)
        normgrad_specBC_iloc = normgrad_specBC.bot(:,i);
        th_specBC_i = th_specBC.bot(:,i);
        xdirsign = [-1 , 1];
        
        gp_i{1}(bndry.bot) =  xdirsign(i)*sin(th_specBC_i).*normgrad_specBC_iloc;
        gp_i{2}(bndry.bot) = cos(th_specBC_i).*normgrad_specBC_iloc;
    end
    
    if ~isempty(in.BCs_specs.top)
        normgrad_specBC_iloc = normgrad_specBC.top(:,i);
        th_specBC_i = th_specBC.top(:,i);
        xdirsign = [1 , -1];
        
        gp_i{1}(bndry.top) =  xdirsign(i)*sin(th_specBC_i).*normgrad_specBC_iloc;
        gp_i{2}(bndry.top) = cos(th_specBC_i).*normgrad_specBC_iloc;
    end
    
    if ~isempty(in.BCs_specs.left)
        normgrad_specBC_iloc = normgrad_specBC.left(:,i);
        th_specBC_i = th_specBC.left(:,i);
        ydirsign = [1 , -1];
        
        gp_i{1}(bndry.left) =  cos(th_specBC_i).*normgrad_specBC_iloc;
        gp_i{2}(bndry.left) = ydirsign(i)*sin(th_specBC_i).*normgrad_specBC_iloc;
    end
    
    if ~isempty(in.BCs_specs.right)
        normgrad_specBC_iloc = normgrad_specBC.right(:,i);
        th_specBC_i = th_specBC.right(:,i);
        ydirsign = [1 , -1];
        
        gp_i{1}(bndry.right) =  cos(th_specBC_i).*normgrad_specBC_iloc;
        gp_i{2}(bndry.right) = ydirsign(i)*sin(th_specBC_i).*normgrad_specBC_iloc;
    end
    
end % func



%% Calc_p_derived_fields_indiv_analytically
function [gradXXp, gradYXp, gradYYp,out_gradsumgradpsq] = Calc_p_derived_fields_indiv_analytically(single_p , single_gp , in_gradsumgradpsq)
    p = single_p;

    global gradXX gradYX gradYY
    gradXXp = gradXX(p);
    gradYXp = gradYX(p);
    gradYYp = gradYY(p);
    out_gradsumgradpsq{1} =  in_gradsumgradpsq{1} + 2*( single_gp{1}.*gradXXp + single_gp{2}.*gradYXp ) ; 
    out_gradsumgradpsq{2} =  in_gradsumgradpsq{2} + 2*( single_gp{1}.*gradYXp + single_gp{2}.*gradYYp ) ; 
end % func Calc_p_derived_fields_indiv_analytically

% %% Calc_p_derived_fields_indiv_numerically
% function [gradXXp, gradYXp, gradYYp,out_gradsumgradpsq] = Calc_p_derived_fields_indiv_analytically(single_p , single_gp , in_gradsumgradpsq)
%     p = single_p;
%     gp = single_gp;
% 
%     global gradXX gradYX gradYY
%     gradXXp = gradXX(p);
%     gradYXp = gradYX(p);
%     gradYYp = gradYY(p);
%     out_gradsumgradpsq{1} =  in_gradsumgradpsq{1} + 2*( gp{1}.*gradXXp + gp{2}.*gradYXp ) ; 
%     out_gradsumgradpsq{2} =  in_gradsumgradpsq{2} + 2*( gp{1}.*gradYXp + gp{2}.*gradYYp ) ; 
% end % func Calc_p_derived_fields_indiv_analytically
%% CalcMagnitudeInterfaceNormal_ij
function [magnitude_of_normal_squared , magnitude_of_normal] = CalcMagnitudeInterfaceNormal_ij(gradientPFi,gradientPFj)
    gpi = gradientPFi;
    gpj = gradientPFj;
    magnitude_of_normal_squared = (gpi{1}-gpj{1}).^2 + (gpi{2}-gpj{2}).^2;
    magnitude_of_normal = sqrt( magnitude_of_normal_squared  ) ; 
end

%% DetermineInnerOuterInterfaceRegion_ij
function [out_interf_ij , inn_interf_ij] = DetermineInnerOuterInterfaceRegion_ij(is_isotropic_or_precycle_ij, intf_normal_magnitude_ij, intf_limval,basic_plotparams)
%     limitvals.iso = 1e-12; 
    dx = basic_plotparams.spacing(1);
     if is_isotropic_or_precycle_ij
         limval = intf_limval.iso; % difference in value of normal field 
%             out_interf_ij = true(simsize);
%             interface = true(simsize);
        out_interf_ij = intf_normal_magnitude_ij > limval/dx;
        inn_interf_ij = intf_normal_magnitude_ij > limval/dx;
%             interface = nnl{k}>filterGradVal/dx;
     else
        outerlim = intf_limval.aniso.outer; 
        innerlim = intf_limval.aniso.inner;
        out_interf_ij = intf_normal_magnitude_ij.^3 > (outerlim/dx)^3;
        inn_interf_ij = intf_normal_magnitude_ij.^3 > (innerlim/dx)^3;
%             interface =  sumpsqpsq>1e-12;
     end
     is_inner_inside_outer = ~(inn_interf_ij&~out_interf_ij); % all ones when interface lays in out_interf_ij
     tstep = basic_plotparams.tstep;
     assert(all(is_inner_inside_outer),['Error in Multip_aniso_conc>DetermineInnerOuterInterfaceRegion_ij @tstep = ' num2str(tstep) ': ''inner'' interface region is not subset of ''outer''.'])
     
%      plotting = basic_plotparams.plotting;
%      size2D = basic_plotparams.size2D;
%      if plotting 
%         figure(13)
%         imagesc(reshape(~inn_interf_ij&out_interf_ij,size2D)), 
%         title([num2str(tstep) 'Interface: in outer & out of inner']),
%         set(gca,'YDir','normal'),
%         daspect([1 1 1])
%      end
    
end % func DetermineInnerOuterInterfaceRegion_ij

%% Calc_p_derived_fields_npairs
function [sumpsqpsq, gradsumpsqpsq,gradsumgradpsq] = Calc_p_derived_fields_npairs(i,p,gp,sumpsqpsq,gradsumpsqpsq,gradsumgradpsq)
    global gradXX gradYY gradYX
    for j = (i-1):-1:1 % cycle runs from i-1 down to 1 because these gp were already calculated 
        assert(i>j,'i>j in Calc_p_derived_fields_npairs')
        sumpsqpsq = sumpsqpsq + (p{i}.^2).*(p{j}.^2);
        gradsumpsqpsq{1} = gradsumpsqpsq{1} + 2*p{i}.*gp{i}{1}.*(p{j}.^2) + (p{i}.^2).*2.*p{j}.*gp{j}{1};
        gradsumpsqpsq{2} = gradsumpsqpsq{2}+ 2*p{i}.*gp{i}{2}.*(p{j}.^2) + (p{i}.^2).*2.*p{j}.*gp{j}{2};
        gradsumgradpsq{1} = gradsumgradpsq{1} + 2*gp{i}{1}.*gradXX(p{i}) + 2*gp{i}{2}.*gradYX(p{i});
        gradsumgradpsq{2} = gradsumgradpsq{2}+ 2*gp{i}{1}.*gradYX(p{i}) + 2*gp{i}{2}.*gradYY(p{i});
    end
end % func Calc_p_derived_fields_npairs

%% Calc_sumpsqpsq
function sumpsqpsq = Calc_sumpsqpsq(i,p,sumpsqpsq)
    for j = (i-1):-1:1 % cycle runs from i-1 down to 1 because these gp were already calculated 
        sumpsqpsq = sumpsqpsq + (p{i}.^2).*(p{j}.^2);
    end
end % Calc_sumpsqpsq
%% CalcGradientOfOrientationField
% indpair = indpairs(k,:)
% intf_norm_mag_sq_ij = nnlsq{k}
% intf_norm_mag_ij = nnl{k}
% intf_norm_ij = nn{k}
function [dnndX_ij,dnndY_ij,doridX_ij,doridY_ij] = CalcGradientOfOrientationField(indpair,out_interf_ij,intf_norm_mag_sq_ij,intf_norm_mag_ij,intf_norm_ij,gradXXp,gradYXp,gradYYp)
    i = indpair(1);
    j = indpair(2);
    systemsize = size(out_interf_ij);
    
    dnndX_ij{1} = zeros(systemsize); % d intf_norm_ij_x /d x = d phi_x/dx
    dnndX_ij{2} = zeros(systemsize); % d intf_norm_ij_y/d x = d phi_y/dx
    dnndY_ij{1} = zeros(systemsize); % d intf_norm_ij_x/d y = d phi_x/dy
    dnndY_ij{2} = zeros(systemsize); % d intf_norm_ij_y/d y = d phi_y/dy
    
    global gradX gradY
    gradXnnlsq = gradX(intf_norm_mag_sq_ij);
    gradYnnlsq  = gradY(intf_norm_mag_sq_ij);
    
    dnndX_ij{1}(out_interf_ij) = ( gradXXp{i}(out_interf_ij) - gradXXp{j}(out_interf_ij) )./intf_norm_mag_ij(out_interf_ij)...
                                            - intf_norm_ij{1}(out_interf_ij).*gradXnnlsq(out_interf_ij)./2./intf_norm_mag_sq_ij(out_interf_ij); % d phi_{i,j}_x d x 
    dnndX_ij{2}(out_interf_ij) = (gradYXp{i}(out_interf_ij) - gradYXp{j}(out_interf_ij))./intf_norm_mag_ij(out_interf_ij)...
                                            - intf_norm_ij{2}(out_interf_ij).*gradXnnlsq(out_interf_ij)./2./intf_norm_mag_sq_ij(out_interf_ij) ;% d phi_{i,j}_y d x
    dnndY_ij{1}(out_interf_ij) = (gradYXp{i}(out_interf_ij) - gradYXp{j}(out_interf_ij))./intf_norm_mag_ij(out_interf_ij)  ...
                                            - intf_norm_ij{1}(out_interf_ij).*gradYnnlsq(out_interf_ij)./2./intf_norm_mag_sq_ij(out_interf_ij) ;% d phi_{i,j}_x d y
    dnndY_ij{2}(out_interf_ij) = (gradYYp{i}(out_interf_ij) - gradYYp{j}(out_interf_ij) )./intf_norm_mag_ij(out_interf_ij)...
                                            - intf_norm_ij{2}(out_interf_ij).*gradYnnlsq(out_interf_ij)./2./intf_norm_mag_sq_ij(out_interf_ij);  % d phi_{i,j}_y d y 
    doridX_ij = zeros(systemsize);
    doridY_ij = zeros(systemsize);
    doridX_ij(out_interf_ij) = (intf_norm_ij{1}(out_interf_ij).*dnndX_ij{2}(out_interf_ij) - intf_norm_ij{2}(out_interf_ij).*dnndX_ij{1}(out_interf_ij)); 
    doridY_ij(out_interf_ij) = (intf_norm_ij{1}(out_interf_ij).*dnndY_ij{2}(out_interf_ij) - intf_norm_ij{2}(out_interf_ij).*dnndY_ij{1}(out_interf_ij));
end % func

%% DispMsgPrecycleTerminated
function DispMsgPrecycleTerminated(precycle,tstep,is_inclination_dependent)
    if (tstep == precycle) && is_inclination_dependent
        disp(['time step ' num2str(tstep) ': Precycle terminated, anisotropic simulation starts.'])
    end
end % func DispMsgPrecycleTerminated

%% FindRowsWith_i_And_j_Larger
function ind_rows_with_i_and_j_larger = FindRowsWith_i_And_j_Larger(indpairs,i)
    j_larger_than_i = max(indpairs,[],2) > i; % bool, all rows with j>i
    rows_with_i = any(indpairs==i,2);
    ind_rows_with_i_and_j_larger =  find(j_larger_than_i & rows_with_i); % find indices of rows containing i where j>i
%     if isempty(ind_rows_with_i_and_j_larger)
%         ind_rows_with_i_and_j_larger = 0;
%     end
end % func FindRowsWith_i_And_j_Larger

% % test, seems to work fine
% nOP = 5; % 10 interfaces
% indpairs = combnk(1:nOP,2); % 10 combinations
% for i = 1:nOP
%     j_larger_than_i = max(indpairs,[],2) > i; % bool, all rows with j>i
%     ind_rows_with_i_and_j_larger =  find(any(indpairs(j_larger_than_i)==i,2)); % find indices of rows containing i where j>i
%     disp(['OP ' num2str(i)])
%     indpairs(ind_rows_with_i_and_j_larger,:)
% end

%% CalcKapGam_iso
function gkL = CalcKapGam_iso(in,k,gkL)
% global simsize
    if strcmp(in.model,'IWc')
        gkL.k_ij = in.kpp0(k); % 
        gkL.g_ij{k} = in.gam0(k);
        
    elseif strcmp(in.model,'IWvG')
        gkL.g_ij{k} = in.gam0(k);
        
    elseif strcmp(in.model,'IWvK')
        gkL.k_ij = in.kpp0(k); % 
    
    end % if model
    
    gkL.L_ij = in.Lij(k);
end% func
%% CalcKapGam_aniso
function [gkL, aniso] = CalcKapGam_aniso(in,modspec, lim_forbb_ang,k, ori_k,gkL, aniso)
            
            cond_ori = ~isnan(ori_k);
            % interface stiffness of f(th) = sig_0(1+d*cos(n*th)) is intf_stiffness = sig_0(1-(n^2-1)*d*cos(n*th))
%             intf_stiffness = in.intf.IE_phases*(1 - (in.intf.params_incl_dep.nfold^2-1)*in.intf.params_incl_dep.soaIE*cos(in.intf.params_incl_dep.nfold*ori_k)); % 
%             th = linspace(-pi/4,pi/4,100);
%             intf_stiffness = in.intf.IE_phases*(1 - (in.intf.params_incl_dep.nfold^2-1)*in.intf.params_incl_dep.soaIE*cos(in.intf.params_incl_dep.nfold*th)); % 
%             figure
%             plot(th,intf_stiffness,'o',th,in.Lfun(th),'x')
            if in.is_inclination_dependent_IE
                if modspec.kpp
                    [aniso.KP{k}, aniso.KPder{k}, aniso.KPdder{k}] = CalcAnisotropyFunc(in,lim_forbb_ang,k,ori_k,in.intf.params_incl_dep.soaIE);

                    if strcmp(in.model,'IWc')
                        gkL.k_ij = in.kpp0(k)*aniso.KP{k};
                        L_anisofactor_IE = 1; % Lcorr
%                         L_anisofactor = aniso.KP{k}(cond_ori); % LcorrEmpiric

                    elseif strcmp(in.model,'IWvK')
                        gkL.k_ij = in.kpp0(k)*aniso.KP{k}.^2;
                        gkL.g_ij = 1.5;
                        [gkL.k_ij, ~] = ReplaceNaNs_kpp_gmm_ij(in.intf.isStrongAniso, gkL.k_ij,in.kpp0(k),0,0);
                        % !!!
                        if ~in.intf.params_incl_dep.override_IW_Lcorr % if Lcorr from variable IW to be applied
                            L_anisofactor_IE = 1./aniso.KP{k}(cond_ori); % Lcorr
    %                         L_anisofactor = 1; % LcorrEmpiric
                        else
                            L_anisofactor_IE  = 1;
                        end

                    end % if IWc or IWvK
                end % if 

                if modspec.gam
                    if in.intf.isStrongAniso 
                        [aniso.GM{k}, aniso.GMder{k}, aniso.GMdder{k}] = CalcAnisotropyFunc(in,lim_forbb_ang,k,ori_k,in.soaGM);
                        gkL.g_ij{k} = CalcGamma_ij(in.intf.isStrongAniso, in.model, in.gam0(k), aniso.GM{k});
                        [gkL.k_ij, gkL.g_ij{k}] = ReplaceNaNs_kpp_gmm_ij(in.intf.isStrongAniso,gkL.k_ij,in.kpp0(k),gkL.g_ij{k},in.gam0(k));
        %                 [kpp_ij, ~] = ReplaceNaNs_kpp_gmm_ij(in.intf.isStrongAniso, kpp_ij,in.kpp0(k),0,0);
                    else  % weak anisotropy
                        if modspec.kpp % IWc
                            % anisoKP is = anisoIE because soaKP=soaIE
                            % CalcGamma_ij distinguishes IWc and IWvG
                            gkL.g_ij{k} = CalcGamma_ij(in.intf.isStrongAniso, in.model, in.A(k), aniso.KP{k});
                            [gkL.k_ij, gkL.g_ij{k}] = ReplaceNaNs_kpp_gmm_ij(in.intf.isStrongAniso, gkL.k_ij,in.kpp0(k),gkL.g_ij{k},in.A(k));
                        else % IWvG
                            % anisofun not calculated yet
                            [aniso.IE{k}, aniso.IEder{k}, aniso.IEdder{k}] = CalcAnisotropyFunc(in,lim_forbb_ang,k,ori_k,in.intf.params_incl_dep.soaIE);
                            % CalcGamma_ij distinguishes IWc and IWvG
                            gkL.g_ij{k} = CalcGamma_ij(in.intf.isStrongAniso, in.model, in.A(k), aniso.IE{k});
                            gkL.k_ij = in.kpp0;
                            [~, gkL.g_ij{k}] = ReplaceNaNs_kpp_gmm_ij(in.intf.isStrongAniso, 0,0,gkL.g_ij{k},in.A(k));
                            % !!!
                            if ~in.intf.params_incl_dep.override_IW_Lcorr % if Lcorr from variable IW to be applied
                                L_anisofactor_IE = aniso.IE{k}(cond_ori); % Lcorr
    %                             L_anisofactor = aniso.KP{k}(cond_ori).^2; % LcorrEmpiric
                            else
                                L_anisofactor_IE  = 1;
                            end

        %                 [kpp_ij, ~] = ReplaceNaNs_kpp_gmm_ij(in.intf.isStrongAniso, kpp_ij,in.kpp0(k),0,0);
                        end % if IWc
                    end % strong vs weak
                end % if IWc or IWvG
            else
                L_anisofactor_IE = 1;
            end % is_inclination_dependent_IE
            
            gkL.L_ij = in.Lij(k)*ones(size(ori_k));
            if in.is_locally_aniso_L(k)
                L_anisofactor_L = in.intf.params_incl_dep.Lfun(ori_k(cond_ori));
            else
                L_anisofactor_L = 1;
            end
            
            gkL.L_ij(cond_ori) = in.Lij(k).*L_anisofactor_L.*L_anisofactor_IE;
            
end % func

%% CalcAnisotropyFunc
function [anisoPAR, danisoPAR, ddanisoPAR] = CalcAnisotropyFunc(in,lim_forbb_ang,k,ori,soa)
    
    % workaround to use functions dealing with anisotropy based on older input structures
    if ~isfield(in.intf.params_incl_dep,'offset_ang') % if offset_ang is not one of fields in params_incl_dep
            in.intf.params_incl_dep.offset_ang = in.misori(k); % implying new input including misorientation
    end
    
    if in.intf.params_incl_dep.Omega <=1
        offset_ang_included = true;
        [fIEijfun,dfIEijfun,ddfIEijfun] = AssignAnisotropyFunction(in.intf.params_incl_dep,offset_ang_included);
        anisoPAR = fIEijfun(ori,soa);
        danisoPAR = dfIEijfun(ori,soa);
    else % Omega>1
%         have th_m global and function [anisoPAR, danisoPAR]
        assert(~in.all_expr_anal,'all_expr_anal = true BUT function ''calc_strong_reg_anisofun'' does not return 2nd derivative needed for analytic expressions.')
        [anisoPAR , danisoPAR] = calc_strong_reg_anisofun(ori, lim_forbb_ang,in.intf.params_incl_dep,soa,false);
    end
    
    if in.all_expr_anal
        ddanisoPAR = ddfIEijfun(ori,soa);
    else
        ddanisoPAR = NaN;
    end
end % func CalcAnisotropyKappaOrStrongGamma
%% CalcGamma_ij
function gamma_ij = CalcGamma_ij(isStrongAniso, model, multiplier, anisofunc)

    if isStrongAniso
        gamma_ij = multiplier*anisofunc;
        
    else % weak anisotropy
        % multiplier = 9/2g(\gamma)^2 = A(k)
        if strcmp(model,'IWc')
            % anisoKP is = anisoIE because soaKP=soaIE
            % g_ij = (-A*anisoIE/2 -1)/(A*anisoIE -2 )
            gamma_ij = (-multiplier*anisofunc/2-1)./(multiplier*anisofunc-2); % weak anisotropy, formula (51) from (Moelans, 2008)
            
        elseif strcmp(model,'IWvG')
            % g_ij = (-A*anisoIE^2/2 -1)/(A*anisoIE^2 -2 )
            gamma_ij = (-multiplier*anisofunc.^2/2-1)./(multiplier*anisofunc.^2-2); % weak anisotropy, MODIFIED formula (51) from (Moelans, 2008)
            
        end %if model
        
    end % if isStrongAniso
    
end % func CalcGamma_ij
%% ReplaceNaNs_kpp_gmm_ij
function [kappa_ij, gamma_ij] = ReplaceNaNs_kpp_gmm_ij(isStrongAniso,kappa_ij,kpp0,gamma_ij,gamma_multip)
    kappa_ij(isnan(kappa_ij)) = kpp0;
    if isStrongAniso
        % gamma_multip = gam0(k) for strong aniso
%         gamij{k}(isnan(gamij{k})) = gam0(k);
        gamma_ij(isnan(gamma_ij)) = gamma_multip; 
    else
        % gamma_multip = A(k) for strong aniso
%         gamij{k}(isnan(gamij{k})) = (-A(k)/2-1)./(A(k)-2); % for gamma isotropic
        gamma_ij(isnan(gamma_ij)) = (-gamma_multip/2-1)./(gamma_multip-2); 
    end
end% func ReplaceNaNs_kpp_gmm_ij
%% CalcKappaGammaL_Add_kth_term
function [kappa, gam, L] = CalcKappaGammaL_Add_kth_Term(kappa, gam, L, gkL, in, p, i, j, k)
% [kappa, gam_sumpsqpsq, L] = CalcKappaGammaL_Add_kth_Term(kappa, gam_sumpsqpsq, L, gk, in, p, i, j, k);
    
    if in.nOP == 2
        L = gkL.L_ij; % in.Lij(k);
    else 
        psqpsq_ij = p{i}.^2.*p{j}.^2 ;
        L = L + gkL.L_ij.*psqpsq_ij;
%         L = L + in.Lij(k).*psqpsq_ij;
    end% if nOP
        
    if strcmp(in.model,'IWc')
        kappa_ij = gkL.k_ij;
        gamma_ij = gkL.g_ij{k};
        
        if in.nOP == 2
            kappa = kappa_ij;
            gam   = gamma_ij;
        else 
            kappa = kappa + kappa_ij.*psqpsq_ij;
            gam   = gam   + gamma_ij.*psqpsq_ij;
        end % if nOP
        
    elseif strcmp(in.model,'IWvG')
        kappa = in.kpp0;
        gamma_ij = gkL.g_ij{k};
        
        if in.nOP == 2
            gam   = gamma_ij;
        else 
            gam   = gam   + gamma_ij.*psqpsq_ij;
        end % if nOP
        
    elseif strcmp(in.model,'IWvK')
        gam   = 1.5;
        kappa_ij = gkL.k_ij;
        
        if in.nOP == 2
            kappa = kappa_ij;
        else 
            kappa = kappa + kappa_ij.*psqpsq_ij;
        end % if nOP
        
    end % if model
    
end % func 
%% CalcSignFromOrderInPair
% prefactor is (-1) when 'i' is in the 2nd column of k-th row of indpairs
% and (1) when it is in the 1st column
function prefactor = CalcSignFromOrderInPair(k,i,indpairs) 
    prefactor = sign(1.5-find(indpairs(k,:)==i)); 
end

%% CalcGderPrefactors
% calculates model-specific factors for gradient derivatives
% Gdermultips ... struct with fields 'kpp', 'gam'
%   - the factors are calculated distinguishing between the 'weak' and
%   'strong' anisotropy approximation (affects way of calculating anisotropy of gamma)
function Gdermultips_ij = CalcGderPrefactors(in,k,aniso)
    
    if strcmp(in.model,'IWc')
        Gdermultips_ij.kpp = in.kpp0(k).*aniso.KPder{k};
        if in.intf.isStrongAniso
            Gdermultips_ij.gam = in.gam0(k).*aniso.GMder{k};
        else
            % IWvG ... d gam/d anisofun = d gam / d a = 2A/(Aa-2)^2
            A = in.A(k);
            Gdermultips_ij.gam = ( 2*A./(A.*aniso.KP{k}-2).^2 ).*aniso.KPder{k}; % for weak aniso: anisoKPder=anisoIEder
        end % if strong aniso
        
    elseif strcmp(in.model,'IWvG')
        Gdermultips_ij.kpp = 0;
        if in.intf.isStrongAniso
            Gdermultips_ij.gam = in.gam0(k).*aniso.GMder{k};
        else
            % IWvG ... d gam/d anisofun = d gam / d a = 4Aa/(Aa^2-2)^2
            A = in.A(k);
            Gdermultips_ij.gam = ( 4*A*aniso.IE{k}./(A.*aniso.IE{k}.^2-2).^2 ).*aniso.IEder{k};
        end % if strong aniso
        
    elseif strcmp(in.model,'IWvK')
        Gdermultips_ij.kpp = 2*in.kpp0(k)*aniso.KP{k}.*aniso.KPder{k};
        Gdermultips_ij.gam = 0;
        
    end % if model
    
    %     remove NaNs originated from angle
    Gdermultips_ij.gam(isnan(Gdermultips_ij.gam)) = 0;
    Gdermultips_ij.kpp(isnan(Gdermultips_ij.kpp)) = 0;
end % func

%% CalcGderivativeOfKappaGamma_ij
% CAREFUL: confusing use of the last in variable daniso_gam (to have strong/weak approx in one function)
function [dkppijdGp_ij,dgmmijdGp_ij] = CalcGderivativeOfKappaGamma_ij(condOuterInterace,nn,nnl,Gdermultips_ij,modspec)
    nn_x = nn{1}; 
    nn_y = nn{2};
    simsize = size(condOuterInterace);
    dxxijdGp{1} = zeros(simsize); 
    dxxijdGp{2} = zeros(simsize); 
    dxxijdGp{1}(condOuterInterace) = -nn_y(condOuterInterace)./nnl(condOuterInterace);
    dxxijdGp{2}(condOuterInterace) =  nn_x(condOuterInterace)./nnl(condOuterInterace);
%     dkppijdGp_ij{1} = kpp0*dxxijdGp{1}.*danisoKP; 
%     dkppijdGp_ij{2} = kpp0*dxxijdGp{2}.*danisoKP; 
    
    for ddim=1:2
        if modspec.kpp
            dkppijdGp_ij{ddim} = dxxijdGp{ddim}.*Gdermultips_ij.kpp; 
        else
            dkppijdGp_ij{ddim} = 0;
        end
        
        if modspec.gam
            dgmmijdGp_ij{ddim} = dxxijdGp{ddim}.*Gdermultips_ij.gam;
        else
            dgmmijdGp_ij{ddim} = 0;
        end
    end
%     if isStrongAniso
%         dgmmijdGp_ij{1} = multip_gam*dxxijdGp{1}.*daniso_gam; 
%         dgmmijdGp_ij{2} = multip_gam*dxxijdGp{2}.*daniso_gam; 
%     else
%         danisoIE = danisoKP;
%         anisoIE = daniso_gam;
%         dgmmijdGp_ij{1} = dxxijdGp{1}.*( 2*multip_gam./(multip_gam.*anisoIE-2).^2 ).*danisoIE; 
%         dgmmijdGp_ij{2} = dxxijdGp{2}.*( 2*multip_gam./(multip_gam.*anisoIE-2).^2 ).*danisoIE; 
% %         dgmmijdGp_ij{1} = dxxijdGp{1}.*( 2*A./(A.*anisoKP{k}-2).^2 ).*danisoKP{k}; 
% %         dgmmijdGp_ij{2} = dxxijdGp{2}.*( 2*A./(A.*anisoKP{k}-2).^2 ).*danisoKP{k}; 
%     end
                
end % func CalcGderivativeOfKappaGamma_ij
%% GetIndexOfTheOther
function j = GetIndexOfTheOther(k,indpairs,i)
    j = indpairs(k,indpairs(k,:)~=i); % index of the other OP in the pair k
end
%% CalcGderivativeOfKappa_i_And_j_AddTerm
% function [dkppdGp_all,dgmmdGp_all] = CalcGderivativeOfKappaGamma_i_And_j_AddTerm(dkppdGp_all,dgmmdGp_all,dkppijdGp_ij,dgmmijdGp_ij, p, nOP, i, j)
function dkppdGp_all = CalcGderivativeOfKappa_i_And_j_AddTerm(dkppdGp_all,dkppijdGp_ij, p, nOP, i, j)
% dkppdGp = CalcGderivativeOfKappa_i_And_j_AddTerm(dkppdGp,dkppijdGp_ij, p, in.nOP, i, j);
    if nOP == 2
        for ddim = 1:2
            dkppdGp_all{i}{ddim}    = dkppijdGp_ij{ddim};
%             dgmmdGp_all{i}{ddim}  = dgmmijdGp_ij{ddim}; % zeros(simsize);
            dkppdGp_all{j}{ddim}    = - dkppijdGp_ij{ddim};
%             dgmmdGp_all{j}{ddim}  = - dgmmijdGp_ij{ddim}; % zeros(simsize);
        end % for ddim
    else % nOP >2 
        psqpsq_ij = p{i}.^2.*p{j}.^2 ;
        ndim = 2;
        for ddim = 1:ndim
            dkppdGp_all{i}{ddim} = dkppdGp_all{i}{ddim} + dkppijdGp_ij{ddim}.*psqpsq_ij;
%             dgmmdGp_all{i}{ddim} = dgmmdGp_all{i}{ddim} + dgmmijdGp_ij{ddim}.*psqpsq_ij;
            
            dkppdGp_all{j}{ddim} = dkppdGp_all{j}{ddim} - dkppijdGp_ij{ddim}.*psqpsq_ij;
%             dgmmdGp_all{j}{ddim} = dgmmdGp_all{j}{ddim} - dgmmijdGp_ij{ddim}.*psqpsq_ij;
        end % for ddim
    end% if nOP
end % func CalcGderivativeOfKappaGamma_i_And_j_AddTerm

%% CalcDrivingForceTerms_Isotropic
% function [df0dp , DFlapp] = CalcDrivingForceTerms_Isotropic(p_i, m, sumgam_ij_psq_i, kappa,lapp_i,condIntoDF)
function [df0dp , DFlapp] = CalcDrivingForceTerms_Isotropic(p_i, m, sumgam_ij_psq_i, kappa,lapp_i)
    p = p_i;
%     df0dp = m*(p.^3 - p + 2*gam.*p.*(sumpsq - p.^2));
    df0dp = m*(p.^3 - p + 2*p.*sumgam_ij_psq_i);
    DFlapp = kappa.*lapp_i;
    
%     figure(2)
%     subplot(1,2,1),imagesc(reshape(df0dp,[100,50])), axis equal, set(gca,'YDir','normal'),colorbar
%     subplot(1,2,2),imagesc(reshape(DFlapp,[100,50])), axis equal, set(gca,'YDir','normal'),colorbar
    
end % func CalcDrivingForceTerms_Isotropic
%% CalcDrivingForceTerms_Anisotropic
% function [DFdivG_i , DFvG_i , DFdivK_i ,  DFvK_i , DFgK_i] = CalcDrivingForceTerms_Anisotropic(condIntoDF,m,p_i, sumdgamdGp_ij_vec_i,sumdgamdGp_ij_vecdotpr_i, sumdgamdGp_ij_div_i, ...
%         sumgradpsq_mod,divdkppdGp_mod_i, grad_sumgrad_over_sumpsqpsq,dkppdGp_sumpsqpsq_i, gp_i,gradkpp)
function [DFdivG_i , DFvG_i , DFdivK_i ,  DFvK_i , DFgK_i] = CalcDrivingForceTerms_Anisotropic(condIntoDF,m, p_i, gp_i, DFin,modspec,i)
    simsize = size(condIntoDF);            
    if modspec.gam
        DFdivG_i = zeros(simsize);
        DFvG_i = zeros(simsize);
        DFdivG_i(condIntoDF) = m*p_i(condIntoDF).^2.*DFin.sumdgamdGp_ij_div_i(condIntoDF);
         % first term of DFvG_i as dot product
        DFvG_i(condIntoDF) = p_i(condIntoDF).*( gp_i{1}(condIntoDF).*DFin.sumdgamdGp_ij_vec_i{1}(condIntoDF) + gp_i{2}(condIntoDF).*DFin.sumdgamdGp_ij_vec_i{2}(condIntoDF) );
        % second term of DFvG_i as element-wise product
        DFvG_i(condIntoDF) = DFvG_i(condIntoDF) + p_i(condIntoDF).^2.*DFin.sumdgamdGp_ij_vecdotpr_i(condIntoDF);
        DFvG_i(condIntoDF) = 2*m*DFvG_i(condIntoDF);
    else
        DFdivG_i = 0;
        DFvG_i = 0;
    end
    
    if modspec.kpp
        DFdivK_i = zeros(simsize);
        DFvK_i = zeros(simsize);
        DFgK_i = zeros(simsize);
        DFdivK_i(condIntoDF) = 0.5*DFin.divdkppdGp_mod{i}(condIntoDF).*DFin.sumgradpsq_mod(condIntoDF);
        DFvK_i(condIntoDF) = 0.5*DFin.dkppdGp_sumpsqpsq{i}{1}(condIntoDF).*DFin.grad_sumgrad_over_sumpsqpsq{1}(condIntoDF) + ... 
                                        0.5*DFin.dkppdGp_sumpsqpsq{i}{2}(condIntoDF).*DFin.grad_sumgrad_over_sumpsqpsq{2}(condIntoDF);
        DFgK_i(condIntoDF) = DFin.gradkpp{1}(condIntoDF).*gp_i{1}(condIntoDF) + DFin.gradkpp{2}(condIntoDF).*gp_i{2}(condIntoDF);
    else
        DFdivK_i = 0;
        DFvK_i = 0;
        DFgK_i = 0;
    end
    
%     DFdivK_i(condIntoDF) = 0.5*divdkppdGp_mod_i(condIntoDF).*sumgradpsq_mod(condIntoDF);
%     DFvK_i(condIntoDF) = 0.5*dkppdGp_sumpsqpsq_i{1}(condIntoDF).*grad_sumgrad_over_sumpsqpsq{1}(condIntoDF) + 0.5*dkppdGp_sumpsqpsq_i{2}(condIntoDF).*grad_sumgrad_over_sumpsqpsq{2}(condIntoDF);
%     DFgK_i(condIntoDF) = gradkpp{1}(condIntoDF).*gp_i{1}(condIntoDF) + gradkpp{2}(condIntoDF).*gp_i{2}(condIntoDF);
    
end % func CalcDrivingForceTerms_Anisotropic
%% dinterpf
function dinterpf_i = dinterpf(p,i,sumpsq)
    dinterpf_i = 2*p{i}.*(sumpsq-p{i}.^2)./sumpsq.^2;
end % func
% %% CalcEnergyDensity_AddTerm
% function Fdens = CalcHomogEnergyDensity_AddTerm(Fdens,p,i,m,gam)
%     Fdens = Fdens + m*(0.25*p{i}.^4 - 0.5*p{i}.^2);
%     for j = (i-1):-1:1
%         Fdens = Fdens + m*gam.*(p{i}.^2).*(p{j}.^2);
%     end
% end% func CalcEnergyDensity_AddTerm
%% CheckAnistropicTermsAreZero
function CheckAnistropicTermsAreZero(DFdivG_i,DFvG_i,DFdivK_i,DFvK_i,DFgK_i,i,tstep)
    DFcheck(1) = ~all(all( DFdivG_i ==0 ));
    DFcheck(2) = ~all(all( DFvG_i==0 ));
    DFcheck(3) = ~all(all( DFdivK_i==0 ));
    DFcheck(4) = ~all(all( DFvK_i==0 ));
    DFcheck(5) = ~all(all( DFgK_i==0 ));
    if any(DFcheck) 
        error(['NONZERO NON-ISOTROPIC term(s): [' num2str(find(DFcheck)) ']  ... OP ' num2str(i) ', tstep ' num2str(tstep)])
    end 
end % func CheckAnistropicTermsAreZero

%% CheckPValuesInInterval
function [pvalUnder_i,pvalAbove_i,is_pvalUnder_i,is_pvalAbove_i] = CheckPValuesInInterval(p_i,pValTol)
        pvalUnder_i = p_i < -pValTol;
        pvalAbove_i = p_i > (1 + pValTol);
        
        if any(any(pvalUnder_i))
            is_pvalUnder_i = true;
        else
            is_pvalUnder_i = false;
        end
        
        if any(any(pvalAbove_i))
            is_pvalAbove_i = true;
        else
            is_pvalAbove_i = false;
        end
end % func CheckPValuesInInterval
%%  CalcEnergyAreas
function [F_ctr, S_ctr] = CalcEnergyAreas(Fdens,p,gam_sumpsqpsq,sumgradpsq,sumpsqpsq,sumpsq,in,tstep,kappa)
% Fdens ... interface energy density

    if length(gam_sumpsqpsq)==1 % gam_sumpsqpsq is only scalar
        gam_sumpsqpsq = sumpsqpsq*gam_sumpsqpsq;
    end

    if in.nOP == 2
        Fdens = Fdens + in.m*(1/4 + gam_sumpsqpsq.*sumpsqpsq);
    else
        Fdens = Fdens + in.m*(1/4 + gam_sumpsqpsq);
    end
    
    
%     F_ctr(1,1) = 2*sum(Fdens(cond_singlePF_IEcalc))*(in.dx)^2; % single grain total IE
    F_ctr(1,1) = 2*sum(Fdens)*(in.dx)^2; % total IE
    
    dx = in.dx;
    % total energy
    % total energy normalized by interface volume
%     F_ctr(1,2) = F_ctr(1,1) / (sum(sumpsqpsq)*(dx)^2 );
    
    % total energy per interface length
    % 2*Fdens_grad./IEfield has _unit_area_ perpendicularly to the intf
    % integration in whole system gives length of the interface
    if strcmp(in.model,'IWvK') || (in.nOP==2 && abs(in.gam0-1.5)<1e-16)
        % field of IE to normalize free energy density Fdens_combined{1}
        % analytic relation holds
        IEfield = sqrt(2/9).*sqrt(kappa*in.m).*ones([in.Ny*in.Nx,1]);
%         plot2D_single_field(IEfield,'IEfield',[in.Ny,in.Nx])
    else
        % polynomial coeff for g(1/gamma)       
        pg = [-0.6130    5.5350  -21.1763   44.7636  -57.2468   45.7228  -22.8980 7.1725 -1.6716    0.7991];
        % field of IE to normalize free energy density Fdens_combined{1}
        if in.nOP == 2
            IEfield = polyval(pg,1./gam_sumpsqpsq).*sqrt(kappa*in.m).*ones([in.Ny*in.Nx,1]);
        else
            IEfield = polyval(pg,sumpsqpsq./gam_sumpsqpsq).*sqrt(kappa*in.m).*ones([in.Ny*in.Nx,1]);
        end
    end    
    
%     F_ctr(1,3) = F_ctr(1,1)/( sum(2*Fdens(cond_singlePF_IEcalc)./IEfield(cond_singlePF_IEcalc))*dx^2 ); % single grain mean IE
    F_ctr(1,2) = F_ctr(1,1)/( sum(2*Fdens./IEfield)*dx^2 );
    
    % mean of interpolation function gives fraction of area of the phase field in the domain
    S_ctr = cellfun(@(x) mean(x.^2./sumpsq),p)';

end % func CalcEnergyAreasArclength

%% plot2D_from_lin_2xn
function plot2D_from_lin_2xn(fignum, cell_with_variables,cell_with_names,size2D,n,nan_alpha)
    assert(length(cell_with_variables)==length(cell_with_names),'Error in plot2D_from_lin_2xn: length(cell_with_variables) ~= length(cell_with_names)')

    figure(fignum)
        for g = 1:(2*n)
            subplot(2,n,g)
                if nan_alpha
                    specs = {'AlphaData',~isnan(reshape(cell_with_variables{g},size2D))};
                else
                    specs = {};
                end
                imagesc(reshape(cell_with_variables{g},size2D),specs{:})
                colorbar
                set(gca,'YDir','normal')
                daspect([1 1 1])
                title(cell_with_names{g})
                
        end
end % func plot2D_from_lin_2xn

%% plot2D_from_lin_2x3
function plot2D_from_lin_2x3(size2D,cell_with_variables,cell_with_names)
    assert(length(cell_with_variables)==length(cell_with_names),'Error in plot2D_from_lin_2x3: length(cell_with_variables) ~= length(cell_with_names)')
    figure(1)
        for g = 1:length(cell_with_variables)
            subplot(2,3,g)
                imagesc(reshape(cell_with_variables{g},size2D))
                colorbar
                set(gca,'YDir','normal')
                daspect([1 1 1])
                title(cell_with_names{g})
        end
end % func plot2D_from_lin_2x3
%% plot2D_single_field
function plot2D_single_field(field_to_plot,name_of_field,size2D)    
    figure(77)
    imagesc(reshape(field_to_plot,size2D))
    colorbar
    set(gca,'YDir','normal')
    daspect([1 1 1])
    title(name_of_field)
end
%% plot2D_from_lin_scalar
function plot2D_from_lin_scalar(size2D,iplot,cell_with_variables,cell_with_names)
    figure(88)
        for g = 1:3
            subplot(1,3,g)
                if iscell(cell_with_variables{g})
                    imagesc(reshape(cell_with_variables{g}{iplot},size2D))
                else
                    imagesc(reshape(cell_with_variables{g},size2D))
                end
                colorbar
                set(gca,'YDir','normal')
                daspect([1 1 1])
                title(cell_with_names{g})
        end
end % func plot2D_from_lin_scalar
%% plot2D_from_lin_dotpr
function plot2D_from_lin_dotpr(size2D,field_in_vector1,title_str1,field_in_vector2,title_str2,DF,title_DF)
    field2D_1x = reshape(field_in_vector1{1},size2D);
    field2D_1y = reshape(field_in_vector1{2},size2D);
    field2D_2x = reshape(field_in_vector2{1},size2D);
    field2D_2y = reshape(field_in_vector2{2},size2D);
    figure(88)
        subplot(121)
            quiver(field2D_1x,field2D_1y,1)
            hold on
            quiver(field2D_2x,field2D_2y,1)
            hold off
            colorbar
            set(gca,'YDir','normal')
            daspect([1 1 1])
            legend(title_str1,title_str2,'location','best')
        subplot(122)
            imagesc(reshape(DF,size2D))
            colorbar
            set(gca,'YDir','normal')
            daspect([1 1 1])
            title(title_DF)
end % func plot2D_from_lin_dotpr
%% plot_basic_info
function plot_basic_info(sumpsq,p,ori_ij,kappa,nnl_ij,size2D,tstep,c)
    RSsumpsq = reshape(sumpsq,size2D);
    RSp_1 = reshape(p{1},size2D);
    RSp_2 = reshape(p{2},size2D);
    RSori = reshape(ori_ij,size2D);
    RSkpp = reshape(kappa,size2D);
%     RSlog10gpl = reshape(log10(nnl_ij),size2D);
    RSc = reshape(c,size2D);

    figure(1)
    set(gcf,'pos',[30,195,1200,555])
    subplot(231)
        imagesc(RSsumpsq),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
%         contour(RSsumpsq,2,'LineColor','k'),set(gca,'YDir','normal'),daspect([ 1 1 1])
        title(['tstep = ' num2str(tstep) ', sum p squared'])
    subplot(232)
        imagesc(RSp_1),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
        title('\xi_1')
    subplot(233)
        imagesc(RSp_2),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
        title('\xi_2')
    subplot(234)
        imagesc(RSori,'AlphaData',~isnan(RSori)),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
        caxis([-pi, pi])
        title('angle of interface normal')
    subplot(235)
        imagesc(RSkpp),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
        title('kappa')   
%         if exist('Lanis')
%             imagesc(reshape(Lanis,size2D)),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
%             title('Laniso')  
%         else                 
%         end
    subplot(236)
        imagesc(RSc),title('c'),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
%         imagesc(RSlog10gpl),title('log_{10}|\nabla\xi_1-\nabla\xi_2|'),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
end % func plot_basic_info

%% DFtoMatrix_i
function [DF,DFname] = DFtoMatrix_i(dt,L,size2D,is_anisotropic_maincalc,df0dp_i,DFlapp_i,DFdivG_i,DFdivK_i,DFvG_i,DFvK_i,DFgK_i)
    DF(:,:,1) = reshape(-dt*L.*df0dp_i,size2D);
    DF(:,:,2) = reshape(dt*L.*DFlapp_i,size2D);
    DFname(1:2) = {'-dt*L*df0dp','dt*L*DFlapp'};
%     if ~is_isotropic && soaIE~=0
    if is_anisotropic_maincalc
        DF(:,:,3) = reshape(dt*L.*DFdivG_i,size2D);
        DF(:,:,4) = reshape(dt*L.*DFdivK_i,size2D);
        DF(:,:,5) = reshape(dt*L.*DFvG_i,size2D);
        DF(:,:,6) = reshape(dt*L.*DFvK_i,size2D);
        DF(:,:,7) = reshape(dt*L.*DFgK_i,size2D);
        DFname(3:7) = {'dt*L*DFdivG','dt*L*DFdivK','dt*L*DFvG','dt*L*DFvK','dt*L*DFgK'};
    end

end % func DFtoMatrix_i

%% plot_RHS_overall
function plot_RHS_overall(in,i,nn_ij,p_i,tstep,DF,DFname,dt,L,df0dp,df0dGp)
    Nx = in.Nx;
    Ny = in.Ny;
    size2D = [Ny,Nx];
    x = 1:Ny;
    y = 1:Nx;
    
    f9 = figure(9);
    f9.Name = 'Homogeneous and total driving force';
        subplot(221)
           quiver(reshape(nn_ij{1},size2D),reshape(nn_ij{2},size2D),1,'Color','k')
           hold on,
           imagesc(x,y,log10(reshape(p_i,size2D)),'AlphaData',0.7*ones(size2D)),colorbar,daspect([1 1 1]),
           hold off
           set(gca,'YDir','normal','xlim',[-2,Ny+2],'ylim',[-2,Nx+2])
           title(['log10(\xi_2) and interf. normal, tstep = ' num2str(tstep) ])
        subplot(222)
            imagesc(DF(:,:,1)),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
            title([DFname{1} ', PF ' num2str(i)])
%             imagesc(DF(:,:,1),'AlphaData',filterCond2D),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])

%             plotlogpm = log10(abs(DF(:,:,1)));
%             plotlogpm(DF(:,:,1)<0) = -plotlogpm(DF(:,:,1)<0);
%             imagesc(plotlogpm,'AlphaData',ones(size2D)),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
%                 title(['log(abs()) ' DFname{1} ', PF ' num2str(i)])
        subplot(223)
                imagesc(reshape(-dt*L*(df0dp{1} - df0dGp{1}),size2D)),colorbar,set(gca,'YDir','normal'),daspect([1 1 1]), title('-dt.L.RHS{1}') 
%                 imagesc(-dt*L*RHS{1}-dt*L*m*sumpsqpsq.*divdgmmdGp{1}),shading flat,colorbar,set(gca,'YDir','normal'),daspect([1 1 1]), title('-dt.L.RHS{1}-dt*L*m*sumpsqpsq.*divdgmmdGp{1}') 
           subplot(224)
                imagesc(reshape(-dt*L*(df0dp{2} - df0dGp{2}),size2D)),colorbar,set(gca,'YDir','normal'),daspect([1 1 1]), title('-dt.L.RHS{2}') 
end % function plot_RHS_overall

%% plot_anisotropic_DF
function plot_anisotropic_DF(DF,DFname,filterCond2D,pind_to_plot)
    for g = 1:3
        f55 = figure(55);
            subplot(1,3,g)
%                 if (g+1)==3
%                     quiverscale = 5;
%                     hold on, quiver(reshape(dgmmdGp{1}{1},size2D),reshape(dgmmdGp{1}{2},size2D),quiverscale),set(gca,'YDir','normal'),daspect([1 1 1]),hold off
%                 end
%                 imagesc(DF(:,:,g+1),'AlphaData',ones(size2D)),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
            imagesc(DF(:,:,g+1),'AlphaData',filterCond2D),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
            title([DFname{g+1} ', PF ' num2str(pind_to_plot)])
        f56 = figure(56);
            subplot(1,3,g)
    % %                 imagesc(DF(:,:,g+4),'AlphaData',ones(size2D)),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
                imagesc(DF(:,:,g+4),'AlphaData',filterCond2D),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
                title([DFname{g+4} ', PF ' num2str(pind_to_plot)])
    end % for g   
            
    f55.Name = 'Scalar-valued anisotropic DF terms';
    f56.Name = 'Dot-product anisotropic DF terms';
            % for plotting DF in log10
            %             plotlogpm = log10(abs(DF(:,:,2)));
            %             plotlogpm(DF(:,:,2)<0) = -plotlogpm(DF(:,:,2)<0);
            %             imagesc(plotlogpm,'AlphaData',ones(size2D)),colorbar,set(gca,'YDir','normal'),daspect([1 1 1])
            %             title(['log(abs()) ' DFname{2} ', PF ' num2str(i)])
end % func plot_anisotropic_DF
%% GetDFInOutliers
function DFsomewhere = GetDFInOutliers(cond,UAswitch,is_anisotropic_maincalc,dt,L,df0dp_i,DFlapp_i,DFdivG_i,DFvG_i,DFdivK_i,DFvK_i,DFgK_i)
    switch UAswitch
        case 'UNDER'
            disp(['Nr of points with value ' UAswitch ': ' num2str(sum(cond))])
        case 'ABOVE'
            disp(['Nr of points with value ' UAswitch ': ' num2str(sum(cond))])
    end
    if all(size(L) == [1,1])
        LL = L;
    else
        LL = L(cond);
    end
        
    DFsomewhere(1,:) = -dt*LL.*df0dp_i(cond);
    DFsomewhere(2,:) = dt*LL.*DFlapp_i(cond) ;
    if is_anisotropic_maincalc
        DFsomewhere(3,:) = dt*LL.*DFdivG_i(cond);
        DFsomewhere(4,:) = dt*LL.*DFdivK_i(cond);
        DFsomewhere(5,:) = dt*LL.*DFvG_i(cond);
        DFsomewhere(6,:) = dt*LL.*DFvK_i(cond);
        DFsomewhere(7,:) = dt*LL.*DFgK_i(cond);
    end
end

%% GetXYCoordOfOutliers
function [xcoor,ycoor] = GetXYCoordOfOutliers(cond,Ny)
    ind = find(cond);
    xcoor = ceil(ind/Ny);
    ycoor = indU-(xcoor-1)*Ny;
end % func GetXYCoordOfOutliers

%% AssignInterpolationFunction
function [interpf,dinterpf] = AssignInterpolationFunction(code)
    if strcmp(code,'p1')
        interpf = @(x) x;
        dinterpf = @(x) ones(size(x));
    elseif strcmp(code,'p3')
        interpf = @(x) x.^2.*(3-2*x);
        dinterpf = @(x) 6*(x-x.^2);
    elseif strcmp(code,'p5')
        interpf = @(x) x.^3.*(6*x.^2-15*x+10);
        dinterpf = @(x) 30*x.^2.*(1-x).^2;
    end % if
end % func

function test_AssignInterpolationFunction(bool)
    if bool
        interpfun_code = {'p1','p3','p5'};
        for k = 1: 3
            [interpf,dinterpf] = AssignInterpolationFunction(interpfun_code{k});
            x = linspace(0,1,100);
            figure(1)
                subplot(1,3,k)
                plot(x,interpf(x),'.',x,dinterpf(x),'.'),grid on,title(interpfun_code{k})
                text(0.5,0.5 ,['\int_0^1 p''(\xi) = ' num2str(trapz(interpf(x))*.01)])
        end
        return
    end % if bool
end % function test