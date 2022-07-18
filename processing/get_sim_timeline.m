%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Benchmarking of different strategies to include anisotropy in a curvature-driven multi-phase-field model”, Mendeley Data, V2, doi: 10.17632/5wrv3ky9pp.2> 
% coupled to publication of the same name by Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754).
% 
%% [t_ctr, t_ctr_p]= get_sim_timeline(resultsentry)
% retrieve time line from input
% t_ctr ... simulation time in timesteps of checkpoint energy and area calculations AND in the timesteps where the simulated phase fields were saved
% t_ctr_p... simulation time in timesteps where the simulated phase fields were saved
function [t_ctr, t_ctr_p]= get_sim_timeline(resultsentry)

    re = resultsentry;

    t_ctr_p = re.Ndt;
    
    if isfield(re,'outputAtAllCtr')
        if re.outputAtAllCtr{1} 
            t_ctr_p = unique([re.in.outputAtAllCtr{2}:re.in.outputAtAllCtr{2}:re.Ndt , t_ctr_p]);
        end
    end
    
    if isfield(re,'outputAtChckpt')
        if ~isempty(re.outputAtChckpt{1})
            if re.outputAtChckpt{1}
                t_ctr_p = unique([t_ctr_p , re.outputAtChckpt{3}]);
            end
        end    
    end
    
    if isfield(re,'ctrcnt')
        if isempty(re.ctrcnt)
            t_ctr = unique([1 re.ctrplot:re.ctrplot:re.Ndt t_ctr_p]);
        else
            t_ctr = unique([re.tsteptsctr t_ctr_p]);
        end
    else
        t_ctr = unique([1 re.ctrplot:re.ctrplot:re.Ndt t_ctr_p]);
    end
    
    t_ctr_p = double(t_ctr_p')*re.dt;
    t_ctr = double(t_ctr')*re.dt;
end