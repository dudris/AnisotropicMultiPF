%% [t_ctr, t_ctr_p]= get_sim_timeline(resultsentry)
% t_ctr ... simulation time in equidistant timesteps with spacing resultsentry.ctrcnt AND in the timesteps where the simulated phase fields were saved
% t_ctr_p... simulation time in timesteps where the simulated phase fields were saved
function [t_ctr, t_ctr_p]= get_sim_timeline(resultsentry)

    re = resultsentry;

    t_ctr_p = re.Ndt;
    if re.outputAtAllCtr{1} 
        t_ctr_p = unique([re.in.outputAtAllCtr{2}:re.in.outputAtAllCtr{2}:re.Ndt , t_ctr_p]);
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