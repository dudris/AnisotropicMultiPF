%% copyright notice
% This file is part of a dataset <Minar, Martin (2022), “Benchmarking of different strategies to include anisotropy in a curvature-driven multi-phase-field model”, Mendeley Data, V2, doi: 10.17632/5wrv3ky9pp.2> 
% coupled to publication of the same name by Minar, Moelans submitted to Physical Review Materials in January 2022.
% Distributed under GPLv3 license.
% This project has received funding from the European Research Council (ERC) under the European Union's 
% Horizon 2020 research and innovation programme (grant agreement n° 714754).
% 
%% GetAllowedAngles
% [th_allwd, ba_allwd] = GetAllowedAngles(border_angles,numpts,plotting)
% uses linspace function to obtain allowed angles between the border points.
% INPUT
%   - border_angles ... as many rows as FORBIDDEN intervals, start and
%   endpoint in he 1st and 2nd column
%   - numpts ... number of output points in every allowed interval
%   - plotting ... bool to turn on or off the visualization of output
% OUTPUT
%   - th_allwd ... vector of allowed angles
%   - ba_allwd ... rearranged border_angles so that in every row are start and end points of ALLOWED interval

function [th_allwd, ba_allwd] = GetAllowedAngles(border_angles,numpts,plotting)
    
    if all(size(border_angles)==[1,1]) % angles are allowed
        th_allwd = linspace(-pi,pi,numpts+1)';
        th_allwd(1)=[];
        ba_allwd= 0;
    
    else
        nfold = size(border_angles,1);

        % __ makes sure  1st column is always smaller than 2nd
        ba_forb(:,1) = min(border_angles,[],2);
        ba_forb(:,2) = max(border_angles,[],2);

        % __ rearrange to identify interval of ALLOWED angles
            % __ the below fails when there should be jump in ALLOWED angles interval
            % __ then the rearranged still marks FORBIDDEN intervals
        ba_allwd = reshape(sort(ba_forb(:)),[2,nfold])';
        negcheck = all(sort(ba_allwd(:,1))==sort(ba_forb(:,1))); % true when failed

        if negcheck 
            % __ assume that allowed interval must have a jump
            clear ba_allwd
            ba_sorted = sort(ba_forb(:));
            % __ largest and smallest value form the interval with jump
            ba_allwd(1,[1,2]) = [max(ba_sorted), (min(ba_sorted)+2*pi) ];
            % __ the rest to be processed as before
            ba_temp = ba_sorted(2:(end-1));
            ba_allwd(2:nfold,:) = reshape(ba_temp,[2,(nfold-1)])';
            negcheck = all(sort(ba_allwd(:,1))==sort(ba_forb(:,1)));
            assert(~negcheck,'GetAllowedAngles in calc_regularized_Wulff_normal_ang failed')
        end

        % __ assumes 1st column is always smaller than 2nd
        int_width = diff(ba_allwd,[],2);
        jump_segment_ind = find((int_width-max(int_width))>pi);
        th_allwd =  nan(nfold*numpts,1);

        indkk = 1:numpts;
        for kk = 1:nfold
            indkk = (1:numpts) + numpts*(kk-1);

            if isempty(jump_segment_ind) || (kk~=jump_segment_ind)
                th_allwd(indkk,1) = linspace(ba_allwd(kk,1),ba_allwd(kk,2),numpts);
            else
                % __ segment with jump 
                negind = 1:floor(numpts/2);
                posind = (floor(numpts/2)+1):numpts;
                th_allwd(indkk(negind),1) = linspace(ba_allwd(kk,1),-pi,numel(negind));
                th_allwd(indkk(posind),1) = linspace(ba_allwd(kk,2),pi,numel(posind));
            end% if
        end% for

        th_allwd = rotate_to_first_inetrval(th_allwd,0);
        th_allwd = unique(th_allwd); % to throuw away duplicaets and sort
    end % if weak anisotropy
   
    if plotting
        figure(66)
        x = 1:numel(th_allwd); plot(x,th_allwd,'o',x(isnan(th_allwd)),th_allwd(isnan(th_allwd)),'rx'), 
        polarplot(th_allwd,ones(size(th_allwd)),'.')
        title('allowed angles intervals')
    end
% sum(isnan(th_mod))

end 