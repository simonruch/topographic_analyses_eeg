function [cluster_stats] = topostats_ReportStats(cfg,stats)
% topostats_ReportStats(stats) reports the key data of all significant
% clusters for any of the the topostats_* analyses (TANOVA, GFP, TCT)
%
% Use as
%   [cluster_stats] = topostats_ReportStats(cfg, stats) ...
%
% ... where the INPUT DATA (stat) ...
%     => is the result from any of the topostats_...-analyses
%
% ... where the CONFIGURATION OPTIONS (cfg) that can be specified are ...
%   cfg.alpha =     alpha value for significant clusters (default = 0.05) set
%                   to cfg.alpha = 1 to plot all clusters
%   cfg.doprint =   whether or not to print the cluster-data to terminal
%                   (default = true)
%
%
% Output:
% table 'cluster_stats' with the variables:
% cluster_stats.time_min:       onset of cluster
% cluster_stats.time_max:       offset of cluster
% cluster_stats.time_dur:       duration of cluster
% cluster_stats.prob:           permutation-based probability of the cluster
% cluster_stats.omnibus_prob:   overall probability of no difference/effect in the data
%                 (the same for all clusters)



%% BASIC SANITY CHECK
% Check if ANY CONFIGURATION options were defined
if isempty(cfg)
    warning('no configuration options were provided (cfg is empty);');   
end

% Check if ANY DATA was provided
if isempty(stats)
    error('ERROR: no data was provided;');
end


%% GET OPTIONAL CONFIGURATION OPTIONS
if ~isfield(cfg,'alpha');               cfg.alpha               = 0.05; end
if ~isfield(cfg,'doprint');             cfg.doprint             = true; end


%% Extract statistics
cluster_stats = struct();

% chek if there are clusters
if isfield(stats,'clusters')
    % get clusters with alpha < cfg.alpha
    if ~isempty(stats.clusters)
        idx_sigclusters = [];
        idx_sigclusters = find([stats.clusters.prob]<=cfg.alpha);
        
        % extract stats
        if ~isempty(idx_sigclusters)
            for ic = 1:numel(idx_sigclusters)
                tmp_clst = struct();
                tmp_clst.time_min       = min(stats.time(stats.clusters(idx_sigclusters(ic)).cluster_idx));
                tmp_clst.time_max       = max(stats.time(stats.clusters(idx_sigclusters(ic)).cluster_idx));
                tmp_clst.time_dur       = tmp_clst.time_max-tmp_clst.time_min;
                tmp_clst.prob           = stats.clusters(idx_sigclusters(ic)).prob;
                tmp_clst.omnibus_prob   = stats.omnibus_prob;
                
                if ic == 1
                    cluster_stats = tmp_clst;
                else 
                    cluster_stats = [cluster_stats; tmp_clst];
                end
                
            end            
        end                                
    end    
end

% Convert to table
cluster_stats = struct2table(cluster_stats);


% Print table to terminal if cfg.doprint = true
if cfg.doprint 
    disp(cluster_stats);
end


