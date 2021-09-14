function [stat] = topostats_SingleConditionTCT(cfg,varargin)
% topostats_SingleConditionTCT performs a TOPOGRAPHIC CONSISTENCY TEST (TCT)
% over data of one single condition
%   =>  tests if topographies are consistent
%       across observations at each time-point
%       observations: single trials or trial averages for different
%       observational units/subjects
%
% Uses permutation statististics to assess significance of TC
%
% TCT:
% 1)compute topographic consistency TC
%    => the global field power (gfp) of the average of normalized
%    topographies
% 2)estimate significance (p) of TC at each single time-point
%    using permutation statistics:
%    a) compute random data samples by randomly
%       permuting electrodes across observations
%    b) compute pseudo-TC for each random sample
%    => The p-value is the percentage of random samples with a larger pseudo-TCs at the
%    given time-point; this is not corrected for multiple comparisons
%    across time!
% 3)estimate OMNIBUS significance (omnibus p): is there any consistence somewhere in
%    the data? Compute the percentage of random samples in which the total 
%    number of time-points that are significant at cfg.alpha is larger
%    than the the sum of significant time-points in the actual data
% 4)estimate significance (p) of clusters of consecutive time-points with
%    significant TCs using permutation statistics;
%    a) find clusters of consecutive time-points with a p-value <
%    cfg.clusteralpha
%    b) compute the percentage of random samples in which the
%    largest cluster is larger (longer duration / higher number of consecutive
%    time-points that ar sig. at cfg.clusteralpha) than the cluster at
%    hand. This percentage is the p-value of the cluster.
%
% Use as
%   [stat] = topostats_SingleConditionTCT(cfg, data1, data2, ...) ...
%
% ... where the INPUT DATA (data1, ...) is ...
%   EITHER  the result from FT_TIMELOCKANALYSIS or FT_TIMELOCKGRANDAVERAGE.
%           => i.e. fieldtrip data struct with fields .avg or .trial
%   OR      raw data matrices with the dimension [n_channels * n_samples]
%           => each matrix contains the data of one observational unit (subject/trial)
%           => can be used if data was not processed with fieldtrip toolbox
%
% ... where the CONFIGURATION OPTIONS (cfg) that can/must be specified are ...
%
%   MANDATORY OPTIONS:      none
%
%   NONMANDATORY OPTIONS:
%
%   cfg.numrandomization  = number of randomizations; default:  5000
%                           if set to 0, computes the topographic
%                           consistency (cfg.stat) without performing
%                           significance tests
%
%   cfg.alpha             = alpha for omnibus test and significance of
%                           clusters (default: 0.05)
%
%   cfg.clusteralpha      = value for thresholding for cluster-formation
%                           (default: 0.05)
%
%   cfg.normalize         = 'no' or 'yes': should topographies be normalized
%                           by global field power at each time-point for 
%                           each observation before TCT is performed?
%                           (default: 'yes');
%
%   cfg.parameter         = string (default = 'trial' or 'avg' or 'raw')
%                           use 'raw' if data matrices are provided
%
%   cfg.randomseed        = 'shuffle': (default: uses current time as seed);
%                           'default': (resets seed restart of matlab)
%                           [numeric]: specific integer between 1 and 2^32
%                                      => for reproducibility
%
% OUTPUT:
%   stat.stat             = topographic consistency: gfp (global field power) 
%                           of average of (normalized) topography at each
%                           time-point
%
%   stat.time             = time-vector
%
%   stat.prob             = permutation p-value at each time point
%                           (i.e. percentage of random samples with larger
%                            TC at this time-point)
%                           => NOT CORRECTED FOR MULTIPLE COMPARISONS
%
%   stat.omnibus_prob     = p-value that there is any TC somewhere in the data
%                           (i.e. percentage of random samples with higher
%                           sum of single significant time-points)
%
%   stat.clusters         = array with structure for each cluster in the
%                           actual data; fields:
%                           'prob': (p-value of cluster)
%                           'clusterstat': size of cluster, i.e. number of 
%                           samples that are significant at
%                           cfg.clusteralpha
%                           'cluster_idx': indices of samples belonging to cluster
%
%   stat.clusterslabelmat = samples are labelled with index of respective
%                           cluster (numbers correspond to nth cluster in
%                           stat.clusters)
%
%   stat.mask             = logic mask indicating time points belonging to 
%                           significant clusters
%
%   stat.cfg              = configuration options that were used for
%                           running the TCT
%
%   stat.permutation.stats = TC for all random samples
%
%   stat.permutation.probs = pseudo-probability of random samples
%
%
% REFERENCES:
% General description of TCT and TANOVA:
% Koenig, T., & Melie-García, L. (2009). Statistical analysis of multichannel scalp field data. In C. M. Michel, D. Brandeis, J. Wackermann, L. R. R. Gianotti, & T. Koenig (Eds.), Electrical Neuroimaging (pp. 169–190). Cambridge University Press. https://doi.org/10.1017/CBO9780511596889.009
%
% Detailed description of randomisation (with Matlab code):
% Koenig, T., & Melie-García, L. (2010). A method to determine the presence of averaged event-related fields using randomization tests. Brain Topography, 23(3), 233–242. https://doi.org/10.1007/s10548-010-0142-1
%
% Detailed description of TANOVA:
% Koenig, T., Kottlow, M., Stein, M., & Melie-García, L. (2011). Ragu: A free tool for the analysis of EEG and MEG event-related scalp field data using global randomization statistics. Computational Intelligence and Neuroscience, 2011, 1–14. https://doi.org/10.1155/2011/938925


%%
fprintf(['\nPERFORMING TCT\n' ] );

%% BASIC SANITY CHECK
% Check if ANY CONFIGURATION options were defined
if isempty(cfg)
    error('ERROR: no configuration options were provided (cfg is empty);');   
end

% Check if ANY DATA was provided
if isempty(varargin)
    error('ERROR: no data was provided;');
end



%% GET OPTIONAL CONFIGURATION OPTIONS
if ~isfield(cfg,'alpha');               cfg.alpha               = 0.05; end
if ~isfield(cfg,'clusteralpha');        cfg.clusteralpha        = 0.05; end
if ~isfield(cfg,'numrandomization');    cfg.numrandomization    = 5000; end
if ~isfield(cfg,'normalize');           cfg.normalize           = 'yes'; end


if ~isfield(cfg,'randomseed');
    cfg.randomseed          = 'shuffle';
elseif isnumeric(cfg.randomseed)
elseif isstr(cfg.randomseed)
    if ~ismember(cfg.randomseed,{'shuffle','default'})
        error('ERROR: invalid input for cfg.randomseed: [''shuffle'',''default'', or numeric value as seed]');
    end
else
    error('ERROR: invalid input for cfg.randomseed: [''shuffle'',''default'', or numeric value as seed]');
end
rng(cfg.randomseed);





%% IDENTIFY PARAMETER THAT SHOULD BE ANALYZED (default; 'raw','avg', or 'trial')
% If no parameter was specified: determine whether it is 'raw', 'trial', or 'avg'
if ~isfield(cfg,'parameter') || (isfield(cfg,'parameter') & isempty(cfg.parameter))
    % if cfg.parameter is not specified or is empty
    if isnumeric(varargin{1})
        % => if first data set that is entered is a numeric matrix => assume raw data
        cfg.parameter = 'raw';
    elseif isstruct(varargin{1})
        % => if first data set that is entered is a struc() => assume fieldtrip data
        if isfield(varargin{1},'trial')
            % => must have trials ...
            cfg.parameter = 'trial';
        elseif isfield(varargin{1},'avg')
            % => or must have trial averages ...
            cfg.parameter = 'avg';
        else
            error('ERROR: unknown parameter / data-type ');
        end
    else
        error('ERROR: unknown data-type; data must be fieldtrip data structures or numeric matrices');
    end
end

% Check if parameter is available in data
if isfield(cfg,'parameter')
    if ~ismember(cfg.parameter, {'raw'})
        if ~isfield(varargin{1},cfg.parameter)
            error(['ERROR: parameter [' cfg.parameter '] not available in data;'] );
        end
    end
end

%% GET TIME
if ismember(cfg.parameter,{'raw'})
    cfg.time = [1:size(varargin{1},2)];
elseif ismember(cfg.parameter,{'trial'})
    cfg.time = varargin{1}.time{1};
else
    cfg.time = varargin{1}.time;
end


%% GET DATA
data = [];
if ismember(cfg.parameter,{'raw'})
    % if data-type is 'raw':
    %   => list all matrices in cells
    for vi = 1:numel(varargin)
        data(vi) = varargin(vi);
    end
    
elseif ismember(cfg.parameter,{'trial'})
    % if data-type is 'trial':
    %   => append all trials of all data sets
    for vi = 1:numel(varargin)
        data = [data varargin{vi}.trial(:)' ];
    end
    
else
    % if data-type is 'avg' or user-defined
    %   => assume that it is a fieldtrip data-struct
    for vi = 1:numel(varargin)
        data{vi} = varargin{vi}.(cfg.parameter);
    end
end



%% CHECK SANITY OF DATA
for di = 1:numel(data)
    if ~isnumeric(data{di})
        error('ERROR: data not numeric');
    end
end

if numel(size(data{1})) ~= 2
    error('ERROR: invalid dimensionality of data (must be n_channels*n_samples');
end

if size(data{1},1) < 2
    error('ERROR: first dimension of data <2; data must have at least two channels');
end

% => check for unequal sizes of data matrices
%    (unequal number of channels or samples across observations)
for di = 1:numel(data)
    if di > 1
        if ~all(size(data{di}) == size(data{di-1}))
            error('ERROR: incoherent input data (unequal number of channels and/or samples across observations)!');
        end
    end
end

%% NORMALIZE DATA => spatial normalization of fieldmaps
if ismember(cfg.normalize,{'yes'})
    for di = 1:numel(data)
        data{di} = (data{di} - mean(data{di},1))./ std(data{di},1,1);
    end
end




%% GET DATA AS MATRIX
% Put data into matrix n_observations*n_channels*n_samples
DATA = cat(3,data{:});          %=> concatenate 2D data in 3rd dimension
%   yields: n_channels*n_samples*n_observations
DATA = permute(DATA,[3,1,2]);   %=> permute dimensions:
%   yields: n_observations*n_channels*n_samples



%% PREPARE OUPTUT
stat = [];
stat.cfg        = cfg;
stat.label      = {'topostats_SingleConditionTCT'};
%stat.input      = varargin;
%stat.data       = data;
stat.time       = cfg.time;
stat.stat       = topographic_consistency(DATA);


%% PERFORM STATISTICAL ANALYSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg.numrandomization > 0
    
    %% DRAW RANDOM SAMPLES
    %%% & GET RANDOM DISTRIBUTION => COMPUTE GFP FOR EACH RANDOM SAMPLE
    disp('Generating random samples via permutation of electrodes within observations!');
    fprintf('\nComputing statistics for each random sample!\n')
      
    stat.permutation.stats = [];
    for di = 1:cfg.numrandomization
        stat.permutation.stats(di,:) = topographic_consistency(permute_data(DATA));        
        progress_plot(di);
    end
    
    
    
    %% COMPUTE PROBABILITY: percentage of random samples that have same or higher gfp than real gfp
    stat.prob = mean(([stat.stat; stat.permutation.stats] - stat.stat)>=0,1);
    
    
    
    
    %% COMPUTE PROBABILITIES FOR PERMUTATION STAT
    fprintf('\nComputing pseudo-probabilities for each random sample!\n');
    
    stat.permutation.probs = zeros(size(stat.permutation.stats));
    
    for pi = 1:size(stat.permutation.stats,1)
        stat.permutation.probs(pi,:) = ...
            mean(...
            (stat.permutation.stats -       ...
            stat.permutation.stats(pi,:))   ...
            >=0,1);
        
        progress_plot(pi);
        
    end
    
    %% COMPUTE OMNIBUS PROBABILITY
    % percentage of random samples with more significant time-points
    stat.omnibus_prob = ...
        (sum(...
        (sum(stat.permutation.probs <= cfg.alpha,2) ...  sum of significant time-points in each random sample
        - sum(stat.prob <= cfg.alpha)) >=0)+1) ...         sum of significant time-points in real data
        /(size(stat.permutation.stats,1)+1);               % number of iterations
    
    %% COMPUTE RANDOM DISTRIBUTION OF SIGNIFICANT CLUSTERS
    
    % =>    for each permuation, find all clusters of coherent time-points that
    %       are significant at cfg.clusterlpha
    cluster_sizes      = [];
    max_cluster_sizes  = [];
    
    for pi = 1:size(stat.permutation.probs,1)
        cluster_labels  = [];
        cluster_labels  = bwlabel(stat.permutation.probs(pi,:)<=cfg.clusteralpha);
        tmp_clust_sizes = [];
        if max(cluster_labels>0)
            for ci = 1:max(cluster_labels)
                tmp_clust_sizes(ci) = sum(cluster_labels == ci);
            end
        else
            tmp_clust_sizes = [0];
        end
        
        
        cluster_sizes{pi}       = tmp_clust_sizes;
        max_cluster_sizes(pi)   = max(tmp_clust_sizes);
        
    end
    
    stat.clustersizedistribution_maxclust   = sort(max_cluster_sizes,'descend');
    stat.clustersizedistribution_allclust   = sort([cluster_sizes{:}],'descend');
    
    
    %%
    cluster_labels  = [];
    cluster_labels  = bwlabel(stat.prob<=cfg.clusteralpha);
    
    clusters = [];
    
    if max(cluster_labels>0)
        % get all availabel clusters in real data
        for ci = 1:max(cluster_labels)
            tmp = [];
            tmp.clusterstat = sum(cluster_labels == ci);
            tmp.cluster_idx = find(cluster_labels == ci);
            clusters = [clusters; tmp];
        end
        
        % sort clusters by their size (descending)
        [~,idx] = sort([clusters.clusterstat],'descend');
        clusters = clusters(idx);
        
        stat.clusters = clusters;
        
        % generate clusterslabelmat:
        % => labels correspond to index of sorted cluster
        %    (1 = largest cluster, 2 = 2nd largest cluster ...)
        stat.clusterslabelmat = zeros(size(cluster_labels));
        for ci = 1:numel(stat.clusters)
            stat.clusterslabelmat(stat.clusters(ci).cluster_idx) = ci;
        end
        
        for ci = 1:numel(stat.clusters)
            stat.clusters(ci).prob = ...
                (sum(...
                stat.clustersizedistribution_maxclust- stat.clusters(ci).clusterstat >=0)+1) ...
                /(numel(stat.clustersizedistribution_maxclust)+1);
        end
        
        
    else
        stat.clusters = [];
        stat.clusterslabelmat = zeros(cluster_labels);
    end
    
    
    %%
    stat.mask = zeros(size(stat.prob));
    
    if ~isempty(stat.clusters) && any([stat.clusters.prob]<= cfg.alpha)
        sig_clusts = find([stat.clusters.prob]<=cfg.alpha);
        for ci = 1:numel(sig_clusts)
            stat.mask(stat.clusters(sig_clusts(ci)).cluster_idx) = 1;
        end
        
    end
    
else
    stat.permutation.stats      = [];
    stat.prob                   = [];
    stat.permutation.probs      = [];
    stat.omnibus_prob           = [];
    stat.clustersizedistribution_maxclust= [];
    stat.clustersizedistribution_allclust= [];
    stat.clusters               = [];
    stat.clusterslabelmat       = [];
    stat.mask                   = [];
end


fprintf(['\nTCT completed with [' num2str(cfg.numrandomization) '] randomizations.\n' ] );

end


%%
function [tc] = topographic_consistency(DATA)
%% TOPOGRAPHIC CONSISTENCY ACROSS OBSERVATIONS
% INPUT:
% - DATA: 3D matrix: n_observations*n_channels*n_samples

n_channels = size(DATA,2);
n_samples  = size(DATA,3);

% COMPUTE GLOBAL AVERAGE ACROSS ALL CONDITIONS
GLOBAL_AVERAGE = [];
GLOBAL_AVERAGE(1,1:n_channels,1:n_samples) = mean(DATA,1);

% COMPUTE GFP OF GLOBAL AVERAGE
tc(1,1:n_samples) = squeeze(std(GLOBAL_AVERAGE,1,2));

end

function [DATA_PERMUTED] = permute_data(DATA)
%% PERMUTE ELECTRDOES:
% => permute electrode labels


UNITS =  size(DATA,1);      %=> number of units/subjects
ELECTRODES = size(DATA,2);  % => number of electrodes

DATA_PERMUTED = DATA;

% loop through observational units/subjects
for ui=1:numel(UNITS)
    % randomly reassign shuffle data of electrodes for each unit
    DATA_PERMUTED(ui,:,:) = DATA(ui,randsample(ELECTRODES,ELECTRODES),:);
end

end

function progress_plot(idx)
if idx == 1
    fprintf('\n .');
end

if mod(idx,10)==0
    fprintf('.');
end

if mod(idx,1000) == 0
    fprintf('\n .');
end

end

