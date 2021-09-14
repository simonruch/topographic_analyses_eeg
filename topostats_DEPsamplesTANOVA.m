function [stat] = topostats_DEPsamplesTANOVA(cfg,varargin)
% topostats_DEPsamplesTANOVA performs topographic analysis of variance for
% dependent samples: different conditions collected in the same
% observational units/subjects
%   =>  tests if different conditions yield systematically different
%       voltage topographies on the scalp across observations.
%       observations: single trials or trial averages for different
%       observational units/subjects; assumes observations are dependent
%       (within-subjects design)
%
% Uses permutation statistics to assess significance of topographic
% dissimilarity
%
% TANOVA:
% 1)compute magnitude of topographic dissimilarity (DISSIM) at each single
%    time-point by calculating the global field power (gfp) of the
%    difference topography between the condition averages of normalized topographies 
%    (multiple conditions: sum of difference between condition averages and
%    global average)
%
% 2)estimate significance (p) of topographic differences at each single time-point
%    using permutation statistics:
%    a) compute random data samples by randomly
%       shuffling condition labels of distinct observatons
%       (conditions are shuffled within subjects)
%    b) compute pseudo-DISSIM for each random sample
%    => The p-value is the percentage of random samples with a larger pseudo-DISSIM at the
%    given time-point; this is not corrected for multiple comparisons
%    across time!
%
% 3)estimate OMNIBUS significance (omnibus p): is there any dissimilarity in
%    the data? Compute the percentage of random samples in which the total 
%    number of time-points that are significant at cfg.alpha is larger
%    than the the sum of significant time-points in the actual data
%
% 4)estimate significance (p) of CLUSTERS of consecutive time-points with
%    significant DISSIMSs using permutation statistics;
%    a) find clusters of consecutive time-points with a p-value <
%    cfg.clusteralpha
%    b) compute the percentage of random samples in which the
%    largest cluster is larger (longer duration / higher number of consecutive
%    time-points that ar sig. at cfg.clusteralpha) than the cluster at
%    hand. This percentage is the p-value of the cluster.
%
% Use as
%   [stat] = topostats_INDEPsamplesTANOVA(cfg, data1, data2, ...)
%
% ... where the INPUT DATA (data 1, ...) is ...
%   EITHER  the result from FT_TIMELOCKANALYSIS or FT_TIMELOCKGRANDAVERAGE.
%           => i.e. fieldtrip data struct with fields .avg or .trial
%   OR      raw data matrices of the dimension n_channels*n_samples
%           => each matrix contains the data of one observational unit (subject/trial)
%           => can be used if data was not processed with fieldtrip toolbox
%
% ... where the CONFIGURATION OPTIONS (cfg) that can/must be specified are:
%   MANDATORY OPTIONS:
%   cfg.ivar              = independent variable (conditions) in cfg.design
%
%   cfg.uvar              = unit-variable (subject-identifier)in cfg.design
%
%   cfg.design            = 2*n matrix;
%                           one row lists the conditions (ivar),
%                           the other row lists the observational units (uvar)
%                           n: number of of observations/data-sets
%                 example with two conditions (1-2) and 5 subjects (1-5)
%                 [1 1 1 1 1 2 2 2 2 2;   => conditions 1-2
%                 [1 2 3 4 5 1 2 3 4 5]   => observational units / subjects 1-5
%
%   NONMANDATORY OPTIONS:
%   cfg.numrandomization  = number of randomizations; default:  5000
%                           if set to 0, computes the the gfp of the map
%                           difference (cfg.stat) without performing
%                           significance tests
%
%   cfg.alpha             = alpha for omnibus test and significance of
%                           clusters; default: 0.05
%
%   cfg.clusteralpha      = value for thresholding for cluster-formation;
%                           default: 0.05
%
%   cfg.normalize         = 'no' or 'yes': should scalp maps be normalized
%                            by global field power at each time-point before
%                           TANOVA is performed? (default: 'yes');
%
%   cfg.parameter         = string (default = 'trial' or 'avg' or 'raw')
%                           use 'raw' if data matrices are provided
%
%   cfg.randomseed        = 'shuffle' (default: uses current time as seed);
%                           'default' (resets seed restart of matlab)
%                           [numeric] => specific integer between 1 and 2^32
%                           => for reproducibility
%
% OUTPUT:
%   stat.stat             = gfp (global field power) of difference between
%                           topographies
%
%   stat.time             = time-vector
%
%   stat.prob             = permutation p-value at each time point
%                           (i.e. percentage of random samples with larger
%                            difference map at this time-point)
%                           => NOT CORRECTED FOR MULTIPLE COMPARISONS
%
%   stat.omnibus_prob     = p-value that there is any difference between
%                           scalp maps somewhere in the data
%                           (i.e. percentage of random samples with higher
%                           sum of single significant time-points)
%
%   stat.clusters         = array with structure for each cluster in the
%                           actual data; fields: 'prob' (p-value of
%                           cluster), 'clusterstat' (size of cluster, i.e.
%                           number of samples that are significant at
%                           cfg.clusteralpha), 'cluster_idx' (indices of
%                           samples belonging to cluster)
%
%   stat.clusterslabelmat = samples are labelled with index of respective
%                           (numbers correspond to nth cluster in
%                           stat.clusters)
%
%   stat.mask             = logic mask indicating significant clusters
%
%   stat.cfg              = configuration options that were used for
%                           running the TANOVA
%
%   stat.cfg              = configuration options
%
%   stat.permutation.designs = designs for drawing random samples
%
%   stat.permutation.stats = gfp of difference maps for all random samples
%
%   stat.permutation.probs = pseudo-probability of permutation samples
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
fprintf(['\nPERFORMING TANOVA\n' ] );

%% BASIC SANITY CHECK
% Check if ANY CONFIGURATION options were defined
if isempty(cfg)
    error('ERROR: no configuration options were provided (cfg is empty);');
end

% Check if ANY DATA was provided
if isempty(varargin)
    error('ERROR: no data was provided;');
end

% Check if ANY DEISGN was specified
if ~isfield(cfg,'design')
    error('ERROR: No design was specified (cfg.design)');
elseif isempty(cfg.design)
    error('ERROR: Empty design specification (cfg.design)');
end

% Check if ivar and uvar are defined in cfg
if ~isfield(cfg,'ivar'); error('ERROR: missing [cfg.ivar]'); end
if ~isfield(cfg,'uvar'); error('ERROR: missing [cfg.uvar]'); end

% Check if ivar and uvar are valid
if ~isnumeric(cfg.ivar); error('ERROR: invalid [cfg.ivar]'); end
if ~isnumeric(cfg.uvar); error('ERROR: invalid [cfg.uvar]'); end
if  ismember(cfg.ivar,cfg.uvar); error('ERROR: [cfg.ivar] and [cfg.uvar] must not be identical'); end



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


%% CHECK SANITY OF DESIGN
% => are dimensions correct?
if numel(size(cfg.design)) ~= 2 |  size(cfg.design,1) ~= 2
    error(['ERROR: design must have 2 diemensions: [2, n_observations]; Current size is: [' num2str(size(cfg.design)) ']']);
end

% => does design match with the input data?
if size(cfg.design,2) ~= numel(data)
    error(['ERROR: mismatch between number of observations in design-specification cfg.design [N=' num2str(size(cfg.design,1)) ...
        '] and number of avaiable observations/datasets [N=' num2str(numel(data)) '].']);
end


%% SUBSET DESIGN AND DATA:
if any(cfg.design(:) == 0)
    warning('Not all observations are used in the design!!');
    data        = data(cfg.design(1,:) ~= 0 & cfg.design(2,:) ~= 0);
    cfg.design  = cfg.design(:,cfg.design(1,:) ~= 0 & cfg.design(2,:) ~= 0 );
end

%% FIX design: ivar = 1; uvar = 2;
cfg.design = cfg.design([cfg.ivar, cfg.uvar],:);
cfg.ivar = 1;
cfg.uvar = 2;

%% CHECK IF DESIGN IS BALANCED
CONDS   = unique(cfg.design(1,:));
UNITS   = unique(cfg.design(2,:));

if numel(UNITS)*numel(CONDS) ~= size(cfg.design,2)
    warning('UNBALANCED DESIGN: missing conditions for some observatonal units! THIS IS IGNORED!');
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
stat.label      = {'topostats_DEPsamplesTANOVA'};
%stat.input      = varargin;
%stat.data       = data;
stat.time       = cfg.time;
stat.stat       = gfp_of_difference_map(cfg.design,DATA);


%% PERFORM STATISTICAL ANALYSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg.numrandomization > 0
    
    %% DRAW RANDOM SAMPLES => permute designs
    disp('Generating random samples via permutation of conditions within subjects!');
    stat.permutation.designs =  permute_designs(cfg.design,cfg.numrandomization);
    
    
    
    %% GET RANDOM DISTRIBUTION => COMPUTE GFP FOR EACH RANDOM SAMPLE
    fprintf('\nComputing statistics for each random sample!\n')
    stat.permutation.stats = [];
    for di = 1:numel(stat.permutation.designs)
        stat.permutation.stats(di,:) = gfp_of_difference_map(stat.permutation.designs{di},DATA);
        
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
    stat.permutation.designs    = [];
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


fprintf(['\nTANOVA completed with [' num2str(cfg.numrandomization) '] randomizations.\n' ] );

end


%%
function [gfp] = gfp_of_difference_map(DESIGN,DATA)
%% COMPUTE GLOBAL FIELD POWER OF DIFFERENCE BETWEEN MAPS
% multiple conditions are allowed: sum of gfp of difference between each
% condition and global mean
%
% Input:
% - DESIGN: 2D matrix: 2*n_observations
%           row 1: condition labels
%           row 2: label of observational units
%                  EXAMPLE:
%                 [1 1 1 1 1 2 2 2 2 2;    ... => conditions 1-2
%                 [1 2 3 4 5 1 2 3 4 5]   ... => observational units / subjects 1-5
% - DATA: 3D matrix: n_observations*n_channels*n_samples

n_channels = size(DATA,2);
n_samples  = size(DATA,3);

% Get unique conditons & observational units
CONDS   = unique(DESIGN(1,:)); %=> list of conditions
UNITS   = unique(DESIGN(2,:)); %=> list of units/subjects

% COMPUTE AVERAGE MAPS FOR EACH CONDITION (unit is irrelevant here)
CONDITION_AVERAGES = [];
for ci = 1:numel(CONDS)
    CONDITION_AVERAGES(ci,1:n_channels,1:n_samples) = mean(DATA(DESIGN(1,:)==CONDS(ci),1:n_channels,1:n_samples),1);
end

% COMPUTE GLOBAL AVERAGE ACROSS ALL CONDITIONS
GLOBAL_AVERAGE = [];
GLOBAL_AVERAGE(1,1:n_channels,1:n_samples) = mean(DATA,1);

% ALTERNATIVE APPROACH
% => first average across conditions for each unit/subject
% => then average across all units/subjects
% for ui = 1:numel(UNITS)
%     GLOBAL_AVERAGE(ui,1:n_channels,1:n_samples) = mean(DATA(DESIGN(2,:)==UNITS(ui),1:n_channels,1:n_samples),1);
% end
% GLOBAL_AVERAGE = mean(GLOBAL_AVERAGE,1);
% => Problem: identical to simply taking 'GLOBAL_AVERAGE = mean(GLOBAL_AVERAGE,1);'
% if design is balanced, but weights observational units with missing
% conditions disproportionally high for unbalanced designs;




% COMPUTE GFP OF DIFFERENCE MAP
% sum of gfp of difference between each condition and global-average
gfp(1,1:n_samples) = squeeze(sum(std(CONDITION_AVERAGES - GLOBAL_AVERAGE,1,2),1));

end

function [design_permutations] = permute_designs(design,n_permutations)
%% PERMUTE DESIGNS: RESCPECT WITHIN-SUBJECT DESIGN
% => permute conditions within obervational units
% => if one condition is missing => permute missing value over conditions
%    (in a subject with only 2 out of 3 conditions, 2 condition-labels are
%    randomly drawn from the list of all conditions for each permutation


UNITS =  unique(design(2,:)); %=> list of units/subjects
CONDITIONS = unique(design(1,:));

design_permutations = {};

% create n_permutations
for pi = 1:n_permutations
    tmp_design = design;
    
    % loop through observational units/subjects
    for ui=1:numel(UNITS)
        % randomly reassign conditions of current observational unit
        tmp_conds = randsample(CONDITIONS,sum((design(2,:)== UNITS(ui))));
        tmp_design(1,(design(2,:)== UNITS(ui))) = tmp_conds;
    end
    design_permutations{pi} = tmp_design;
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

