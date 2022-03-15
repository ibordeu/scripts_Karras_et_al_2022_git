%%-----------------------------------------------------------------------%%
% This script loads the clone sizes and output from script_biexponential_fit
% (summarised in Table 3), and uses those parameters to run stochastic 
% simulations of the two-compartment SP model. 
% This script produces plots comparing the empirical cumulative distribution 
% of clone sizes with the obtained from the numerical results (Extended data
% Figure 5f-j). 
%
%%-----------------------------------------------------------------------%%
close all; clear all; clc;

% LOAD DATA 
% read clone sizes a cell array: 5 samples, 2 channels each
samples = 1:5; channels =1:2; channel_labels = {'RFP','YFP'};
clone_sizes = read_spreadsheet('clone_sizes.xlsx',samples,channels);
% read "result_table" output from script_biexponential_fit
load('result_table.mat');
% ignore_singlets : flag to ignore single cell clones, set to 1 to ignore, 
% 0 to keep. Singles are ignored in only to consider only clones that 
% originated from cells that were proliferative at the moment of induction. 
ignore_singlets = 1;
%
N_max = Inf;
n_realis = 10000; % to speed up, lower this number to 1000 or 100;

% We run through each dataset and run 100 repeats of the
pvalues = [];
for nsample = samples
    figure;
    all_sim_clone_sizes = [];
    for nchannel = channels
        clones = clone_sizes{nsample,nchannel};
        if ignore_singlets == 1
            clones(clones == 1) = [];
        end
        nclones = length(clones);
        % extract parameters from result_table
        params = result_table(2*nsample-1 + nchannel-1,:);
        t0 = params.t0;
        r = params.r;
        sigma = params.sigma;
        Delta_p =  params.Delta_p;
        fs =  params.fs;
        % Run numerical simulations and plot the result for each channel in
        % a separate subplot
        subplot(1,length(channels),nchannel)
        binedges = 0.5:1:10000.5; % define bins
        cumNsims = zeros(n_realis,length(binedges)-1); % preallocate array
        % run n_realis realizations
        for n = 1:n_realis
            % each realisation simulates nclones clones
            clone_info = function_SP_model(sigma,r,Delta_p,nclones,fs,t0,N_max);
            % clone size = number of S cells + number of P cells in clone
            sim_clones_size = clone_info.numS + clone_info.numP;
            % remove singlets if ignore_singlets == 1
            if ignore_singlets == 1
                sim_clones_size(sim_clones_size == 1) = [];
            end
            % construct cdf for every realisation
            [~,cumNsim,bincents] = cumulativeProb(sim_clones_size,binedges);
            cumNsims(n,:) = cumNsim;
            % store the size of all clones (for KS-test later)
            all_sim_clone_sizes = [all_sim_clone_sizes;sim_clones_size]; 
        end
        % plot numerical cdf (average +/- SD) 
        cdfMean = mean(cumNsims,1); cdfSD = std(cumNsims,[],1);
        errorbar(bincents,cdfMean,cdfSD,'Color',[0.9 0.9 0.9],'HandleVisibility','off')
        hold on
        plot(bincents,mean(cumNsims,1),'-','Color',[0.7 0.7 0.7])
        % plot empirical cdf
        [~,cumN,bincents] = cumulativeProb(clones,binedges);
        plot(bincents,cumN,'k','linewidth',1.5); % plot empirical CDF
        hold off
        % set y-axis to log scale. Comment out this line for linear-linear
        set(gca,'YScale','log')
        % add labels and adjust axes for better visualisation
        [minval,mind] = min(cumN(cumN~=0)); 
        axis([0 bincents(mind) 10^-2 1]); % adjust margins
        xlabel('Clone size (cells)'); ylabel('Cumulative probability')
        title(channel_labels{nchannel}); % add title RFP/YFP
        legend({'model','data'}); legend boxoff
        % perform KS-test between empirical and numerical distribution of
        % clones sizes.
        [~,p] = kstest2(all_sim_clone_sizes,clones);
        pvalues(end+1,1) = p;
    end
    suptitle(['Sample ',num2str(nsample)])
    set(gcf,'color','w');
end
% store p-value in the result_table and update it (save).
result_table.pvalues = pvalues;
save('result_table.mat','result_table')

%% Functions
function clone_info = function_SP_model(sigma,r,Delta_p,nclones,prop_s,t_max,N_max)
% Simulations of the two compartment SP model, where a population of self-r
% enewing stem cells (S) that produce self-renewing progenitors (P).
% Supplemental Note Eq. (1)
%
% S : stem cells
% P : Progetnitor cells
% processes:
% sigma   : S -> S+S Prob. r (renewal through symmetric division)
% sigma   : S -> P Prob 1-r (loss through differentiation)
% Delta_p : P -> P+P (renewal through symmetric division)
% with expansion rates:
% Delta_s = sigma(2r-1), for S cells
% Delta_p for P cells
% ----------------------------------------------------------------------- %
% Inputs:
% sigma   : S cycling rate
% r       : S duplication probability
% Delta_p : P expansion rate
% nclones : Number of realizations (each realization simulates a single clone)
% prop_s  : Proportion of Nreps that are initialised with a single S cell
%           i.e., S(t=0)=1, P(t=0)=0, the remaining 1-prop_s realisations are 
%           initialized with a single progenitor i.e., S(t=0)=0, P(t=0)=1.
% t_max   : Maximum simulation time (equated to the chase time)
% N_max   : Upper bound to the total number of particles S+P in the system,
%           can be set to Inf.
%
% A single realisation terminates if t_max or N_max are reached, or if the
% total particle count S+P is zero.
% If S+P is zero, then the realization in re-run.
%
% Outputs:
% clone_sizes : table where each row shows the results for a single 
%               realisation. Columns show total simulation time, initial 
%               number of S (numS_init) and P (numP_init) cells, and final 
%               number of S (numS) and P (numP) cells. The clone size is 
%               then numS+numP.
% -------------------------------------------------------------------------
clone_sizes = zeros(nclones,5); % initialize array to store results
rep = 0;
while rep < nclones % run until Nreps clones are recorded
    % initial condition:
    if rep <= prop_s*nclones
        numS_init = 1; numP_init = 0;
    else
        numS_init = 0; numP_init = 1;
    end
    % initialize variable
    numS = numS_init;  % number of stem cells
    numP = numP_init;  % number of progenitors
    t = 0;             % initial time
    % run realisation until reaching maximum time or number of particles.
    % or the clone dies (numS + numP == 0).
    while t < t_max && (numS + numP) < N_max && (numS + numP) > 0
        % calculate propensity functions:
        ws1 = sigma*r*numS;
        ws2 = sigma*(1-r)*numS;
        wp = Delta_p*numP;
        w = ws1 + ws2 + wp;
        % update time until the next event
        dt =  - log(1-rand()) / w;
        t = t + dt;
        % choose an action to perform
        ran = rand();
        if ran < ws1/w              % S -> S + S
            numS = numS + 1;
        elseif ran < (ws1+ws2)/w    % S -> P
            numS = numS - 1;
            numP = numP + 1;
        elseif ran < (ws1+ws2+wp)/w % P -> P + P
            numP = numP + 1;
        end  
    end
   
    if (numS + numP) > 0 % if the clone is "alive"
        % store clone size, and initial condition information
        rep = rep + 1;
        clone_sizes(rep,:) = [t, numS_init, numP_init, numS, numP];
    end
end
% return table:
clone_info = array2table(clone_sizes,...
    'VariableNames',{'t','numS_init','numP_init','numS','numP'});
end

function [N,cumN,bincents] = cumulativeProb(data,binedges)
% Returns the probabilities N and cumulative probabilities cumN, and
% corresponding bincentres for the Nx1 data array. 
if nargin == 1
binedges = 0.5:1:max(data)*1+0.5; 
end
bincents = (binedges(2:end) + binedges(1:end - 1))/2;
N = histcounts(data,binedges,'Normalization','pdf');
cumN = 1 - cumsum(N)./sum(N);
end

function clone_sizes = read_spreadsheet(file_path,samples,channels) 
% read spreadsheet and turn data into cell array
clone_sizes_xls =  xlsread(file_path); % load numeric values from excel sheet
clone_sizes_xls(:,1) = []; % first row has just clone indices;
clone_sizes = {};
for nsample = samples
    for nchannel = channels
        clone_sizes_sample = clone_sizes_xls(:,2*nsample-1+nchannel-1);
        % remove nans create by loading empty cells.
        clone_sizes_sample(~isfinite(clone_sizes_sample)) = []; 
        clone_sizes{nsample,nchannel} = clone_sizes_sample;
    end
end
end