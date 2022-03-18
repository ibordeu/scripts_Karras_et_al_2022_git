%%-----------------------------------------------------------------------%%
%%- PARAMETER FITS FOR THE TWO-COMPARTMENT (SP) MODEL -------------------%%
%%-----------------------------------------------------------------------%%
% This script loads the clone sizes, construct the cumulative distribution
% of clone sizes and estimates the parameters of the bi-exponential
% distribution, Supplemental Note Eq. (10), following the strategy
% described in the "Model fits" section.
%
% In this script we generate:
% Figure 2F, Supplementary Figures 5b-5f, and Tables 3-4.
%%-----------------------------------------------------------------------%%
close all; clear all; clc;
%-DATA AND PARAMETERS-----------------------------------------------------%
% Read spreadsheet as a cell array: 5 samples, 2 channels each
samples = 1:5; channels =1:2; % 1-RFP, 2-YFP
clone_sizes = read_spreadsheet('clone_sizes.xlsx',samples,channels);
% we provide the chase times:
chase_time = [30,25,35,70,70]; % days
empirical_tumour_expansion = [10.2, 7.8, 10.8, 14.1, 9.8];
% fix a value of the duplication probability, r, to break the fit
% degeneracy, the value must be larger that 1/2 and lower than 1
r = 0.75;
% Optional flag
% ignore_singlets : to ignore single cell clones, set to 1.  0 to keep them.
% Singles are ignored to consider only clones that originated from cells
% that were proliferative at the moment of induction.
ignore_singlets = 1;
% set show_plots to 1 to show results.
show_plots = 1;
%- SCRIPT ----------------------------------------------------------------%
% preallocate result table of parameters
variables = {'sample','channel','t0','num_clones','num_singlets',...
    'size_threshold','n0int','nbar','nbar_p','r','sigma','Delta_s',...
    'Delta_p','fs','predicted_tumour_volume','empirical_tumour_expansion',...
    'iT1_mT1','R2','S'};
result_table = array2table(zeros(length(samples)*length(channels),...
    length(variables)),'VariableNames',variables);
channel_labels = {'RFP','YFP'};
% run the analysis on each sample
counter = 0;
sensitivity_r = {}; dDelta_s = []; dsigma = []; dfs = [];
for nsample = samples
    if show_plots == 1
        figure;
        hold on;
    end
    mincdf = 1; % this is for adjusting the y-limis of the figure
    % run through each channel
    for nchannel = channels
        % clones to analyse
        sample_clones = clone_sizes{nsample,nchannel};
        num_clones = length(sample_clones); % number of clones
        num_singlets = sum(sample_clones == 1); % number of singlets
        % delete singles if ignore_singlets == 1
        if ignore_singlets == 1
            sample_clones(sample_clones == 1) = [];
        end
        % Optimization step to find the bet size_threshold in order to fit
        % the long term exponential decay
        initial_size_threshold = 20; % initial guess
        size_threshold_min = minimize(sample_clones,ignore_singlets);
        % Fit the biexponential model, Eq. (10), using the  value for
        % size_threshold found in the previous step, and extract the
        % parameters: nbar,nbar_p and p0p
        [gof,nbar,nbar_p,n0int,T1p] = fit_biexponential(sample_clones,size_threshold_min,ignore_singlets);
        % obtain the empirical and theoretical cdfs
        [emp_cdf,bincents] = empirical_cdf(sample_clones,show_plots,nchannel);
        Cn = theoretical_cdf(nbar,nbar_p,n0int,ignore_singlets,show_plots);
        % compute the goodness-of-fit measures: R2,S
        [R2,S] = gof_estimation(bincents,emp_cdf,Cn);
        
        %% ----------------------------------------------------------------
        % here we store all the fitted parameters in the results table:
        counter = counter+1;
        result_table.sample(counter) = nsample;
        result_table.channel(counter) = nchannel;
        t0 = chase_time(nsample); % chase time
        result_table.t0(counter) = t0;
        % number of clones
        result_table.num_clones(counter) = num_clones;
        % number of single cell clones
        result_table.num_singlets(counter) = num_singlets;
        % optimal clone size threshold
        result_table.size_threshold(counter) = size_threshold_min;
        % fitted parameters; n=0 intersect p0/(1-s(t)), nbar and nbar_p
        result_table.n0int(counter) = n0int; % p0/(1-s(t))
        result_table.nbar(counter) = nbar;
        result_table.nbar_p(counter) = nbar_p;
        % we now estimate other model parameters based on the theoretical
        % expression described in the "Model fits" section
        % calculate r-dependent parameters
        [Delta_s,sigma,fs,V] = estimate_remaining_parameters(r,t0,nbar,nbar_p,n0int,T1p);
        % SC cycling rate
        result_table.r(counter) = r;
        % stem cell expansion rate
        result_table.Delta_s(counter) = Delta_s;
        % stem cell cycling rate
        result_table.sigma(counter) = sigma;
        % stem cell expansion rate
        Delta_p = log(nbar_p)/t0;
        result_table.Delta_p(counter) = Delta_p;
        % T1p being the n=1 intersect of the small clone exponential decay
        result_table.fs(counter) = fs;
        % infered T1 - measured T1
        result_table.iT1_mT1(counter) = num_singlets/num_clones - T1p;
        % tumor volume, Eq. (11)
        result_table.predicted_tumour_volume(counter) = V;
        % empirical_tumour_expansion
        result_table.empirical_tumour_expansion(counter) = empirical_tumour_expansion(nsample);
        % goodness-of-fit estimates
        result_table.R2(counter) = R2;
        result_table.S(counter) = S;
        % -- Study sensitivity to variations of "r"
        rcounter = 0;
        for rvals = 0.6:0.01:1
            rcounter = rcounter + 1;
            [Delta_s,sigma,fs,~] = estimate_remaining_parameters(rvals,t0,nbar,nbar_p,n0int,T1p);
            sensitivity_r{nsample,nchannel}(rcounter,:) = [rvals,Delta_s,sigma,fs];
        end
        % lets calculate the relative deviation of Delta_s,sigma,fs given a
        % variation of 10% in the value of r
        r_range = [0.75,0.9]; % 10% around 0.75
        [dDelta_s(end+1,1),dsigma(end+1,1),dfs(end+1,1)] = variation_of_parameters(r,r_range,t0,nbar,nbar_p,n0int,T1p);
        % -----------------------------------------------------------------
        mincdf = min([mincdf,emp_cdf(emp_cdf>0)]);
    end
    if show_plots == 1 % plot adjustments
        hold off
        % adjust axes
        axis([0 max(vertcat(clone_sizes{nsample,1:2})) min([mincdf,10^-2]) 1])
        xlabel(['Clone size (cells)'])
        ylabel('Cumulative probability')
        title(['sample ',num2str(nsample)])
        set(gca,'YScale','log'); % set y-axis to log scale
        legend(channel_labels)
        legend boxoff
        set(gcf,'color','w');
    end
end
%%
% calculate and show average deviation of the parameters
mean([dDelta_s,dsigma,dfs],1)

if show_plots == 1 % plot sensitivity to variations of "r"
    figure;
    r_ref = 0.75;
    legs = {}; 
    markers = {'o','s','^','v','x'}; lstyle = {'-','--'};
    for nsample = samples
        for nchannel = channels
            [c,~] = candm(nchannel);
            rs = sensitivity_r{nsample,nchannel}(:,1);
            mls = [lstyle{nchannel},markers{nsample}];
            % plot SC expansion rate
            subplot(1,3,1); hold on
            plot(rs,sensitivity_r{nsample,nchannel}(:,2),mls,'Color',c); 
            add_vref_line(r_ref); hold off;
            xlabel('SC duplication Prob. (r)'); ylabel('SC expansion rate (\Delta_s)'); 
            xlim([min(rs),max(rs)]);
            % plot SC cycling rate
            subplot(1,3,2); hold on
            plot(rs,sensitivity_r{nsample,nchannel}(:,3),mls,'Color',c); 
            add_vref_line(r_ref); hold off;
            xlabel('SC duplication Prob. (r)'); ylabel('SC cycling rate (\sigma)'); 
            xlim([min(rs),max(rs)]); 
            % plot SC fraction 
            subplot(1,3,3); hold on
            plot(rs,sensitivity_r{nsample,nchannel}(:,4),mls,'Color',c); 
            add_vref_line(r_ref); hold off;
            xlabel('SC duplication Prob. (r)'); ylabel('SC fraction (f_S)'); 
            xlim([min(rs),max(rs)]);
            set(gcf,'color','w');
            
            legs{end+1} = ['S',num2str(nsample),'-',channel_labels{nchannel}];
        end
    end
    legend(legs); legend boxoff
end
%%
% save results
save('result_table.mat','result_table')
disp('Data saved as result_table.mat')
%% FUNCTIONS
function [size_threshold_min, gof] = minimize(clone_sizes,ignore_singlets)
% Input:
% clone_sizes: list of clone sizes for a given sample
% find optimal size_threshold_min in the range [1,50] cells. This range was
% chosen as it was obseved than in all samples the transition between the
% two exponential regimes lied within this range.
[size_threshold_min, gof] = fminbnd(@(size_threshold) fit_biexponential(clone_sizes,size_threshold,ignore_singlets),1,50);
end

function [gof,nbar,nbar_p,p0p,T1,bincents] = fit_biexponential(clone_sizes,size_threshold,ignore_singlets)
    % ---------------------------------------------------------------------
    % This function fits the bi-exponential CDF function to the empirical
    % data as described in the "Model fits" section of the supplementary
    % Note
    % ---------------------------------------------------------------------
    % construct empirical CDF
    [emp_cdf,bincents] = empirical_cdf(clone_sizes,0);
    % extract values above size_threshold (long term dependence)
    cdf_long = emp_cdf(bincents > size_threshold);
    bincents_long = bincents(bincents > size_threshold);
    % ignore zeros in the CDF
    bincents_long(cdf_long == 0) = []; cdf_long(cdf_long == 0) = [];
    
    % define exponential function for fitting
    exp_fun = @(c, xi) c(1)*exp(-xi./c(2));
    
    % 1) fit exponential decay to the long term dependence
    % initial guess
    x0 = [1,size_threshold];
    % run fitting routine
    % provide reasonable lower (lb) and upper (ub) bounds to improve
    % convergence
    lb = [0 1]; ub = [1 max(clone_sizes)];
    c_long = lsqcurvefit(exp_fun, x0, bincents_long, cdf_long,lb,ub);
    % evaluate exponential function with the fitted parameters
    y_long = exp_fun(c_long,bincents_long);
    % ---------------------------------------------------------------------
    % extract nbar
    nbar = c_long(2);
    % extract the composite parameter p0p=fs*(2r-1)/r) from the n=0 intersect 
    p0p = c_long(1);
    % ---------------------------------------------------------------------
    % 2) Substract the long-term dependence from the empirical CDF and
    % estimate the short term decay
    y_long = exp_fun(c_long,bincents);
    cdf_short = emp_cdf - y_long;
    % ignore zeros in the CDF
    bincents_short = bincents(cdf_short > 0); cdf_short(cdf_short <= 0) = [];
    % extract values below size_threshold (short term dependence)
    % initial guess
    x0 = [1,size_threshold];
    % run fitting routine
    % provide reasonable lower (lb) and upper (ub) bounds to improve
    % convergence
    lb = [0 1]; ub = [1 max(clone_sizes)];
    [c_short] = lsqcurvefit(exp_fun, x0, bincents_short, cdf_short,lb,ub);
    % evaluate exponential function with the fitted parameters
    y_short = exp_fun(c_short,bincents_short);
    T1 = 1-y_short(1); % y_short(1) corresponds to the Prob. that a clone has than 1 cell.
    % ---------------------------------------------------------------------
    % extract nbar_p
    nbar_p = c_short(2);
    % ---------------------------------------------------------------------
    % evaluate the theoretical cdf, and extract the goodness-of-fit
    % parameter S.
    Cn = theoretical_cdf(nbar,nbar_p,p0p,ignore_singlets,0);
    [~,S] = gof_estimation(bincents,emp_cdf,Cn);
    gof = S;
end 

function [R2,S] = gof_estimation(bincents,emp_cdf,theo_cdf)
    % we just compare clone sizes where we do have experimental data
    theo_cdf = theo_cdf(isfinite(emp_cdf));
    emp_cdf = emp_cdf(isfinite(emp_cdf));
    ybar = nanmean(emp_cdf);
    Stot = nansum((emp_cdf - ybar).^2);
    Sres = nansum((emp_cdf - theo_cdf).^2);
    R2 = 1 - Sres/Stot;
    S = sqrt(Sres/length(bincents));
end

function Cn = theoretical_cdf(nbar,nbar_p,n0int,ignore_singlets,show_plot)
    % Here we contuct the theoretical cdf based on the input parameters
    % nbar, nbar_p and p0p.
    % If ignore_singlets == 1, then T{n<2} is removed from Tn(t).
    % provide range of clone sizes to evaluate
    bincents = 0.5:1:10^5;
    % compute T(n), Eq. (10), and normalize
    Tn = @(xi) n0int*exp(-xi./nbar)./nbar+(1-n0int)*exp(-xi./nbar_p)./nbar_p;
    if ignore_singlets == 1
        Tn_list = Tn(bincents)/(1-sum(Tn(0:1)));
    else
        Tn_list = Tn(bincents);
    end
    % compute cdf, C(n):
    Cn = 1 - cumsum(Tn_list)./sum(Tn_list);
    if show_plot == 1
       plot(bincents,Cn,'--k','HandleVisibility','off') 
    end
end

function [cdf,bincents] = empirical_cdf(clone_sizes,show_plots,nchannel)
% Here we contuct the empirical cdf
% define bins (1 cell spacing)
binedges = 0.5:1:max(clone_sizes)*1+0.5; bincents = (binedges(2:end) + binedges(1:end - 1))/2;
% compute pdf
Nnorm = histcounts(clone_sizes,binedges,'Normalization','pdf');
% compute cdf
cdf = 1 - cumsum(Nnorm)./sum(Nnorm);
if show_plots == 1
    [c,m] = candm(nchannel);
    plot(bincents,cdf,m,'color',c,'MarkerFaceColor','w','MarkerSize',5)
end
end

function clone_sizes = read_spreadsheet(file_path,samples,channels) % read spreadsheet
clone_sizes_xls =  xlsread(file_path); % load numeric values from excel sheet
clone_sizes_xls(:,1) = []; % first row has just clone indices;
clone_sizes = {};
for nsample = samples
    for nchannel = channels
        clone_sizes_sample = clone_sizes_xls(:,2*nsample-1+nchannel-1);
        clone_sizes_sample(~isfinite(clone_sizes_sample)) = []; % this is to remove nan create by loading empty cells.
        clone_sizes{nsample,nchannel} = clone_sizes_sample;
    end
end
end

function [c,m] = candm(n)
% define colours and markers for plots
cols = {[1,0,0],[0.9290, 0.6940, 0.1250]};
mkrs = {'-o','-s'};
c = cols{n};
m = mkrs{n};
end

function [Delta_s,sigma,fs,V] = estimate_remaining_parameters(r,t0,nbar,nbar_p,n0int,T1p)
% calculate parameters for a given value of "r":
factor = (2*r-1)/r;
% stem cell expansion rate
Delta_s = (1/t0)*log(factor^2*nbar);
% stem cell cycling rate
sigma = Delta_s/(2*r-1);
% T1p being the n=1 intersect of the small clone exponential decay
fs = n0int*(1-T1p)/factor; % from fsp = fs/(1-s(t))
% tumor volume, Eq. (11)
V = tumour_volume_SP_model(r,nbar,nbar_p,fs);
end

function V = tumour_volume_SP_model(r,nbar,nbar_p,fs)
% calculate tumour volume from fitted parameters nbar,nbar_p,fs, using the
% fixed value of r, Eq. (11):
V = fs*(2*r-1)/r*nbar + (1-(2*r-1)/r*fs)*nbar_p;
end

function add_vref_line(x)
ys = ylim();
xs = [x,x];
plot(xs,ys,'k--','HandleVisibility','off')
end

function [dDelta_s,dsigma,dfs] = variation_of_parameters(r,r_range,t0,nbar,nbar_p,n0int,T1p)
[Delta_s,sigma,fs,~] = estimate_remaining_parameters(r,t0,nbar,nbar_p,n0int,T1p);
[Ds_sup,s_sup,fs_sup,~] = estimate_remaining_parameters(r_range(1),t0,nbar,nbar_p,n0int,T1p);
[Ds_inf,s_inf,fs_inf,~] = estimate_remaining_parameters(r_range(2),t0,nbar,nbar_p,n0int,T1p);
% relative deciation
dDelta_s = abs(Ds_sup-Ds_inf)/Delta_s;
dsigma = abs(s_sup-s_inf)/sigma;
dfs = abs(fs_sup-fs_inf)/fs;
end
