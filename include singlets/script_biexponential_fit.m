%%-----------------------------------------------------------------------%%
%%- PARAMETER FITS FOR THE TWO-COMPARTMENT (SP) MODEL -------------------%%
%%-----------------------------------------------------------------------%%
% This script loads the clone sizes, construct the cumulative distribution
% of clone sizes and estimates the parameters of the bi-exponential
% distribution, Supplemental Note Eq. (10), following the strategy
% described in the "Model fits" section.
%
% In this script we generate:
% Figure 2F, Supplementary Figures 5b-5f, and Table S4.
%%-----------------------------------------------------------------------%%
close all; clear all; clc;
%-DATA AND PARAMETERS-----------------------------------------------------%
% Read spreadsheet as a cell array: 5 samples, 2 channels each
samples = 1:5; channels =1:2; % 1-RFP, 2-YFP
clone_sizes = read_spreadsheet('clone_sizes.xlsx',samples,channels);
% we provide the chase times:
chase_time = [30,25,35,70,70]; % days
empirical_tumour_expansion = [10.2, 7.8, 10.8, 14.1, 9.8];
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
    'size_threshold','p0','nbar','nbar_p','r','sigma','Delta_s',...
    'Delta_p','fs','predicted_tumour_volume','empirical_tumour_expansion','R2','S'};
result_table = array2table(zeros(length(samples)*length(channels),...
    length(variables)),'VariableNames',variables);
channel_labels = {'RFP','YFP'};
% run the analysis on each sample
counter = 0;
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
        size_threshold_min = minimize(sample_clones);
        % Fit the biexponential model, Eq. (10), using the  value for
        % size_threshold found in the previous step, and extract the
        % parameters: nbar,nbar_p and p0
        [gof,nbar,nbar_p,p0] = fit_biexponential(sample_clones,size_threshold_min);
        % obtain the empirical and theoretical cdfs
        [emp_cdf,bincents] = empirical_cdf(sample_clones,show_plots,nchannel);
        Cn = theoretical_cdf(nbar,nbar_p,p0,show_plots);
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
        % fitted parameters p0, nbar and nbar_p
        result_table.p0(counter) = p0;
        result_table.nbar(counter) = nbar;
        result_table.nbar_p(counter) = nbar_p;
        % we now estimate other model parameters based on the theoretical
        % expression described in the "Model fits" section
        % first we fix a value of the duplication probability, r.
        r = 0.75; factor = (2*r-1)/r;
        result_table.r(counter) = r;
        % stem cell expansion rate
        Delta_s = (1/t0)*log(factor^2*nbar);
        result_table.Delta_s(counter) = Delta_s;
        % stem cell cycling rate
        sigma = Delta_s/(2*r-1);
        result_table.sigma(counter) = sigma;
        % stem cell expansion rate
        Delta_p = log(nbar_p)/t0;
        result_table.Delta_p(counter) = Delta_p;
        % fraction fo induces stem cells
        fs = p0/factor;
        result_table.fs(counter) = fs;
        % tumor volume, Eq. (11)
        V = tumour_volume_SP_model(r,nbar,nbar_p,fs);
        result_table.predicted_tumour_volume(counter) = V;
        % empirical_tumour_expansion
        result_table.empirical_tumour_expansion(counter) = empirical_tumour_expansion(nsample);
        % goodness-of-fit estimates
        result_table.R2(counter) = R2;
        result_table.S(counter) = S;
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
% save results
save('result_table.mat','result_table')
disp('Data saved as result_table.mat')
%% FUNCTIONS
function [size_threshold_min, gof] = minimize(clone_sizes)
% Input:
% clone_sizes: list of clone sizes for a given sample
% find optimal size_threshold_min in the range [1,50] cells. This range was
% chosen as it was obseved than in all samples the transition between the
% two exponential regimes lied within this range.
[size_threshold_min, gof] = fminbnd(@(size_threshold) fit_biexponential(clone_sizes,size_threshold),1,50);
end

function [gof,nbar,nbar_p,p0,bincents] = fit_biexponential(clone_sizes,size_threshold)
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
    % extract the composite parameter p0=fs*(2r-1)/r) from the n=0 intersect 
    p0 = c_long(1);
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
    % ---------------------------------------------------------------------
    % extract nbar_p
    nbar_p = c_short(2);
    % ---------------------------------------------------------------------
    % evaluate the theoretical cdf, and extract the goodness-of-fit
    % parameter S.
    Cn = theoretical_cdf(nbar,nbar_p,p0,0);
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

function Cn = theoretical_cdf(nbar,nbar_p,p0,show_plot)
    % Here we contuct the theoretical cdf based on the input parameters
    % nbar, nbar_p and p0.
    % provide range of clone sizes to evaluate
    bincents = 0.5:1:10^5;
    % compute T(n), Eq. (10), and normalize
    Tn = @(xi) p0*exp(-xi./nbar)./nbar+(1-p0)*exp(-xi./nbar_p)./nbar_p;
    Tn_norm = Tn(bincents)./sum(Tn(bincents));
    % compute cdf, C(n):
    Cn = 1 - cumsum(Tn_norm);
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

function V = tumour_volume_SP_model(r,nbar,nbar_p,fs)
% calculate tumour volume from fitted parameters nbar,nbar_p,fs, using the
% fixed value of r, Eq. (11):
V = fs*(2*r-1)/r*nbar + (1-(2*r-1)/r*fs)*nbar_p;
end


