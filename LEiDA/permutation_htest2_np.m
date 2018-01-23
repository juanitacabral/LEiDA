function [stats] = permutation_htest2_np(data,design,niter,pthr,htest)
% PERMUTATION_HTEST2_NP - A "non parametric" two-sample hypotesis test that, instead of
% relying on the test-type standard distribution, uses permutations of group labels to
% estimate the null distribution. The null distribution is computed
% independently for each data point (= row), i.e. we do not assume the same
% distribution for each datapoint. However, we do assume that the data
% points are comparable (e.g. they correspond to the same location
% collected across all subjects)
%
% USAGE:
%   stats = bramila_ttest2_np(data,design,niter)
% INPUT:
%   data   - a matrix where each column is a subject and each row is a
%            data-point for example a voxel intensity in fMRI, a node level
%            value in a network, etc. NaN values will be ignored.
%   design - a row vector containing the numbers 1 and 2 for the two groups
%   niter  - number of permutations (recommended 5000)
%   htest  - hypothesis test used to compare populations. The script is
%            prepared to run the ttest2, kstest2, and ranksum tests. 
%
%
% OUTPUT:
%   stats is a struct with the following subfields:
%       pvals - p-values for each datapoint; it returns in order the p-values
%               for the right tail and for the left tail
%       tvals - test statistic values for datapoint, positive tvals mean 
%               group 1 > group 2
%
% Notes: the null distribution is estimated using the matlab function
% ksdensity by interpolating the permuted data. The distribution is
% estimated over 200 points if niter<=5000, otherwise it is estimated over
% round(200*niter/5000) points, for greater precision.

%  Henrique Fernandes 2014
%  adapted from: Enrico Glerean 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT VALIDATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nsubj = size(data,2);   % number of subjects
if(size(design,2) ~= Nsubj)
    error('Mismatched number of subjects: the number of columns of data variable  should match the number of columns of the design variable.')
end
if(size(design,1) ~= 1)
    error('The design variable should only contain 1 row')
end

g1 = find(design==1);
g2 = find(design==2);
if((length(g1)+length(g2))~=Nsubj)
    error('The design variable should only contain numbers 1 and 2.')
end

if(niter<=0)
    disp('The variable niter should be a positive integer, function will continue assuming niter=5000.')
    niter=5000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYPOTHESIS TESTING (for each row/area)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stats.tvals=tt_np(data,g1,g2); % similar to ttest2

NC    = size(data,1); % number of comparisons
tvals = zeros(NC,1);
means = zeros(NC,1);

switch htest
    case 'ttest' 
        % - the population means are not equal. (alternative hypothesis)
        % - the two groups are derived from normal distributions with unknown and unequal variances.
        for t=1:NC
            [H,P,CI,STATS] = ttest2(data(t,g1)',data(t,g2)',pthr,'both','unequal');
            tvals(t,:)     = STATS.tstat;
            diffs(t,:)      = mean(data(t,g1))-mean(data(t,g2));
        end
    case 'kstest'
        for t=1:NC
            [H,P,STATS]=kstest2(data(t,g1)',data(t,g2)',pthr);
            tvals(t,:)=STATS;
            diffs(t,:)      = mean(data(t,g1))-mean(data(t,g2));
        end
    case 'ranksum'
        for t=1:NC
            [P,H,STATS]=ranksum(data(t,g1)',data(t,g2)','alpha',pthr);
            tvals(t,:)=STATS;
            diffs(t,:)      = mean(data(t,g1))-mean(data(t,g2));
        end
    otherwise
        error('\n-------------------------------\n\nHypothesis test %s not recognized. \n\n-------------------------------\n',htest)
end
    
stats.tvals         = tvals;
% tvals(isnan(tvals)) = 0;   % or tvals(tvals~=tvals) = 0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERMUTATION TESTING (for each row/area)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outptus the pval (from the computed null distribution using permutation 
% testing) given the tstat previously calculated.
% each comparison is treated independently

pvals=zeros(NC,2);

parfor n=1:NC
    if median(data(n,g1))~=0 || median(data(n,g2))~=0  % Exclude tests where all (tstat=NaN) or most of the population (median=0) as a null value.
        pvals(n,:) = tt_np_pval(data(n,:),g1,g2,niter,tvals(n),pthr);
    else
        pvals(n,:) = [NaN NaN];
    end
end

stats.pvals=pvals;

stats.diffs=diffs;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tval=tt_np(data,g1,g2)
    % helper function similar to matlab function ttest2.m for the case of
    % groups with difference variance
    
    xnans = isnan(data(:,g1));
    if any(xnans(:))
        nx = sum(~xnans,2);
    else
        nx = size(data(:,g1),2); 
    end
    ynans = isnan(data(:,g2));
    if any(ynans(:))
        ny = sum(~ynans,2);
    else
        ny = size(data(:,g2),2); % a scalar, => a scalar call to tinv
    end

    difference = nanmean(data(:,g1),2) - nanmean(data(:,g2),2);
    
    s2x = nanvar(data(:,g1),[],2);
    s2y = nanvar(data(:,g2),[],2);
    s2xbar = s2x ./ nx;
    s2ybar = s2y ./ ny;
    se = sqrt(s2xbar + s2ybar); % standard error of the estimated coefficient 'difference'.
    if(any(se == 0) || any(isnan(se)))
        error('Group variance seems to be null or NaN, please check your data')
    end
    tval = difference ./ se;

end

function pval = tt_np_pval(data,g1,g2,niter,tval,pthr)
    
    outiter=zeros(niter,1);
    ND=length(data);
    for iter=1:niter
        perm=randperm(ND);
        % one could add a test to see that they are indeed permuted
        temp=data(perm);
        [H,P,CI,STATS]=ttest2(temp(:,g1)',temp(:,g2)',pthr,'both','unequal');
        outiter(iter)=STATS.tstat;
%          outiter(iter)=tt_np(temp,g1,g2); % t-stats for this permutation test
    end
    
    NCDF=200;
    if(niter>5000)
        NCDF=round(200*niter/5000);
    end
    [fi xi]=ksdensity(outiter,'function','cdf','npoints',NCDF); % estimated cumulative distribution function

%     figure;
%     subpfd=subplot(1,2,1); hold on;    vline(tval,'g','tval');
%     subcdf=subplot(1,2,2); hold on;    vline(tval,'g','tval');
%     [h, fhat, xgrid] = kde(outiter);
%     [fi2 xi2]=ksdensity(outiter,'npoints',NCDF); % estimated probability density function
%     plot(subpfd,xgrid, fhat, 'linewidth', 1.7, 'color', 'blue');
%     plot(subpfd,outiter,zeros(size(outiter,2),1), 'r+');
%     set(subpfd,'box','off');
%     axis(subpfd,[min(xi2),max(xi2),0,max(fi2)+0.05])
%     xlabel(subpfd,'t-stats permuted values')
%     ylabel(subpfd,'Probability density function')
%     plot(subcdf,xi,fi,'linewidth', 1.7, 'color', 'blue')
%     axis(subcdf,[min(xi2),max(xi2),0,1]);
%     xlabel(subcdf,'t-stats permuted values')
%     ylabel(subcdf,'Cumulative distribution function')

    % trick to avoid NaNs, we approximate the domain of the CDF between
    % -Inf and Inf using the atanh function and the eps matlab precision
    % variable
    
    pval_left=interp1([atanh(-1+eps) xi atanh(1-eps)],[0 fi 1],tval); % G1 > G2
    pval_right=1-pval_left; % G1 < G2
    pval=[pval_right pval_left];
end
 