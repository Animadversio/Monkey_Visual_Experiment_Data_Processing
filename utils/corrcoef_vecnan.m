function [r,p] = corrcoef_vecnan(x,y)
% function [r,p] = corrcoef_vecnan(x,y)
% estimate correlation coefficient between a vector x and the columns of a
% matrix y
% Example: x is of size (81,1), y of size (81, 6000). 
% adapted from corrcoef_vec to support nan values

tail = 'b';
%%
n = length(x);
nanmskx = isnan(x);
nanmsky = isnan(y);
xc = x - nanmean(x,1); % Remove mean 
dx = sqrt(nansum(abs(xc).^2,1));
%xc = bsxfun(@minus,x,sum(x,1)/n);  % Remove mean
%dx = sqrt(sum(abs(xc).^2, 1)); % sqrt first to avoid under/overflow

%%
yc = y - nanmean(y,1); % Remove mean
dy = sqrt(nansum(abs(yc).^2, 1));
% yc = bsxfun(@minus,y,sum(y,1)/n);  % Remove mean
% dy = sqrt(sum(abs(yc).^2, 1)); % sqrt first to avoid under/overflow

r = xc' * yc; % 1/(n-1) doesn't matter, renormalizing anyway
r = bsxfun(@rdivide,r,dx');
r = bsxfun(@rdivide,r,dy); % r = r ./ dx'*dy;

if nargin>1
    %PVALPEARSON Tail probability for Pearson's linear correlation.
    t = r.*sqrt((n-2)./(1-r.^2)); % +/- Inf where r == 1
    switch tail
        case 'b' % 'both or 'ne'
            p = 2*tcdf(-abs(t),n-2);
        case 'r' % 'right' or 'gt'
            p = tcdf(-t,n-2);
        case 'l' % 'left' or 'lt'
            p = tcdf(t,n-2);
    end
end

