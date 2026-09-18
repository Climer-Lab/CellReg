function progress_tick(k,n)
% This function updates the text progress bar (display_progress_bar) for
% iteration k of n, but only when the displayed percentage changes, so that
% tight loops are not dominated by printing.

% Inputs:
% 1. k - current iteration (1..n)
% 2. n - number of iterations

persistent last_percent
if k<=1
    last_percent=-1;
end
percent=floor(100*k/n);
if percent~=last_percent || k==n
    last_percent=percent;
    display_progress_bar(100*k/n,false)
end

end
