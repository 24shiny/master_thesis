%% *A. HS extrapolation functions*
% Given conditions

r = 100; %resolution
tolval = 0.00000001;
tol = tolval*ones(r,1); %tolerance
tau = 2;
% Define variables

e = zeros(r,1);
for idx = 1:r
    e(idx) = -1 + 2*(idx-1)/(r-1);
end

t = zeros(r,1);
for idx = 1:r
    t(idx) = (idx-1)/(r-1);
end

expt1 = zeros(r);
for idx = 1:r
    for jdx = 1:r
        expt1(idx, jdx) = exp(-1*1i*t(jdx)*e(idx));
    end
end

expt2 = exp(-1i*tau*e);
% Define and solve the problem

cvx_begin
    variables c(r);
    minimize(norm(c,1));
    subject to
        abs(expt1*c - expt2) <= tol;
cvx_end
% Plot the solution

scatter(1:r,c)
xlabel('j')
ylabel('c_j')
label1 = sprintf('delta = %.8f, tau=%.1f',tolval,tau);
title(label1);
cvx_optval
%% 
% additional data for checking values of e

% LHS = abs(expt1*c - expt2);
% sz = 10*ones(1,r);
% scatter(e,LHS,sz); hold on
% for idx = 1:r
%     if LHS(idx) == tolval
%         scatter(e(idx), LHS(idx), 13,"Red",'filled'); hold on;
%     elseif LHS(idx) > tolval*0.99999
%         scatter(e(idx), LHS(idx), 13, "Red",'filled'); hold on;
%     end
% end
% plot(e,tol,'r-'); hold off
% xlabel('e');
% ylabel('|del(E)|')
% label1 = sprintf('delta = %.4f, tau=%.1f',tolval,tau);
% title(label1);