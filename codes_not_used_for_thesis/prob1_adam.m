delta = 1;
tau = 2;
r = 101; % to check e=0
t = zeros(r,1);

for idx = 1:r
    t(idx) = (idx-1)/(r-1);
end

e = zeros(r,1);
for idx = 1:r
    e(idx) = -1 + 2*(idx-1)/(r-1);
end

expt1 = zeros(r);
for idx = 1:r
    for jdx = 1:r
        expt1(idx, jdx) = exp(-1*1i*t(jdx)*e(idx));
    end
end

expt2 = exp(-1i*tau*e);

A = expt1;
b = -expt2;
%%
cvx_begin
    variables c(r);
    minimize(norm(c,1)+delta*max(abs(A*c+b)));
cvx_end
%% 
% subgradient method

f_min = cvx_optval;
fprintf(1,'Optimal value is %0.4f.\n\n',f_min);

% initial point
x_1 = ones(r,1);
Rtrue = norm( x_1 - c );
R = 10; % we use a much bigger R than it is needed (heuristic)

% constant step length
TOL = 1e-3;
MAX_ITERS = 3000; 
gammas = [.05 .005];
label1 = sprintf('%.3f',gammas(1));
label2 = sprintf('%.3f',gammas(2));

% run subgradient method with constant step length for different gammas
[x,hist1] = sgm_pwl_const_step_length(A,b,x_1,R,gammas(1),TOL,MAX_ITERS); 
[x,hist2] = sgm_pwl_const_step_length(A,b,x_1,R,gammas(2),TOL,MAX_ITERS); 
%%
%********************************************************************
% generate plots
%********************************************************************
% setup plot data
iters = [1:MAX_ITERS];
f1 = hist1{1}; fbest1 = hist1{2}; lbest1 = hist1{3};
f2 = hist2{1}; fbest2 = hist2{2}; lbest2 = hist2{3};
%%
% plots
iter_sm = 1000;
semilogy(2:iter_sm, f1(2:iter_sm)-f_min, 'r-','LineWidth',1.5 ); hold on;
semilogy(2:iter_sm, f2(2:iter_sm)-f_min, 'g-','LineWidth',1.5 ); hold off 
xlabel('k');
ylabel('f - fmin');
legend(label1, label2,'Location','northeast')