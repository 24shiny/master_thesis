delta = 0.01;
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
%% 
% CVX - why fail?

A = expt1;
b = -expt2;

cvx_begin
    variable c(r); 
    minimize(norm(c,1)+delta*max(abs(A*c+b)));
cvx_end
%% 
% solve non diff LP


%% 
% subgradient method

f_min = cvx_optval;
fprintf(1,'Optimal value is %0.4f.\n\n',f_min);

x_1 = 0.5*ones(r,1); % start

MAX_ITERS = 3000; 

% run subgradient method with diminishing step sizes
% consider different step sizes
as = [.1 .08 .07];
%[x,hist1] = sgm_pwl_nonsum_dimin(A,b,x_1,as(1),MAX_ITERS);
[x,hist1] = sgm_pwl_sqrsum_nonsum(A,b,x_1,as(1),MAX_ITERS); 
[x,hist2] = sgm_pwl_sqrsum_nonsum(A,b,x_1,as(2),MAX_ITERS); 
[x,hist3] = sgm_pwl_sqrsum_nonsum(A,b,x_1,as(3),MAX_ITERS); 
label1 = sprintf('%.2f/k',as(1));
label2 = sprintf('%.2f/k',as(2));
label3 = sprintf('%.2f/k',as(3));

% run subgradient method with Polyak's optimal step
[x,histo] = sgm_pwl_optimal_step(A,b,x_1,f_min,MAX_ITERS); 

% setup plot data
iters = [1:MAX_ITERS];
f1 = hist1{1}; fbest1 = hist1{2};
f2 = hist2{1}; fbest2 = hist2{2};
f3 = hist3{1}; fbest3 = hist3{2};
fo = histo{1}; fbesto = histo{2};
%%
set_iters = 2000;
% plots
figure(1), clf
set(gca, 'FontSize',18);
plot(2:set_iters, fbest1(2:set_iters)-f_min, 'b','LineWidth',1.5 ); hold on; %semilogy
plot(2:set_iters, fbest2(2:set_iters)-f_min, 'r','LineWidth',1.5 ); 
plot(2:set_iters, fbest3(2:set_iters)-f_min, 'g','LineWidth',1.5 ); 
hold off
xlabel('k');
ylabel('fbest - fmin');
legend(label1, label2, label3);
%print -depsc pwl_dimin_step_fbest