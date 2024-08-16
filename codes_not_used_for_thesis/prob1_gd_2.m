%% A. HS unconstrained
% variables

omega = 500;
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

% expt1 = zeros(r);
% for idx = 1:r
%     for jdx = 1:r
%         expt1(idx, jdx) = exp(-1*1i*t(jdx)*e(idx));
%     end
% end
% 
% expt2 = exp(-1i*tau*e);
%% 
% further variables

m1 = zeros(r); M1 = cell(r,1);
m4 = zeros(r,1); M4 = cell(r,1); m5 = zeros(r); M5 = cell(r,1);

for idx = 1:r %e
    temp_e = e(idx);
    for jdx = 1:r %j
        for kdx = 1:r %k
            m1(jdx,kdx) = cos(temp_e*(t(jdx)-t(kdx)));
        end
    end
    M1{idx} = m1;
    m1 = zeros(r);
end

for idx = 1:r %e
    for jdx = 1:r %j
        m5(idx, jdx) = cos(e(idx)*(t(jdx)-tau));
    end
end
%% 
% cvx

cvx_begin
    variables c(r); 
    for idx = 1:r
        M4{idx} = c'*M1{idx}*c;
    end

    for idx = 1:r
        M5{idx} = 2*m5(idx,:)*c;
    end
    
    maxobj = cellfun(@minus,M4,M5,'Un',0); %M4-2*M5+1
    maxobj = cellfun(@(x) x+1,maxobj,'un',0);
    
    minimize(norm(c,1)+omega*max([maxobj{:}]));
    
cvx_end

delta = max([maxobj{:}])
[fval,ind] = max([maxobj{:}]);
Estar = e(ind);

% trials and errors
%m3 = c*c'; % trial1 not accepted
%m3 = c'*M1{idx}*c; % trail2 cannot put cvx objs into an array
%m4(idx) = m3;
%maxobjmat = cell2mat(maxobj); % trial3 cell2mat does not support cells containing objects
%cellfun(@(x) max(cell2mat(x)),A)
%B = [A{:}] % another way to transform cell to array
%% 
% subgradient method

f_min = cvx_optval;
fprintf(1,'Optimal value is % 0.4f.\n\n',f_min);

x_1 = 0.001*ones(r,1); % start

MAX_ITERS = 3000; 

% run subgradient method with diminishing step sizes
% consider different step sizes
as = [.003 .002 .001];
%[x,hist1] = sgm_pwl_nonsum_dimin(A,b,x_1,as(1),MAX_ITERS);
[x,hist1] = sgm_pwl_sqrsum_nonsum(M1,m5,e,t,r,tau,omega,c,x_1,as(1),MAX_ITERS); 
[x,hist2] = sgm_pwl_sqrsum_nonsum(M1,m5,e,t,r,tau,omega,c,x_1,as(2),MAX_ITERS); 
[x,hist3] = sgm_pwl_sqrsum_nonsum(M1,m5,e,t,r,tau,omega,c,x_1,as(3),MAX_ITERS); 
label1 = sprintf('%.3f/k',as(1));
label2 = sprintf('%.3f/k',as(2));
label3 = sprintf('%.3f/k',as(3));

% run subgradient method with Polyak's optimal step
%[x,histo] = sgm_pwl_optimal_step(A,b,x_1,f_min,MAX_ITERS); 

% setup plot data
iters = [1:MAX_ITERS];
f1 = hist1{1}; fbest1 = hist1{2};
f2 = hist2{1}; fbest2 = hist2{2};
f3 = hist3{1}; fbest3 = hist3{2};
%fo = histo{1}; fbesto = histo{2};
%%
% plot
set_iters = 3000;
figure(1), clf
set(gca, 'FontSize',18);
semilogy(2:set_iters, fbest1(2:set_iters)-f_min, 'b','LineWidth',1.5 ); hold on; %semilogy/plot
semilogy(2:set_iters, fbest2(2:set_iters)-f_min, 'r','LineWidth',1.5 ); 
semilogy(2:set_iters, fbest3(2:set_iters)-f_min, 'g','LineWidth',1.5 ); 
hold off
xlabel('k');
ylabel('fbest - f*');
legend(label1, label2, label3);