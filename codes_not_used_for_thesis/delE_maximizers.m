%% 
% variables

omega = 1000;
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

m1 = zeros(r); M1 = cell(r,1);
m5 = zeros(r); M5 = cell(r,1);

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
% gradient method

x1 = 0.1*ones(r,1); % start

MAX_ITERS = 1000; 

% with diminishing step sizes
a = [0.001, 0.0001, 0.00001];
[c1,delta1,f1,ind1] = gradient_descent_max(M1,m5,e,t,r,tau,omega,x1,a(1),MAX_ITERS); 
[c2,delta2,f2,ind2] = gradient_descent_max(M1,m5,e,t,r,tau,omega,x1,a(2),MAX_ITERS); 
[c3,delta3,f3,ind3] = gradient_descent_max(M1,m5,e,t,r,tau,omega,x1,a(3),MAX_ITERS); 

label1 = sprintf('%.3f/k',a(1));
label2 = sprintf('%.4f/k',a(2));
label3 = sprintf('%.5f/k',a(3));
%% 
% plot

set_iters = 1:MAX_ITERS-1;
semilogy(set_iters, f1(set_iters), 'b','LineWidth',1.5 ); hold on; %semilogy/plot
semilogy(set_iters, f2(set_iters), 'r','LineWidth',1.5 ); 
semilogy(set_iters, f3(set_iters), 'g','LineWidth',1.5 ); 
hold off
xlabel('k');
ylabel('f(c)');
legend(label1, label2, label3);
%% 
% min(f(c))

A = [f1(end), f2(end), f3(end)]
B = [delta1(end), delta2(end), delta3(end)]
[obj_omega, idx] = min(A)
c_ary = [c1 c2 c3]; c = c_ary(:,idx)
c_1norm = norm(c,1)
delE = min(B)