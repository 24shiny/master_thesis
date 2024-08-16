%% 
% variables

omega = 100;
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

m1 = zeros(r); m5 = zeros(r,1);

for jdx = 1:r %j
   for kdx = 1:r %k
       m1(jdx,kdx) = cos(e(r)*(t(jdx)-t(kdx)));
    end
end

for jdx = 1:r %j
        m5(jdx) = cos(e(r)*(t(jdx)-tau));
end
%% 
% gradient method

x_1 = 0.1*ones(r,1); % start

MAX_ITERS = 3000; 

% with diminishing step sizes
as = [0.00001, 0.000001, 0.00005];
[c1,hist1,fval1] = gradient_descent_Emax(m1,m5,omega,x_1,as(1),MAX_ITERS); 
[c2,hist2,fval2] = gradient_descent_Emax(m1,m5,omega,x_1,as(2),MAX_ITERS); 
[c3,hist3,fval3] = gradient_descent_Emax(m1,m5,omega,x_1,as(3),MAX_ITERS); 

label1 = sprintf('%.3f/k',as(1));
label2 = sprintf('%.4f/k',as(2));
label3 = sprintf('%.5f/k',as(3));

% setup plot data
f1 = hist1{1};
f2 = hist2{1};
f3 = hist3{1};
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

a = [fval1, fval2, fval3]; 
[obj_omega, idx] = min(a)
c_ary = [c1 c2 c3]; c = c_ary(:,idx)
c_1norm = norm(c,1)
delE = (obj_omega - c_1norm)/omega