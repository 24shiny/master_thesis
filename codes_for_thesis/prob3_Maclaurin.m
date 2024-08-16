%% 
% variables

% set N(high) and h(low)
N = 5;
r = N+1;
h = 1/5;
tau = 0.7;
E_max = 1;
T = 2;

t = zeros(N+1,1);
for idx = 1:N+1
    t(idx) = ((idx-1)*h)*T;
end

t

e = zeros(r,1);
for idx = 1:r
    e(idx) = (-1 + 2*(idx-1)/(r-1))*E_max;
end

c_temp = zeros(N+1);
for jdx = 1:N+1 % j runs from 0 to N (N+1) terms
    for kdx = 1:jdx % k
        c_temp(jdx,kdx) = abs(tau^(jdx-1)/(factorial(jdx-1)*h^(kdx-1))*nchoosek(jdx-1,kdx-1)*(-1)^(kdx+jdx-2));
    end    
end

c = sum(c_temp,2)
%%
% a set of maximizers
[M1,ind,Eset] = set_maximizers(r,tau,e,t,c);
% compute opt elements;
[gradient_max] = prob3_gradient(r,tau,e,t,c,M1,ind); % the coeffient 2 included
%% 
% opt problem

cvx_begin
    variables s(r) epsilon(r) Lambda(length(ind));
    minimize(norm(s,1));
    constraint_obj = zeros(r,1);
    for idx = 1:length(ind)
        constraint_obj_temp = Lambda(idx)*gradient_max{idx};
        constraint_obj = constraint_obj + constraint_obj_temp;
    end
    subject to
        (epsilon(r)+s).*sign(c) + constraint_obj == 0;
        sum(Lambda) == 1;
        epsilon >= 0;
        Lambda >= 0;
cvx_end

cvx_optval
%%
% plot(epsilon);
% ylabel({'error ϵ_j'});
% xlabel({'time t_j'});
% grid on
% title_spec = sprintf('N = %.1f, h = %.1f',N, h);
% title(title_spec);