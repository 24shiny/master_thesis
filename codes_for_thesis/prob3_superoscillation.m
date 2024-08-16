%% 
% variables

N = 19; % >> E*max(T,tau)
r = N+1;
tau = 0.7;
E_max = 1;
T = 1;

t = zeros(N+1,1);
for idx = 1:N+1
    t(idx) = (idx-1)*T/N;
end

t

e = zeros(r,1);
for idx = 1:r
    e(idx) = (-1 + 2*(idx-1)/(r-1))*E_max;
end

c= zeros(N+1,1);
for jdx = 1:N+1 % j runs from 0 to N (N+1) terms
    c(jdx) = power((tau/T),N)*nchoosek(N,jdx-1)*power((T/tau-1),N-jdx+1);  
end
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