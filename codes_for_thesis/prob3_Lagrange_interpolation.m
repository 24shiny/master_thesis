%% 
% variables

r = 101; 
T = 1;
tau = 2;
E_max = 1;

e = zeros(r,1);
for idx = 1:r
    e(idx) = (-1 + 2*(idx-1)/(r-1))*E_max;
end

t = zeros(r,1);
for idx = 1:r
    t(idx) = ((idx-1)/(r-1))*T;
end

c_temp = zeros(r);
for idx = 1:r % assign values to c
    for jdx = 1:r
        c_temp(idx,jdx) = (tau-t(jdx))/(t(idx)-t(jdx));
    end
end

for idx = 1:r
    c_temp(idx,idx) = 1;
end

c = prod(c_temp,2)*power(10,-70); % too large
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
        epsilon >= 0.00001;
        Lambda >= 0;
cvx_end
%% 
% plot

% plot(e);
% ylabel({'error ϵ_j'});
% xlabel({'time t_j'});
% grid on
% title_spec = sprintf('T = %.1f, tau = %.1f, Emax = %d',T, tau, E_max);
% title(title_spec);