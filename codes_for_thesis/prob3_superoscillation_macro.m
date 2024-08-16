%% 
% variables

N = 19; % >> E*max(T,tau)
r = N+1;
E_max = 1;
T = 1;
for_parameter = 80;
s_values = [];
time_vectors = [];
epsilons = zeros(r,for_parameter);
s_vectors = zeros(r,for_parameter);
c_vectors = zeros(r,for_parameter);

for kdx = 1:for_parameter
    tau = 0.05*kdx;
    time_vectors(end+1) = tau;

    t = zeros(N+1,1);
    for idx = 1:N+1
        t(idx) = (idx-1)*T/N;
    end

    e = zeros(r,1);
    for idx = 1:r
        e(idx) = (-1 + 2*(idx-1)/(r-1))*E_max;
    end

    c= zeros(N+1,1);
    for jdx = 1:N+1 % j runs from 0 to N (N+1) terms
        c(jdx) = power((tau/T),N)*nchoosek(N,jdx-1)*power((T/tau-1),N-jdx+1);  
    end
    
    [M1,ind,Eset] = set_maximizers(r,tau,e,t,c);
    [gradient_max] = prob3_gradient(r,tau,e,t,c,M1,ind); % the coeffient 2 included
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

    s_values(end+1) = cvx_optval;
    epsilons(:,kdx) = epsilon';
    s_vectors(:,kdx) = s;
    c_vectors(:,kdx) = c;

end
%%
% combine result matrices
zero_row = ones(1,for_parameter);
result = [time_vectors; s_values; zero_row; s_vectors; zero_row; epsilons; zero_row; c_vectors]