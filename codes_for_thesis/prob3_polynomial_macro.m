%% 
% variables

r = 10;
T = 1;
E_max = 1;
s_values = [];
for_parameter = 81;
epsilons = zeros(r,for_parameter);
s_vectors = zeros(r,for_parameter);
c_vectors = zeros(r,for_parameter);

for kdx = 1:for_parameter
    tau = 0+ 0.05*kdx;

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

    c = prod(c_temp,2);

    c_vectors(:,kdx) = c;

    [M1,ind,Eset] = set_maximizers(r,tau,e,t,c);
    % compute opt elements;
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
            epsilon >= 0.00001;
            Lambda >= 0;
    cvx_end

    s_values(end+1) = cvx_optval;
    epsilons(:,kdx) = epsilon';
    s_vectors(:,kdx) = s;
    c_vectors(:,kdx) = c;
end