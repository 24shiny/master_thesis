omegas = []; obj_omega = []; c_norm = []; delE = []; %c_gather = zeros(r,1);

for idx = 1 : 1
    omega_var = 10000;
    x_1 = zeros(r,1); % start

    MAX_ITERS = 1000; 

    % with diminishing step sizes
    as = [0.001, 0.0001, 0.00001];
    [c1,hist1] = gradient_descent(M1,m5,e,t,r,tau,omega_var,x_1,as(1),MAX_ITERS); 
    [c2,hist2] = gradient_descent(M1,m5,e,t,r,tau,omega_var,x_1,as(2),MAX_ITERS); 
    [c3,hist3] = gradient_descent(M1,m5,e,t,r,tau,omega_var,x_1,as(3),MAX_ITERS); 

    f1 = hist1{1};
    f2 = hist2{1};
    f3 = hist3{1};

    a = [f1(MAX_ITERS-1) f2(MAX_ITERS-1) f3(MAX_ITERS-1)]; 
    [obj_omega_temp, ind] = min(a);
    c_ary = [c1 c2 c3]; c = c_ary(:,ind);
    c_norm_temp = norm(c,1);
    delE_temp = (obj_omega_temp - c_norm_temp)/omega;

    omegas(end+1) = omega_var;
    obj_omega(end+1) = obj_omega_temp;
    c_norm(end+1) = c_norm_temp;
    delE(end+1) = delE_temp;
    %c_gather(:,end+1) = c';
end
%%
omegas = omegas';
obj_omega = obj_omega';
c_norm = c_norm';
delE = delE';
%%
result = [omegas, obj_omega, c_norm, delE];