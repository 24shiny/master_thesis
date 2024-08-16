%% A. HS unconstrained
% variables

omega = 10000;
tau = 2;
r = 101; % to check e=0
t = zeros(r,1);
d = -0.0001 + (0.0002)*rand(1,r) % random perturbation

for idx = 1:r
    t(idx) = (idx-1)/(r-1);
end

e = zeros(r,1);
for idx = 1:r
    e(idx) = -1 + 2*(idx-1)/(r-1);
end
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

    minimize(norm(c,1)+omega*max([maxobj{:}]));  %+d*c
cvx_end
%%
scatter(1:r,c)
xlabel('j')
ylabel('c_j')
label1 = sprintf('omega = %d, tau=%.1f',omega,tau);
title(label1);
cvx_optval
norm(c,1)
del = sqrt((cvx_optval - norm(c,1))/omega)
%%
% trials and errors
%m3 = c*c'; % trial1 not accepted
%m3 = c'*M1{idx}*c; % trail2 cannot put cvx objs into an array
%m4(idx) = m3;
%maxobjmat = cell2mat(maxobj); % trial3 cell2mat does not support cells containing objects
%cellfun(@(x) max(cell2mat(x)),A)
%B = [A{:}] % another way to transform cell to array