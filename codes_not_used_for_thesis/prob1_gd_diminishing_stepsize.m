%% check values of E
% define variables

delta = 10;
tau = 2;
r = 101; % not 100 to check e=0
t = zeros(r,1);
for idx = 1:r
    t(idx) = (idx-1)/(r-1);
end

vals = zeros(r,1);
% For loop for chaning e and cvx

for idx = 1:r
    e = -1 + 2*(idx-1)/(r-1); % runs from -1 through 0 to 1

    expt1 = zeros(r,1);
    for jdx = 1:r
        expt1(jdx) = exp(-1*1i*t(jdx)*e);
    end

    expt2 = exp(-1i*tau*e);
    
    cvx_begin
        variable c(r);
        minimize(norm(c,1)+delta*abs(sum(expt1.*c) - expt2));
    cvx_end

    vals(idx) = cvx_optval;
end