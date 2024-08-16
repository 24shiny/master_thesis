function [delE_max] = delE_maximizers(r,tau,e,t,c,ind)

expt1 = zeros(r);
for idx = 1:r
    for jdx = 1:r
        expt1(idx, jdx) = exp(-1*1i*t(jdx)*e(idx));
    end
end

expt2 = exp(-1i*tau*e);

delE_max = [];
for idx = 1:length(ind)
    delE_temp = abs(expt1(idx,:)*c - expt2(ind(idx)));
    delE_max(end+1) = delE_temp;
end