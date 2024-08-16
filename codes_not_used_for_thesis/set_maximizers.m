function [M1,ind,Eset] = set_maximizers(r,tau,e,t,c)

m1 = zeros(r); M1 = cell(r,1);
m5 = zeros(r); 

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

M4 = cell(r,1); M5 = cell(r,1);

for idx = 1:r
        M4{idx} = c'*M1{idx}*c; % change variable c -> x
end

for idx = 1:r
        M5{idx} = 2*m5(idx,:)*c;
end
maxobj1 = cellfun(@minus,M4,M5,'Un',0); %M4-2*M5+1
maxobj1 = cellfun(@(x) x+1,maxobj1,'un',0);

max_val = max([maxobj1{:}]);
[nu,ind] = find([maxobj1{:}]==max_val);
Eset = e(ind);

%[M,I] = max(A)
% maximum = max(max(A));
% [x,y]=find(A==maximum)
end