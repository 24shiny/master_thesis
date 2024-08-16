lambda = 
tau = 2;

r = 101; % to check e=0
t = zeros(r,1);

for idx = 1:r
    t(idx) = (idx-1)/(r-1);
end

e = zeros(r,1);
for idx = 1:r
    e(idx) = -1 + 2*(idx-1)/(r-1);
end

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
        M4{idx} = x'*M1{idx}*x; % change variable c -> x
end

for idx = 1:r
        M5{idx} = 2*m5(idx,:)*x;
end
maxobj = cellfun(@minus,M4,M5,'Un',0); %M4-2*M5+1
maxobj = cellfun(@(x) x+1,maxobj,'un',0);