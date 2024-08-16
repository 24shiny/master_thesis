r = 10;

t = zeros(r,1);
for idx = 1:r
    t(idx) = (idx-1)/(r-1);
end

e = zeros(r,1);
for idx = 1:r
    e(idx) = -1 + 2*(idx-1)/(r-1);
end

test_mat = zeros(r);
collect_eig = zeros(r);
for kdx = 1:r
    for idx = 1:r
        for jdx = 1:r
            test_mat(idx,jdx) = cos(e(kdx)*(t(idx)-t(jdx)));   
        end
    end
    collect_eig(:,kdx) = eig(test_mat);
end