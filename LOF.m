function lof = LOF(X,k)
% X is a matrix with the data points
% k is the number of nearest neighbors to use

    [n,p] = size(X);
    D = pdist2(X,X);
    kDist = zeros(n,n);
    for i = 1:n
        kDist(i,:) = sort(D(i,:), 'MissingPlacement', 'last');
    end
    kDist = kDist(:,k);
    lrd = zeros(n,1);
    for i = 1:n
        idx = D(i,:)<=kDist(i);
        rd = kDist(idx);
        lrd(i) = 1/(sum(rd)/k);
    end
    lrd(isnan(lrd))=0;
    lof = zeros(n,1);
    for i = 1:n
        idx = D(i,:)<=kDist(i);
        lof(i) = sum(lrd(idx))/lrd(i)/k;
    end
    lof = lof';
end
