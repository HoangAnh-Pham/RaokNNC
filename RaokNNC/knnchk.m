%% Find k nearest neighbors
function [n,D] = knnchk(X,P,k)
    NP = length(P(:,1));
    Bound = max(P)*(1+0.000001)-min(P);
    for i=1:NP
        Y = P(i,:);
        dis(i) = (sum(((X-Y)./Bound).^2))^0.5;
    end
    [D,id] = sort(dis); 
    n = id(1:k);
    D = D(1:k);
end
