function [D,V] = PPFeature(G,numLmk, BNN, idx, dist)
% [BV,~] = FindBoundaries(G);
% [bidx, ~] = knnsearch(BV,BV,'K',6);
% BV = [BV,bidx];
N = 5;
omega = sqrt((0.2 * 2* pi .*(1:N))/2);
% o = (0.2 * 2* pi .*(1:N));
%omega = sqrt(Lambda * 0.5);
alpha = linspace(0, 0.5*pi, N);
num = size(G,1);
fullPhi = sparse(G.nV,G.nV);
% for i = 1:size(dist, 1)
%     mink = min(dist(i, 2:end));
%     dist(i,2:end) = (dist(i,2:end)-mink)/max(dist(i, 2:end));
% end
normalized_dist = dist;
normalized_dist(:,1) = 0;
%normalized_dist = normalize(dist,2,'range');
%normalized_dist(:,2:end) = (normalized_dist(:,2:end) - min(min(normalized_dist(:,2:end))))./max(max(normalized_dist(:,2:end)));
%normalized_dist(:,2:end) = (normalized_dist(:,2:end))./max(max(normalized_dist(:,2:end)));
for oi = 1:N
    % fullPhi = fullPhi + sparse(repmat(1:G.nV,1,BNN+1),idx,(0.25*(1/pi)).*(1./normalized_dist).*cos(omega.*dist + Lambda * alpha(oi)).*exp(-omega.*dist),G.nV,G.nV);
    fullPhi = fullPhi + sparse(repmat(1:G.nV,1,BNN+1),idx,(0.25*(1/pi)).*(1./dist).*abs(cos(omega(oi).*dist + alpha(oi))).*exp(-omega(oi).*dist),G.nV,G.nV);
    %fullPhi = fullPhi + sparse(repmat(1:G.nV,1,BNN+1),idx,(0.25*(1/pi)).*cos(omega(oi).*dist+ alpha(oi)).*exp(-omega(oi).*dist),G.nV,G.nV);
    %fullPhi = fullPhi + sparse(repmat(1:G.nV,1,BNN+1),idx, exp(cos(omega(oi).*dist)),G.nV,G.nV);
    %fullPhi = fullPhi + sparse(repmat(1:G.nV,1,BNN+1),idx,0.2.*exp(-2*sqrt(5).*dist)+(1+2*sqrt(5).*dist+(1/0.15).*(dist.^2)),G.nV,G.nV);
    % fullPhi = (fullPhi+fullPhi')/2;
    fullPhi(1:(G.nV+1):G.nV*G.nV) = 0;
end
Lambda = sum(fullPhi,2)/BNN;
%fullPhi = 0.2 .* fullPhi;
%Lambda = (0.25*sqrt(pi)).*(1./mean(dist(:,2:end).^2,2));
fullPhi = (fullPhi+fullPhi')/2;
%K_full = fullPhi * sparse(1:dim,1:dim,Lambda,dim,dim) * fullPhi;
%fullMatProd = fullPhi * diag(Lambda) * fullPhi;
%fullMatProd = fullPhi * fullPhi;
%fullMatProd = fullPhi;
%fullMatProd = fullPhi + diag(Lambda);
% PDistMat = squareform(pdist(G.V'));
% fullPhi = exp(-PDistMat.^2/bandwidth);

%disp('Constructing full kernel......');
%tic;
fullMatProd = fullPhi * sparse(1:G.nV,1:G.nV,Lambda,G.nV,G.nV) * fullPhi;
%fullMatProd = fullPhi * fullPhi;
%disp(['full kernel constructed in ' num2str(toc) ' sec.']);
[V, D] = eigs(fullMatProd,51,'sm');
V = V(:,2:end);
D = D(2:end, 2:end);
D = diag(D);
D=abs(real(D));

end

