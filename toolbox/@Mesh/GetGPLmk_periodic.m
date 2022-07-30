function [GPLmkIdx,ptuq] = GetGPLmk_periodic(G,numLmk, BNN, idx, dist)

N = 5;
omega = sqrt((0.2 * 2* pi .*(1:N))/2);
% o = (0.2 * 2* pi .*(1:N));
%omega = sqrt(Lambda * 0.5);
alpha = linspace(0, 0.5*pi, N);

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
    % fullPhi = fullPhi + sparse(repmat(1:num,1,BNN+1),idx,(0.25*(1/pi)).*(1./dist).*abs(cos(omega(oi).*dist + alpha(oi))).*exp(-omega(oi).*dist),num,num);
    %fullPhi = fullPhi + sparse(repmat(1:G.nV,1,BNN+1),idx,(0.25*(1/pi)).*cos(omega(oi).*dist+ alpha(oi)).*exp(-omega(oi).*dist),G.nV,G.nV);
    fullPhi = fullPhi + sparse(repmat(1:G.nV,1,BNN+1),idx, exp(cos(omega(oi).*dist)),G.nV,G.nV);
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

KernelTrace = diag(fullMatProd);
% GPLmkIdx = zeros(1,numLmk);
GPLmkIdx = [];
true_GPLmkIdx = [];
invKn = zeros(numLmk);

cback = 0;

j = 1;
while(length(GPLmkIdx)<numLmk)    
%for j=1:numLmk
%j = length(GPLmkIdx)+1;
    for cc=1:cback
        fprintf('\b');
    end
    cback = fprintf('Landmark: %4d\n',j);
    
    if j == 1
        ptuq = KernelTrace;
    else
        if j == 2
            invKn(1:(j-1),1:(j-1)) = 1/fullMatProd(GPLmkIdx(1),GPLmkIdx(1));
            ptuq = KernelTrace - sum(fullMatProd(:,GPLmkIdx(1:(j-1)))'...
                .*(invKn(1:(j-1),1:(j-1))*fullMatProd(GPLmkIdx(1:(j-1)),:)),1)';
        else
            p = fullMatProd(GPLmkIdx(1:(j-2)),GPLmkIdx(j-1));
            mu = 1./(fullMatProd(GPLmkIdx(j-1),GPLmkIdx(j-1))-p'*invKn(1:(j-2),1:(j-2))*p);
            invKn(1:(j-2),1:(j-1)) = invKn(1:(j-2),1:(j-2))*[eye(j-2)+mu*(p*p')*invKn(1:(j-2),1:(j-2)),-mu*p];
            invKn(j-1,1:(j-1)) = [invKn(1:(j-2),j-1)',mu];
            productEntity = invKn(1:(j-1),1:(j-1))*fullMatProd(GPLmkIdx(1:(j-1)),:);
            ptuq = KernelTrace - sum(fullMatProd(:,GPLmkIdx(1:(j-1)))'...
                .*productEntity,1)';
        end
    end
    [~,maxUQIdx] = max(ptuq);
%     ptuq =  (ptuq - min(ptuq))/max(ptuq)-10;
%     figure(j);
%     pp= patch('Faces',G.F','Vertices',G.V','FaceVertexCData',full(ptuq)*100000,'FaceColor','interp','EdgeColor','none');
%     axis off;
%     view(103, -90);
%     camroll(-180);
    GPLmkIdx = [GPLmkIdx, maxUQIdx]; 
%     if ismember(maxUQIdx, BV)
%         true_GPLmkIdx = [true_GPLmkIdx, 0];
%     end
     j = j + 1;
end

p = fullMatProd(GPLmkIdx(1:(end-1)),GPLmkIdx(end));
mu = 1./(fullMatProd(GPLmkIdx(end),GPLmkIdx(end))-p'*invKn(1:(end-1),1:(end-1))*p);
invKn(1:(end-1),:) = invKn(1:(end-1),1:(end-1))*[eye(numLmk-1)+mu*(p*p')*invKn(1:(end-1),1:(end-1)),-mu*p];
invKn(end,:) = [invKn(1:(end-1),end)',mu];
ptuq = KernelTrace - sum(fullMatProd(:,GPLmkIdx)'...
    .*(invKn*fullMatProd(GPLmkIdx,:)),1)';

end

