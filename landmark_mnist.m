clear;close all;
dataset = convertMNIST( 'G:/train-images.idx3-ubyte', 'G:/train-labels.idx1-ubyte');
image = dataset.images(:,:,1);
padding = zeros(30, 30);
padding(2:29,2:29) = image;
image = padding;
imshow(padding);
N = 5;
numLmk = 10;
omega = sqrt((0.2 * 2* pi .*(1:N))/2);
alpha = linspace(0, 0.5*pi, N);
K = zeros(30, 30);
for i = 2:29
    for j = 2:29
        c=image(i, j);
        up = image(i-1,j);
        down = image(i+1, j);
        left = image(i,j-1);
        right = image(i, j+1);
        for oi = 1:N
            dist = [norm(c-up), norm(c-down),norm(c-left),norm(c-right)];
            K(i, j) = K(i, j) + sum((0.25*(1/pi)).*abs(cos(omega(oi).*dist + alpha(oi))).*exp(-omega(oi).*dist));
        end
    end
end
K = max(max(K))-K(2:29, 2:29);
Lambda = sum(K,2);
%fullPhi = 0.2 .* fullPhi;
%Lambda = (0.25*sqrt(pi)).*(1./mean(dist(:,2:end).^2,2));
K = (K+K')/2;
%K_full = fullPhi * sparse(1:dim,1:dim,Lambda,dim,dim) * fullPhi;
%fullMatProd = fullPhi * diag(Lambda) * fullPhi;
%fullMatProd = fullPhi * fullPhi;
%fullMatProd = fullPhi;
%fullMatProd = fullPhi + diag(Lambda);
% PDistMat = squareform(pdist(G.V'));
% fullPhi = exp(-PDistMat.^2/bandwidth);

%disp('Constructing full kernel......');
%tic;
fullMatProd = K * diag(Lambda)* K;
%fullMatProd = fullPhi * fullPhi;
%disp(['full kernel constructed in ' num2str(toc) ' sec.']);

KernelTrace = fullMatProd;
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
    [row, col] = find(ismember(ptuq, max(ptuq(:))));
    maxUQIdx(j,:) = [row, col];
%     ptuq =  (ptuq - min(ptuq))/max(ptuq)-10;
%       figure(j);
%       pp= patch('Faces',G.F','Vertices',G.V','FaceVertexCData',full(ptuq)*100000,'FaceColor','interp','EdgeColor','none');
%       axis off;
%      view(103, -90);
% %     camroll(-180);
%     GPLmkIdx = [GPLmkIdx, maxUQIdx]; 
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


