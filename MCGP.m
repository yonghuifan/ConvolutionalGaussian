clear;close all;

base_path = 'G:\SHREC14\Real\Data\';
% base_path = 'D:\YonghuiFan\ICCV2019\SHREC18_protein';
% base_path = 'D:\YonghuiFan\Data\ModelNet10\bathtub\train';

params.preprocess.dataFileType = 'off';
numLmk = 500;
BNN = 50;

% exe_dir = 'C:\Users\yfan61\Documents\MATLAB';
% meshfix_dir = [exe_dir,'\meshfix.exe'];
% FilterMesh_dir = [exe_dir,'\FilterMesh.exe'];
% MeshSimplify_dir = [exe_dir,'\MeshSimplify.exe'];

% Now we are doing 1000 landmark with 1 convolution operation
num_lmk = 250;

for i = 1:400
    %if ~exist([base_path,int2str(i-1),'_periodic_landmark.txt'],'file')
    lmk = load(fullfile(base_path, [int2str(i-1),'_fan_landmark.txt']));
    wks = load(fullfile(base_path, [int2str(i-1),'_WKS.txt']));
    wks_dist = load(fullfile(base_path,[int2str(i-1),'_wks_dist.mat']));
    wks_dist = wks_dist.dist;
    
    lmk = lmk(1:num_lmk);
    
    wks_dist_lmk = wks_dist(lmk,:);
    
    
    %dist = load(fullfile(base_path,[int2str(i-1),'_geodesicDist.txt']));
    idx = load(fullfile(base_path,[int2str(i-1),'_geodesicIdx.txt']));
    idx_lmk = idx(lmk,:);
    wks_lmk = wks(idx,:);
     filename = fullfile(base_path, [int2str(i-1),'_new.off']);
%     display(filename);
     G = Mesh(params.preprocess.dataFileType,filename);
     G.Normalize();
%     % reorient (outward facing normals)
     [~,~,flip] = G.ComputeNormal();
     if flip
         G.F = G.F([1 3 2],:);
     end
    %[idx, dist] = knnsearch(G.V',G.V','K',BNN+1);
    %dist = dist(:,1:BNN+1);
    idx = idx(:,1:BNN+1);
    wks_dist = wks_dist(:,1:BNN+1);
    %Lmks = GetGPLmk1(G, numLmk, BNN, idx, dist);
    
    %[BV,~] = FindBoundaries(G);
    %[~,curvature] = findPointNormals(G.V',10);
    %Lambda = G.Aux.VertArea.*curvature/sum(curvature);
%     Lambda = 0.5;
%     if size(G.E,1) > 2
%         [I,J] = find(tril(G.E));
%         G.E = ([I,J])';
%     end
%     EdgeIdxI = G.E(1,:);
%     EdgeIdxJ = G.E(2,:);
    
    N = 5;
    omega = sqrt((0.2 * 2* pi .*(1:N))/2);
    % o = (0.2 * 2* pi .*(1:N));
    %omega = sqrt(Lambda * 0.5);
    alpha = linspace(0, 0.5*pi, N);
    fullPhi = sparse(G.nV,G.nV);
    
    for oi = 1:N
        fullPhi = fullPhi + sparse(repmat(1:G.nV,1,BNN+1),idx,(0.25*(1/pi)).*(1./wks_dist).*abs(cos(omega(oi).*wks_dist + alpha(oi))).*exp(-omega(oi).*wks_dist),G.nV,G.nV);
        fullPhi(1:(G.nV+1):G.nV*G.nV) = 0;
    end
    Lambda = sum(fullPhi,2)/BNN;
    fullPhi = (fullPhi+fullPhi')/2;
    
    fullMatProd = fullPhi * sparse(1:G.nV,1:G.nV,Lambda,G.nV,G.nV) * fullPhi;
    
    conv1_feature = fullMatProd * wks;
    conv2_feature = fullMatProd * conv1_feature;
    
    conv2_feature_lmk = reshape(pca(conv2_feature(lmk,:)),1,[]);
    parsave(fullfile(base_path, [int2str(i-1),'_250_conv2.mat']),'conv2_feature_lmk');
    
%     KernelTrace = diag(fullMatProd);
%    
%     % GPLmkIdx = zeros(1,numLmk);
%     GPLmkIdx = [];
%     true_GPLmkIdx = [];
%     invKn = zeros(numLmk);
%     
%     cback = 0;
%     
%     j = 1;
%     while(length(GPLmkIdx)<numLmk)
%         %for j=1:numLmk
%         %j = length(GPLmkIdx)+1;
%         for cc=1:cback
%             fprintf('\b');
%         end
%         cback = fprintf('Landmark: %4d\n',j);
%         
%         if j == 1
%             ptuq = KernelTrace;
%         else
%             if j == 2
%                 invKn(1:(j-1),1:(j-1)) = 1/fullMatProd(GPLmkIdx(1),GPLmkIdx(1));
%                 ptuq = KernelTrace - sum(fullMatProd(:,GPLmkIdx(1:(j-1)))'...
%                     .*(invKn(1:(j-1),1:(j-1))*fullMatProd(GPLmkIdx(1:(j-1)),:)),1)';
%             else
%                 p = fullMatProd(GPLmkIdx(1:(j-2)),GPLmkIdx(j-1));
%                 mu = 1./(fullMatProd(GPLmkIdx(j-1),GPLmkIdx(j-1))-p'*invKn(1:(j-2),1:(j-2))*p);
%                 invKn(1:(j-2),1:(j-1)) = invKn(1:(j-2),1:(j-2))*[eye(j-2)+mu*(p*p')*invKn(1:(j-2),1:(j-2)),-mu*p];
%                 invKn(j-1,1:(j-1)) = [invKn(1:(j-2),j-1)',mu];
%                 productEntity = invKn(1:(j-1),1:(j-1))*fullMatProd(GPLmkIdx(1:(j-1)),:);
%                 ptuq = KernelTrace - sum(fullMatProd(:,GPLmkIdx(1:(j-1)))'...
%                     .*productEntity,1)';
%             end
%         end
%         [~,maxUQIdx] = max(ptuq);
%         %     ptuq =  (ptuq - min(ptuq))/max(ptuq)-10;
%         %     figure(j);
%         %     pp= patch('Faces',G.F','Vertices',G.V','FaceVertexCData',full(ptuq)*100000,'FaceColor','interp','EdgeColor','none');
%         %     axis off;
%         %     view(103, -90);
%         %     camroll(-180);
%         GPLmkIdx = [GPLmkIdx, maxUQIdx];
%         if ismember(maxUQIdx, BV)
%             true_GPLmkIdx = [true_GPLmkIdx, 0];
%         end
%         j = j + 1;
%     end
%     
%     p = fullMatProd(GPLmkIdx(1:(end-1)),GPLmkIdx(end));
%     mu = 1./(fullMatProd(GPLmkIdx(end),GPLmkIdx(end))-p'*invKn(1:(end-1),1:(end-1))*p);
%     invKn(1:(end-1),:) = invKn(1:(end-1),1:(end-1))*[eye(numLmk-1)+mu*(p*p')*invKn(1:(end-1),1:(end-1)),-mu*p];
%     invKn(end,:) = [invKn(1:(end-1),end)',mu];
%     ptuq = KernelTrace - sum(fullMatProd(:,GPLmkIdx)'...
%         .*(invKn*fullMatProd(GPLmkIdx,:)),1)';
    
    
    %end
end
