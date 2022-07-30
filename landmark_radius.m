%%% init
clear
close all;
%initialize;
% yalmip('clear');
rng(1);

%base_path = 'C:\Users\yfan61\Documents\MATLAB\CPsurfcomp\DATA\teeth\meshes';
base_path = 'G:\CPsurfcomp\DATA\radius\meshes\';
%base_path = 'D:\YonghuiFan\Data\CPsurfcomp\DATA\radius\meshes';
% base_path = 'D:\YonghuiFan\ICCV2019\SHREC14\Real\Data';
% base_path = 'C:\Users\yfan61\Documents\MATLAB\GPLmkBDMatch-master\sample_data';
% base_path = 'D:\YonghuiFan\Data\ModelNet10\bathtub\train';

dire = dir(fullfile(base_path, '*.off'));

params.preprocess.dataFileType = 'off';
numLmk = 100;

BNN = 50;

% exe_dir = 'C:\Users\yfan61\Documents\MATLAB';
% meshfix_dir = [exe_dir,'\meshfix.exe'];
% FilterMesh_dir = [exe_dir,'\FilterMesh.exe'];
% MeshSimplify_dir = [exe_dir,'\MeshSimplify.exe'];

for i = 14:length(dire)% 13, 22, 23, 44, 34
    
    %if ~exist(outputfilename,'file')
    display(dire(i).name);
    G = Mesh(params.preprocess.dataFileType,fullfile(base_path,[dire(i).name]));
    %node = G.V;
    %node = PerformMeshSmoothing(G,node);
    %G.V = node';
    G.Normalize();
    [~,~,flip] = G.ComputeNormal();
    if flip
        G.F = G.F([1 3 2],:);
    end
    
    %dist = load(fullfile(base_path,[dire(i).name,'_geodesicDist.txt']));
    %idx = load(fullfile(base_path,[dire(i).name,'_geodesicIdx.txt']));
    %dist = dist(:,1:BNN+1);
    %idx = idx(:, 1:BNN+1);
    [idx, dist] = knnsearch(G.V',G.V','K',BNN+1);
    Lmks = GetGPLmk_Fan(G, numLmk, BNN, idx, dist);
    % Lmks = GetGPLmk1(G, numLmk, BNN, idx, dist);
    % Lmks = GetGPLmk(G, numLmk);
    % Lmks = GetGPLmk_Euclidean(G, numLmk);
    %Lmks = GetGPLmk_NoCurvature(G, numLmk);
    % Lmks = GetGPLmk_periodic(G, numLmk, BNN, idx, dist);
    Lmks = reshape(Lmks,1,[]);
    
    % [node,face] = readOFF(fullfile(base_path,[dire(i).name]));
    node = G.V';
    face = G.F';
     [BV,BE] = FindBoundaries(G);
     [bidx, ~] = knnsearch(BV',BV','K',31);
     BV = [BV,unique(bidx)'];
%     %[C, ia, ib] = intersect()
     tmp = setdiff(Lmks, sort(unique(BV)),'stable');
    %tmp(Lmk_needed+1:end) = [];
    %tmp = Lmks;
    %parsave([fullfile(base_path,[dire(i).name]),'_periodic_landmark1000.txt'],tmp);
    CMap=colormap(parula(10));
    % [Cgauss,Cmean,Umin,Umax,Cmin,Cmax,Normal] = ComputeCurvature(G);
    [GC, MC]= surfcurvatures(node(:,1),node(:,2),node(:,3),face);
    %MC = MC;
    GC = normalize(GC);
    % MC_min = min(MC);
    % MC_max = max(MC);
    % GC_min = min(GC);
    % GC_max = max(GC);
    node(:,3) = -node(:,3);
    figure(1);
    clf 
    %set(img, 'Position', [300 400 400 400]); 
    hold on
    axis equal
    %set(gca, 'Layer', 'top')
    pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',abs(GC*400),'FaceColor','interp','EdgeColor','none','BackFaceLighting','lit');
    %pp= patch('Faces',face,'Vertices',node,'FaceColor','interp','EdgeColor','none');
    %view(136, -13);
    view(92,-82);
    %ax = gca;
    %camroll(-45);
    %ax.CameraPosition = [0 10 0];
    caxis([0,400]);
    axis off;
    %colormap jet
    %colorbar
    hold on
    sphereSize = 0.018; 
    c = [1 0 0];
    
    drawSpheres(node(tmp,:),sphereSize,c);
%     color = ones(size(node,1),3);
%     color(:,1) = 227;
%     color(:,2) = 218;
%     color(:,3) = 201;
%     color(BV,1) = 255;
%     color(BV,2) = 0;
%     color(BV,3) = 0;
    %writeOFF('bone_vis.off', node, face, [], color);
    % end
end
