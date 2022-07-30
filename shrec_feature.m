%%% init
clear
close all;
%initialize;
yalmip('clear');
%rng(1);

base_path = 'G:\SHREC14\Real\Data\';
% base_path = 'D:\YonghuiFan\ICCV2019\SHREC18_protein';
% base_path = 'D:\YonghuiFan\Data\ModelNet10\bathtub\train';

params.preprocess.dataFileType = 'off';
numLmk = 120;
BNN = 800;

% exe_dir = 'C:\Users\yfan61\Documents\MATLAB';
% meshfix_dir = [exe_dir,'\meshfix.exe'];
% FilterMesh_dir = [exe_dir,'\FilterMesh.exe'];
% MeshSimplify_dir = [exe_dir,'\MeshSimplify.exe'];

for i = 6:400
    %if ~exist([base_path,int2str(i-1),'_fan_landmark1000.txt'],'file')
    filename = fullfile(base_path, [int2str(i-1),'_new.off']);
    display(filename);
    G = Mesh(params.preprocess.dataFileType,filename);
    G.Normalize();
    % reorient (outward facing normals)
    [~,~,flip] = G.ComputeNormal();
    if flip
        G.F = G.F([1 3 2],:);
    end
    [node,face] = readoff(fullfile(base_path,[int2str(i-1),'_new.off']));
    pp= patch('Faces',face,'Vertices',node);
    node = G.V';
    face = G.F';
    %dist = load(fullfile(base_path,[int2str(i-1),'_geodesicDist.txt']));
    %idx = load(fullfile(base_path,[int2str(i-1),'_geodesicIdx.txt']));
    [idx, dist] = knnsearch(G.V',G.V','K',BNN+1);
    %dist = dist(:,1:BNN+1);
    %idx = idx(:,1:BNN+1);
    [eigvalue,eigvect] = PPFeature(G, numLmk, BNN, idx, dist);
    %Lmks = GetGPLmk(G, numLmk);    
    Lmks = reshape(Lmks,1,[]);
    
    %[BV,BE] = FindBoundaries(G);
    %tmp = setdiff(sort(Lmks), sort(BV));
    tmp = Lmks;
    %parsave([base_path,int2str(i-1),'_fan3_landmark1000.txt'],'tmp');
    CMap=colormap(parula(10));
    [Cgauss,Cmean,Umin,Umax,Cmin,Cmax,Normal] = ComputeCurvature(G);
    %Cmean = normalize(Cmean,'range',[-80 80]);
    %Cgauss = normalize(Cgauss,'range',[0,100]);
    %[GC, MC]= surfcurvatures(node(:,1),node(:,2),node(:,3),face);
    %MC = MC;
    
    % MC_min = min(MC);
    % MC_max = max(MC);
    % GC_min = min(GC);
    % GC_max = max(GC);
    node(:,3) = -node(:,3);
    img=figure(2);
    clf 
    set(img, 'Position', [300 100 600 600]); 
    hold on
    %axis equal
    set(gca, 'Layer', 'top')
    pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',abs(normalize(Cgauss))*150,'FaceColor','interp','EdgeColor','none');
    %pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',Cmean,'FaceColor','interp','EdgeColor','none');
    %view(90,180);
    view(180, -90);
    %camroll(90)
    ax = gca;
    %ax.CameraPosition = [0 10 0];
    caxis([0,400]);
    axis off;
%     %colormap jet
%     %colorbar
    hold on
    sphereSize = 0.008; 
    c = [1 0 0];
    drawSpheres(node(tmp,:),sphereSize,c);
    % drawSpheres(node(Lmks{j},:),sphereSize,[1 0 1]);
    %end
end
