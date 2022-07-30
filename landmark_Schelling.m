%%% init
clear
close all;
initialize;
yalmip('clear');
rng(1);

base_path = 'G:\SchellingData\Meshes';

params.preprocess.dataFileType = 'off';
numLmk = 200;
BNN =500;

for i =255:400
    filename = fullfile(base_path, [int2str(i),'.off']);
    G = Mesh(params.preprocess.dataFileType,filename);
     G.Normalize();
%     % reorient (outward facing normals)
     [~,~,flip] = G.ComputeNormal();
     if flip
         G.F = G.F([1 3 2],:);
     end
    %dist = load(fullfile(base_path,[dire(i).name,'_geodesicDist.txt']));
    %idx = load(fullfile(base_path,[dire(i).name,'_geodesicIdx.txt']));
    [idx, dist] = knnsearch(G.V',G.V','K',BNN+1);
    Lmks = GetGPLmk_Fan(G, numLmk, BNN, idx, dist);
    %Lmks = GetGPLmk(G, numLmk);
    % Lmks = GetGPLmk_Euclidean(G, numLmk);
    %Lmks = GetGPLmk_NoCurvature(G, numLmk);    
    Lmks = reshape(Lmks,1,[]);
    % [node,face] = readOFF(fullfile(base_path,[dire(i).name]));
    node = G.V';
    face = G.F';
    %[BV,BE] = FindBoundaries(G);
    %tmp = setdiff(sort(Lmks), sort(BV));
    tmp = Lmks;
    %save(fullfile(base_path, [int2str(i),'_Fan_landmark500.txt']),'tmp','-ascii');
%     CMap=colormap(parula(10));
%     [Cgauss,Cmean,Umin,Umax,Cmin,Cmax,Normal] = ComputeCurvature(G);
%     %[GC, MC]= surfcurvatures(node(:,1),node(:,2),node(:,3),face);
%     %MC = MC;
%     
%     % MC_min = min(MC);
%     % MC_max = max(MC);
%     % GC_min = min(GC);
%     % GC_max = max(GC);
%     %node(:,3) = -node(:,3);
%     img=figure(1);
%     clf 
%     set(img, 'Position', [300 100 600 600]); 
%     hold on
%     %axis equal
%     set(gca, 'Layer', 'top')
%     pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',Cmean,'FaceColor','interp','EdgeColor','none');
%     %pp= patch('Faces',face,'Vertices',node,'FaceVertexCData',Cmean,'FaceColor','interp','EdgeColor','none');
%     %view(90,180);
%     ax = gca;
%     %ax.CameraPosition = [0 10 0];
%     caxis([0,400]);
%     axis off;
%     %colormap jet
%     %colorbar
%     hold on
%     sphereSize = 0.008; 
%     c = [1 0 0];
%     drawSpheres(node(tmp,:),sphereSize,c);
    % drawSpheres(node(Lmks{j},:),sphereSize,[1 0 1]);
end
