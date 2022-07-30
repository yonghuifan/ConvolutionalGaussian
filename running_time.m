clear;close all;

%%% init
clear
close all;
% 
% yalmip('clear');
% rng(1);
% 
% base_path = 'G:\SHREC14\Real\Data\';
% 
% dire = dir(fullfile(base_path, '*.off'));
% 
% params.preprocess.dataFileType = 'off';
% numLmk = 500;
% 
% BNN = 100;
% 
% idx = 0;
% %for i = 1:length(dire)
% G = Mesh(params.preprocess.dataFileType,fullfile(base_path,[int2str(idx),'_new.off']));
% node = G.V;
% face = G.F;
% node = PerformMeshSmoothing(G,node);
% G.V = node';
% G.Normalize();
% [~,~,flip] = G.ComputeNormal();
% if flip
%     G.F = G.F([1 3 2],:);
% end
% [idx, dist] = knnsearch(G.V',G.V','K',BNN+1);
% computing_time2 = zeros(7,length(20:50:1000));
% count = 1;
% for numLmk = 20:50:1000
%     tic;
%     [idx, dist] = knnsearch(G.V',G.V','K',BNN+1);
%     Lmks = GetGPLmk_Fan(G, numLmk, BNN, idx, dist);
%     computing_time(1, count) = toc;
%     tic;
%     Lmks = GetGPLmk(G, numLmk);
%     computing_time(2, count) = toc;
%     tic;
%     Lmks = GetGPLmk_Euclidean(G, numLmk);
%     computing_time(3, count) = toc;
%     tic;
%     Lmks = GetGPLmk_NoCurvature(G, numLmk);
%     computing_time(4, count) = toc;
%     tic;
%     [idx, dist] = knnsearch(G.V',G.V','K',BNN+1);
%     Lmks = GetGPLmk1(G, numLmk, BNN, idx, dist);
%     computing_time(5, count) = toc;
%     tic;
%     Mesh.v = node;
%     Mesh.f = face;
%     [meshSaliency, az, el, az2, el2] = meshSaliencyPipeline(Mesh);
%     [B,Idx] = maxk(meshSaliency, numLmk);
%     computing_time(6, count) = toc;
%     tic;
%     [idx, dist] = knnsearch(G.V',G.V','K',BNN+1);
%     Lmks = GetGPLmk1(G, numLmk, BNN, idx, dist);
%     computing_time(7, count) = toc;
%     count = count + 1;
% end
width = 5.5;     % Width in inches
height = 5.5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 10;      % Fontsize
lw = 2;      % LineWidth
msz = 8;       % MarkerSize
computing_time = load('computing_time.mat');
%computing_time2 = load('computing_time.mat');
f=figure(1);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
%acc_cur_mark = log(acc_cur_mark);
acc_cur_mark = [computing_time.computing_time(:,1:20)];
%x = 20:50:1000;
x = 20:50:1000;
plot(x,acc_cur_mark(2,:)','Color','#D95319','LineWidth',2); % WGP
hold on;
plot(x,acc_cur_mark(3,:)','Color','#77AC30','LineWidth',2); % RBF
hold on;
plot(x,acc_cur_mark(5,:)','Color','#EDB120','LineWidth',2); % SMK
hold on;
plot(x,acc_cur_mark(6,:)','Color','#4DBEEE','LineWidth',2); % MS
hold on;
plot(x,acc_cur_mark(4,:)','Color','#0000FF','LineWidth',2); % matern
hold on;
plot(x,acc_cur_mark(1,:)','Color','#FF0000','LineWidth',2); % ours

% hold on;
% plot(acc_cur_mark(7,:)','Color','k','LineWidth',2); % PK-GP
%ylim([1,4.2]);
%x=xlabel('Number of salient points','FontSize',12);
%set(x,'position',get(x,'position')+[0 10 0]);
%t=title('\fontsize{16}(d) Average running time','FontName','Times New Roman');
%set(t,'position',get(t,'position')-[0 1110 0]);
xlabel('Number of salient points','FontSize',16,'FontWeight','bold');
ylabel('Running time/second','FontSize',16,'FontWeight','bold');
%legend('Weighted','HK-GP','SMK-GP','MeshSaliency','PK-GP','PPD-GP','Location','northwest','FontSize',12);
legend('W-GP','RBF-GP','SMK-GP','MeshSaliency','Matern-GP','Ours','Location','northwest','FontSize',16,'FontName','Times New Roman','FontWeight','bold');
set(gca,'XTick',0:200:1000,'FontSize',16,'FontWeight','bold');
set(gca,'YTick',0:200:1000,'FontSize',16,'FontWeight','bold');
%set(gca,'XTick',0:50:300,'FontSize',16,'FontWeight','bold');
%set(gca,'YTick',0:50:300,'FontSize',16,'FontWeight','bold');
%t=title('(d) Average running time','FontSize',14,'FontWeight','bold');
%set(t,'position',get(x,'position')-[0 85 0]);
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print(f,'f','-depsc2','-r800');
system('gswin64c -o -q -sDEVICE=png256 -dEPSCrop -r800 -o D:\YonghuiFan\NIPS2020\runnimg_time_eps.png C:\Users\yfan61\Documents\MATLAB\f.eps');