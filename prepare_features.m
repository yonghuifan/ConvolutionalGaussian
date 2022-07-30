clear;close all;

base_path = 'G:\SHREC14\Real\Data\';
dire = dir([base_path, '*.off']);
numLmk = 1000;
data_hks = [];
data_sihks = [];
data_wks = [];
count = 1;
num_points = 300;%30  50  80  120  200 300
full_feature = [];
for class = 1: 10
    %     fan_sihks_tmp = [];
    fan_wks_tmp = [];
    %     gao_sihks_tmp = [];
    %     gao_wks_tmp = [];
    %     ec_sihks_tmp = [];
    %     ec_wks_tmp = [];
    %     nc_sihks_tmp = [];
    %     nc_wks_tmp = [];
    %     sm_sihks_tmp = [];
    %     sm_wks_tmp = [];
    %     ms_sihks_tmp = [];
    %     ms_wks_tmp = [];
    %     sihks_tmp = [];
    %     wks_tmp = [];
    %     p_sihks_tmp = [];
    %     p_wks_tmp = [];
    for i = (class-1):10:399
        fan_lmk = load(fullfile(base_path,[int2str(i),'_Fan_landmark400.txt']));% 1
        %         gao_lmk = load(fullfile(base_path,[int2str(i),'_Gao_landmark.txt']));% 2
        %         Euclidean_lmk = load(fullfile(base_path,[int2str(i),'_Euclidean_landmark.txt']));% 3
        %         NCurvature_lmk = load(fullfile(base_path,[int2str(i),'_NoCurvature_landmark.txt']));% 4
        %         SM_lmk = load(fullfile(base_path,[int2str(i),'_SM_landmark.txt']));% 5
        %         %MaxGaussian_lmk = load(fullfile(base_path,[int2str(i),'_MaxGaussian_landmark1000.txt']));% 5
        %         %MaxMean_lmk = load(fullfile(base_path,[int2str(i),'_MaxMean_landmark.txt']));% 6
        %         p_lmk = load(fullfile(base_path,[int2str(i),'_periodic_landmark.txt']));% 6
        %         MeshSaliency_lmk = load(fullfile(base_path,[int2str(i),'_MeshSaliency_landmark1000.txt']));% 7
        
        fan_lmk = fan_lmk(1:num_points);
        %         gao_lmk = gao_lmk(1:num_points);
        %         Euclidean_lmk = Euclidean_lmk(1:num_points);
        %         NCurvature_lmk = NCurvature_lmk(1:num_points);
        %         SM_lmk = SM_lmk(1:num_points);
        %         MeshSaliency_lmk = MeshSaliency_lmk(1:num_points);
        %         p_lmk = p_lmk(1:num_points);
        
        %data_hks = [data_hks; load(fullfile(base_path,[int2str(i),'_hks.txt']))];
        %         data_sihks = load(fullfile(base_path,[int2str(i),'_sihks.txt']));
        data_wks = load(fullfile(base_path,[int2str(i),'_WKS.txt']));
        % acc_cur_mark = zeros(7,length(1:10:500));
        % fan_sihks_tmp = [fan_sihks_tmp; reshape(pca(data_sihks(fan_lmk,:)),1,[])];
        fan_wks_tmp = [fan_wks_tmp; reshape(pca(data_wks(fan_lmk,:)),1,[])];%,'NumComponents',4
        
        % fan_wks_tmp = [fan_wks_tmp, pca(data_wks(fan_lmk,:))];
        % fan_wks_tmp = cat(3, fan_wks_tmp, pca(data_wks(fan_lmk,:)));
        %         gao_sihks_tmp = [gao_sihks_tmp; reshape(pca(data_sihks(gao_lmk,:)),1,[])];
        %         gao_wks_tmp = cat(3, gao_wks_tmp, pca(data_wks(gao_lmk,:)));
        %         %
        %         ec_sihks_tmp = [ec_sihks_tmp; reshape(pca(data_sihks(Euclidean_lmk,:)),1,[])];
        %         ec_wks_tmp = [ec_wks_tmp; reshape(pca(data_wks(Euclidean_lmk,:)),1,[])];
        %         %
        %         nc_sihks_tmp = [nc_sihks_tmp; reshape(pca(data_sihks(NCurvature_lmk,:)),1,[])];
        %         nc_wks_tmp = [nc_wks_tmp; reshape(pca(data_wks(NCurvature_lmk,:)),1,[])];
        %         %
        %         sm_sihks_tmp = [sm_sihks_tmp; reshape(pca(data_sihks(SM_lmk,:)),1,[])];
        %         sm_wks_tmp = [sm_wks_tmp; reshape(pca(data_wks(SM_lmk,:)),1,[])];
        %
        %         ms_sihks_tmp = [ms_sihks_tmp; reshape(pca(data_sihks(MeshSaliency_lmk,:)),1,[])];
        %         ms_wks_tmp = [ms_wks_tmp; reshape(pca(data_wks(MeshSaliency_lmk,:)),1,[])];
        %         %
        %         sihks_tmp = [sihks_tmp; reshape(pca(data_sihks),1,[])];
        %         wks_tmp = [wks_tmp; reshape(pca(data_wks),1,[])];
        %         p_sihks_tmp = [p_sihks_tmp; reshape(pca(data_sihks(p_lmk,:)),1,[])];
        %         p_wks_tmp = [p_wks_tmp; reshape(pca(data_wks(p_lmk,:)),1,[])];
    end
    % fan_category{class}.sihks = fan_sihks_tmp;
    % fan_category{class}.wks = fan_wks_tmp;
    % save(fullfile(base_path,['gao_category_',int2str(class),'_',int2str(num_points),'.mat']),'gao_wks_tmp');
    %save(fullfile(base_path,['fan_',int2str(class),'_',int2str(num_points),'.mat']),'fan_wks_tmp');
    full_feature = [full_feature;fan_wks_tmp];
    %     gao_category{class}.sihks = gao_sihks_tmp;
    %     gao_category{class}.wks = gao_wks_tmp;
    %
    %     euclidean_category{class}.sihks = ec_sihks_tmp;
    %     euclidean_category{class}.wks = ec_wks_tmp;
    %
    %     nc_category{class}.sihks = nc_sihks_tmp;
    %     nc_category{class}.wks = nc_wks_tmp;
    %
    %     SM_category{class}.sihks = sm_sihks_tmp;
    %     SM_category{class}.wks = sm_wks_tmp;
    %
    %     ms_category{class}.sihks = ms_sihks_tmp;
    %     ms_category{class}.wks = ms_wks_tmp;
    %
    %     category{class}.sihks = sihks_tmp;
    %     category{class}.wks = wks_tmp;
    %     p_category{class}.sihks = p_sihks_tmp;
    %     p_category{class}.wks = p_wks_tmp;
end
% save(fullfile(base_path,['fan_category_',int2str(num_points),'.mat']),'fan_category');
% save(fullfile(base_path,['gao_category_',int2str(num_points),'.mat']),'gao_category');
% save(fullfile(base_path,['ec_category_',int2str(num_points),'.mat']),'euclidean_category');
% save(fullfile(base_path,['nc_category_',int2str(num_points),'.mat']),'nc_category');
% save(fullfile(base_path,['sm_category_',int2str(num_points),'.mat']),'SM_category');
% save(fullfile(base_path,['ms_category_',int2str(num_points),'.mat']),'ms_category');
% save(fullfile(base_path,['p_category_',int2str(num_points),'.mat']),'p_category');
% save(fullfile(base_path,['category_',int2str(num_points),'.mat']),'category');