clear;clc;
% plot brainnet viewer
% cpm_path_result = 'Z:\NN_FU\gPPI_Glasser_yihe\CPMresults\ICV\BL_0.1\output_LSO\IRV_CPM.mat';
rootdir = 'Z:\NN_FU\gPPI_Glasser_yihe\CPMresults\SSRT\B*';
CPM_name = 'SSRT_CPM.mat';
rootfolder = dir(rootdir);
% rootfolder(1:2,:) = [];
for g = 1: length(rootfolder)
    group_dir = dir(fullfile(rootfolder(g).folder,rootfolder(g).name,'CPM*'));
 
    for f = 1:length(group_dir)
        wk_dir = fullfile(group_dir(f).folder,group_dir(f).name);
        cpm_mat_file = fullfile(wk_dir,CPM_name);
        threshold = [1,0.98,0.95];
        plot_folder = fullfile(wk_dir,'plot');
        if exist(plot_folder,'dir')
            rmdir(plot_folder,'s')
        end
        for i = 1:length(threshold)
            % op_dir
            brainnet_op_dir =  fullfile(wk_dir,'plot',['brainnet_',num2str(threshold(i))]);
            % extract network under different threholds, save as a pos/neg_mask.txt(358 regions)
            threhold_cpm_network_weight(cpm_mat_file,wk_dir,brainnet_op_dir,threshold(i),0)
            % create pos/neg mask (360 regions) for connviwer plot
            create_connviewer_threshold(brainnet_op_dir,brainnet_op_dir,threshold(i));

        end
    end
end
