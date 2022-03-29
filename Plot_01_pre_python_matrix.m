clear;
clc;
% brainnet1 folder path 

rootdir = 'Z:\NN_FU\gPPI_Glasser_yihe\CPMresults\SSRT\B*';
group = 'BL';
% network label
if strcmp(group,'BL')
    % BL Glasser labels
    network_labels = 'Z:\NN_FU\gPPI_Glasser_yihe\CPMresults\GlasserLabel_for_CPM\BL_1.1\Glasser_networklabels.csv';
else  
    % FU_glasser network
    network_labels = 'Z:\NN_FU\gPPI_Glasser_yihe\CPMresults\GlasserLabel_for_CPM\FU2_1.08\Glasser_networklabels.csv';
end

threshold = [1,0.98,0.95];
labels = importdata(network_labels).data;
labels = labels(:,2);

% plot
rootfolder = dir(rootdir);
for g = 1: length(rootfolder)
    group_dir = dir(fullfile(rootfolder(g).folder,rootfolder(g).name,'CPM*'));

    for f = 1:length(group_dir)
        wk_dir = fullfile(group_dir(f).folder,group_dir(f).name);
        plot_folder = fullfile(wk_dir,'plot');
        
        for t = 1:length(threshold)

            mask_dir =  fullfile(plot_folder,['brainnet_',num2str(threshold(t))]);

            pos_matrix = importdata(fullfile(mask_dir,['Threshold',num2str(threshold(t)),'_pos_mask.txt']));
            neg_matrix = importdata(fullfile(mask_dir,['Threshold',num2str(threshold(t)),'_neg_mask.txt']));

            % calculate all edges for each modularity
            m = max(labels);
            pos_mod_matrix = zeros(m,m);
            neg_mod_matrix = zeros(m,m);
            combined_mod_matrix = zeros(m,m);
            for i = 1:m
                for n = 1:m
                    id_i = find(labels == i);
                    id_n = find(labels == n);
                    pos_m = pos_matrix(id_i,id_n);
                    pos_mod_matrix(i,n) = sum(sum(pos_m));

                    neg_m = neg_matrix(id_i,id_n);
                    neg_mod_matrix(i,n) = sum(sum(neg_m));  

                    combined_mod_matrix = pos_mod_matrix + neg_mod_matrix;
                end
            end

            writematrix(pos_mod_matrix,fullfile(mask_dir,'py_pos_mod_matrix.csv'),'Delimiter','tab');
            writematrix(neg_mod_matrix,fullfile(mask_dir,'py_neg_mod_matrix.csv'),'Delimiter','tab');
            writematrix(combined_mod_matrix,fullfile(mask_dir,'py_combined_mod_matrix.csv'),'Delimiter','tab');

        end
    end
end
