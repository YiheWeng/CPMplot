function create_connviewer_threshold(mask_path,op_dir,freq_thresh)
%% create pos/neg mask (360 regions) for connviwer plot
% mask_path -- the folder where pos/neg_mask.txt are saved 
% freq_thresh -- threshold range:0-1
% op_dir = 'Z:\NN_FU\gPPI_Glasser_yihe\CPM_FU2\IRV\FU_IRV\output\Relate0.01_node154\relate_plot\brainnet';


    matrix_pos = importdata(fullfile(mask_path, ['Threshold',num2str(freq_thresh),'_pos_mask.txt']));
    matrix_neg = importdata(fullfile(mask_path, ['Threshold',num2str(freq_thresh),'_neg_mask.txt']));

     matrix_pos_final = [matrix_pos(1:92,:);zeros(1,358);matrix_pos(93:end,:)];
     matrix_pos_final = [matrix_pos_final(1:272,:);zeros(1,358);matrix_pos_final(273:end,:)];
     matrix_pos_final = [matrix_pos_final(:,1:92),zeros(360,1),matrix_pos_final(:,93:end)];
     matrix_pos_final = [matrix_pos_final(:,1:272),zeros(360,1),matrix_pos_final(:,273:end)];
    % writematrix(matrix_pos_final,'Fail_pos_matrix.txt');
    writematrix(matrix_pos_final,fullfile(op_dir,['Wholebrain360_Threshold',num2str(freq_thresh),'_pos_mask.txt']));


     matrix_neg_final = [matrix_neg(1:92,:);zeros(1,358);matrix_neg(93:end,:)];
     matrix_neg_final = [matrix_neg_final(1:272,:);zeros(1,358);matrix_neg_final(273:end,:)];
     matrix_neg_final = [matrix_neg_final(:,1:92),zeros(360,1),matrix_neg_final(:,93:end)];
     matrix_neg_final = [matrix_neg_final(:,1:272),zeros(360,1),matrix_neg_final(:,273:end)];
    % writematrix(matrix_neg_final,'Fail_neg_matrix.txt');
    writematrix(matrix_neg_final,fullfile(op_dir,['Wholebrain360_Threshold',num2str(freq_thresh),'_neg_mask.txt']));

   
end 