function [return_value] = after_regis_image_index_update(move_r,move_c,process_ID_master,process_ID_slave)
%% 更新配准后的主辅图像坐标范围。
fid1 = fopen(['..\output\',process_ID_master,'/SSC_infor_2.txt']);
result = textscan(fid1,'%s %f');
fclose(fid1);
row1   = result{2}(1,1);
rowend1 = result{2}(2,1);
col1   = result{2}(3,1);
colend1 = result{2}(4,1);
fid1 = fopen(['..\output\',process_ID_slave,'/SSC_infor_2.txt']);
result = textscan(fid1,'%s %f');
fclose(fid1);
row2   = result{2}(1,1);
rowend2 = result{2}(2,1);
col2   = result{2}(3,1);
colend2 = result{2}(4,1);


% 更新主辅图像行坐标
if move_r>0 %辅图像已向上搬移
   rowend1 = rowend1 - move_r;
   row2 = row2 + move_r;
else
   row1 = row1 - move_r;    
   rowend2 = rowend2 + move_r;
end
% 更新主辅图像列坐标
if move_c>0 %辅图像已向左搬移
   colend1 = colend1 - move_c;
   col2 = col2 + move_c;
else
   col1 = col1 - move_c;    
   colend2 = colend2 + move_c;
end
% 写入文件  
fid1 = fopen(['..\output\',process_ID_master,'/SSC_infor_3.txt'],'wt');
fprintf(fid1,'%s\n%i \n','裁剪起始行: ',row1);
fprintf(fid1,'%s\n%i \n','裁剪终止行: ',rowend1);
fprintf(fid1,'%s\n%i \n','裁剪起始列: ',col1);
fprintf(fid1,'%s\n%i \n','裁剪终止列: ',colend1);
fclose(fid1);
fid1 = fopen(['..\output\',process_ID_slave,'/SSC_infor_3.txt'],'wt');
fprintf(fid1,'%s\n%i \n','裁剪起始行: ',row2);
fprintf(fid1,'%s\n%i \n','裁剪终止行: ',rowend2);
fprintf(fid1,'%s\n%i \n','裁剪起始列: ',col2);
fprintf(fid1,'%s\n%i \n','裁剪终止列: ',colend2);
fclose(fid1);
return_value = ['已更新配准后的图像位置索引'];
end