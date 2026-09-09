function [return_value] = cut_image_index_update(nr1,nc1,nr,nc,process_ID_master)
delta_r    = nr1 - nr;
delta_c    = nc1 - nc;
fid1 = fopen(['..\output\',process_ID_master,'/SSC_infor.txt']);
result = textscan(fid1,'%s %f');
fclose(fid1);
row1   = result{2}(3,1);
rowend = result{2}(4,1);
col1   = result{2}(5,1);
colend = result{2}(6,1);
rowend = rowend - delta_r;
colend = colend - delta_c;
fid1 = fopen(['..\output\',process_ID_master,'/SSC_infor_2.txt'],'wt');
fprintf(fid1,'%s %i \n','裁剪起始行: ',row1);
fprintf(fid1,'%s %i \n','裁剪终止行: ',rowend);
fprintf(fid1,'%s %i \n','裁剪起始列: ',col1);
fprintf(fid1,'%s %i \n','裁剪终止列: ',colend);
fclose(fid1);
return_value = ['已更新图像位置索引'];
end