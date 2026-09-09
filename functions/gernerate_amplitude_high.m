%---------------产生高频幅度图---------------
display('产生幅度图......');
in_folder = '.\Kernel\Input\high\';
out_folder = '.\Kernel\Input\high\';
load([in_folder 'mast.mat']);
load([in_folder 'slave.mat']);
mast = mapminmax(abs(mast), 0, 1);
slave = mapminmax(abs(slave), 0, 1);
imwrite(mast, [out_folder 'mast_amp.jpg'], 'jpg');
imwrite(slave, [out_folder 'slave_amp.jpg'], 'jpg');
display('完成！');