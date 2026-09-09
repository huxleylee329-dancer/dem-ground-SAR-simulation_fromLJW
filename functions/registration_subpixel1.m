function [image_slave_regis2] = registration_subpixel1(image_master,image_slave_regis)
% 精配准(亚像素级)
blocksize = 96; %子块大小128*128
[nr,nc] = size(image_master);
nrnew = floor(nr/blocksize) * blocksize;%向下取整
ncnew = floor(nc/blocksize) * blocksize;
% nrnew = 2^nextpow2(nr);
% ncnew = 2^nextpow2(nc);
% temp1 = zeros(nrnew,ncnew);               %主
% temp2 = zeros(nrnew,ncnew);               %辅
temp1 = image_master(1:nrnew,1:ncnew);      %主
temp2 = image_slave_regis(1:nrnew,1:ncnew); %辅,使块能取整
% 分块
image_master_cell = mat2cell(temp1,ones(nrnew/blocksize,1)*blocksize,ones(ncnew/blocksize,1)*blocksize);
image_slave_cell = mat2cell(temp2,ones(nrnew/blocksize,1)*blocksize,ones(ncnew/blocksize,1)*blocksize);
%mat2cell为矩阵分块，image_master_cell{1，1}表示第一块
% 子块总数
[nsubr,nsubc] = size(image_master_cell);
nsub   = nsubr * nsubc;

sub_r = zeros(nsub,1);    sub_c = zeros(nsub,1);
m = zeros(nsub,1);        n = zeros(nsub,1);     indx = zeros(nsub,1);

% 按列顺序依次提取各子块
for nn = 1:nsub
    image1 = cell2mat(image_master_cell(nn));
    image2 = cell2mat(image_slave_cell(nn));
    % 复图像插值,补零插值,可增采样
    interptimes = 32;
    [image_interp1] = interp_paddingzero(image1,interptimes);
    [image_interp2] = interp_paddingzero(image2,interptimes);%image_interp2的大小为原图像长和宽的interptimes倍
    % 实相关函数求取偏移量
    [move_r,move_c,~] = real_coherent(image_interp1,image_interp2,0);
    % 亚像素偏移量
    sub_r(nn) = move_r/interptimes;
    sub_c(nn) = move_c/interptimes;
    
    jj = floor(nn/nsubr);    ii = rem(nn,nsubr); % rem是取余，floor向负无穷大方向取整
    n(nn) = jj*blocksize + blocksize/2;%n为每块的中心像素的行数
    if ii==0
        m(nn) = nsubr*blocksize + blocksize/2;%m为每块中心像素的列数？？认为加号应为减号
    else
        m(nn) = (ii-1)*blocksize + blocksize/2;
    end
    %图像第一次插值，找出相关系数高的子块
    image2_interp = interp_cubic(image2,sub_r(nn),sub_c(nn));
    [~,coherence] = Calculation_Coherence_Coefficient(image1,image2_interp);
    if mean(mean(coherence)) > 0.5
        indx(nn) = 1;
    end
end
idx = (indx > 0);
m = m(idx);    n = n(idx);
sub_r = sub_r(idx);
sub_c = sub_c(idx);
[a,b,c,d,e,f] = all_subpixel_move(m,n,sub_r,sub_c);
% 图像插值
image_slave_regis2 = interp_cubic1(image_slave_regis,a,b,c,d,e,f);
%% 输出
end