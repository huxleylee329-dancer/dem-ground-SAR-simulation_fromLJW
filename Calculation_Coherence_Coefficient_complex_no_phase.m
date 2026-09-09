function [coherence] = Calculation_Coherence_Coefficient_complex_no_phase(image_master,image_slave)
% 计算相干系数
%% 相干系数图尺寸
wa = 3;            %窗口方位尺寸
wr = 3;            %窗口距离尺寸
[Na,Nr] = size(image_master);
win_a = (wa-1)/2;  %方位窗半径
win_r = (wr-1)/2;  %距离窗半径
Na_new = Na - 2 * win_a;
Nr_new = Nr - 2 * win_r;
coherence = zeros(Na_new+1,Nr_new+1);
%----------------- display information -------------------%
% display([repmat('-', 1, 3),' coherence coefficient ',repmat('-', 1, 4)])
% display(['window size:        ',num2str(wa),'*',num2str(wr)])
% display(repmat('-', 1, 30))
%% 像素叠加,实部虚部分别叠加平均
j = sqrt(-1);
for ii = 1+win_a : Na-win_a
    for jj = 1+win_r : Nr-win_r
        up = abs(sum(sum(image_master(ii-win_a : ii+win_a,jj-win_r : jj+win_r).*conj(image_slave(ii-win_a : ii+win_a,jj-win_r : jj+win_r)))));
        down = sqrt(sum(sum(abs(image_master(ii-win_a : ii+win_a,jj-win_r : jj+win_r)).^2)).*sum(sum(abs(image_slave(ii-win_a : ii+win_a,jj-win_r : jj+win_r)).^2)));
        coherence(ii,jj) = up/down;
    end
end
%function registration end
coherence(1,:)=[];
coherence(:,1)=[];
end