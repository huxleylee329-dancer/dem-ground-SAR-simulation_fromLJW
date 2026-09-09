function [rou,coherence] = Calculation_Coherence_Coefficient(image_master,image_slave)
% 计算相干系数
%% 相干系数图尺寸
wa = 3;            %窗口方位尺寸
wr = 3;            %窗口距离尺寸
[Na,Nr] = size(image_master);
win_a = (wa-1)/2;  %方位窗半径
win_r = (wr-1)/2;  %距离窗半径
Na_new = Na - 2 * win_a;
Nr_new = Nr - 2 * win_r;
coherence = zeros(Na_new,Nr_new);
rou = zeros(Na_new,Nr_new);
%----------------- display information -------------------%
% display([repmat('-', 1, 3),' coherence coefficient ',repmat('-', 1, 4)])
% display(['window size:        ',num2str(wa),'*',num2str(wr)])
% display(repmat('-', 1, 30))
%% 像素叠加,实部虚部分别叠加平均

for ii = 1+win_a : Na-win_a
    for jj = 1+win_r : Nr-win_r
        s1 = abs(image_master(ii-win_a : ii+win_a,jj-win_r : jj+win_r));
        s2 = abs(image_slave(ii-win_a : ii+win_a,jj-win_r : jj+win_r));
        rou(ii-win_a,jj-win_r)  = sum(sum(s1.^2.*s2.^2))/sqrt(sum(sum(s1.^4))*sum(sum(s2.^4)));
        coherence(ii-win_a,jj-win_r) = rou(ii-win_a,jj-win_r);
%         if rou(ii-win_a,jj-win_r)  >= 0.5
%            coherence(ii-win_a,jj-win_r) = sqrt(2*rou(ii-win_a,jj-win_r) -1);
%         else
%            coherence(ii-win_a,jj-win_r) = 0;
%         end
    end
end
%function registration end
end