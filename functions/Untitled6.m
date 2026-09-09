% count = 0;
% accu = 0;
% for ii = 1:260
%     for jj = 1:260
%         if abs(delta(ii, jj)) >= 2*pi  
%             count = count + 1;
%         else
%             accu = accu + abs(delta(ii, jj));
%         end
%             
%     end
% end
% mean = accu/(260*260 - count)
% sum_s = 0;
% for ii = 1:260
%     for jj = 1:260
%         if abs(delta(ii, jj)) < 2*pi  
%             sum_s = sum_s + (abs(delta(ii, jj)) - mean)^2;
%         end
%             
%     end
% end
% sum_s = sqrt(sum_s/(260*260 - count))
% figure;imagesc(abs(delta/pi));axis('off');colorbar;
% ylabel(colorbar,'¾ø¶ÔÏàÎ»Îó²î(\pi)');
% set(gca,'LooseInset',get(gca,'TightInset'))

P_ground = zeros(260,260,3);
P_ground(:,:,1) = 6378931.02332382;
P_ground(:,:,2) = -993.003690212264;
P_ground(:,:,3) = 1030.37745707323;
