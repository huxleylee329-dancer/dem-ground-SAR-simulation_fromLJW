function [image_master_regis,image_slave_regis,move_r,move_c] = registration_pixel(image_master,image_slave,nr,nc)
% ´ÖÅä×¼(ÏñËØ¼¶)
% ÊµÏà¹ØÏµÊıº¯ÊıÇóÈ¡Æ«ÒÆÁ¿
[move_r,move_c,~] = real_coherent(image_master,image_slave,1);
%
image_slave_temp = image_slave;
% ĞĞÆ«ÒÆ£¨ÊúÖ±ÒÆ¶¯£©
if move_r >= 0
   % ¸¨Í¼ÏñÏòÉÏ°áÒÆ
   image_slave_temp(1:nr-move_r,:) = image_slave(move_r+1:nr,:);
   % ²Ã¼ô
   image_master_mid = image_master(1:nr-move_r,:);
   image_slave_mid  = image_slave_temp(1:nr-move_r,:);
elseif move_r < 0
   % ¸¨Í¼ÏñÏòÏÂ°áÒÆ
   image_slave_temp(1-move_r:nr,:) = image_slave(1:nr+move_r,:);
   % ²Ã¼ô
   image_master_mid = image_master(1-move_r:nr,:);
   image_slave_mid  = image_slave_temp(1-move_r:nr,:);
end
%
image_slave_temp_mid = image_slave_mid;
% ÁĞÆ«ÒÆ£¨Ë®Æ½ÒÆ¶¯£©
if move_c >= 0
   % ¸¨Í¼ÏñÏò×ó°áÒÆ
   image_slave_temp_mid(:,1:nc-move_c) = image_slave_mid(:,move_c+1:nc);
   % ²Ã¼ô
   image_master_regis = image_master_mid(:,1:nc-move_c);
   image_slave_regis  = image_slave_temp_mid(:,1:nc-move_c);
elseif move_c < 0
   % ¸¨Í¼ÏñÏòÓÒ°áÒÆ
   image_slave_temp_mid(:,1-move_c:nc) = image_slave_mid(:,1:nc+move_c);
   % ²Ã¼ô
   image_master_regis = image_master_mid(:,1-move_c:nc);
   image_slave_regis  = image_slave_temp_mid(:,1-move_c:nc);
end
end