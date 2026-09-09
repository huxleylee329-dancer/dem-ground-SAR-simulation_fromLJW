function [image_slave_regis] = interp_cubic1(image_slave,a,b,c,d,e,f)
[nr,nc] = size(image_slave);
%扩充矩阵
NEW_image_slave = zeros(nr+3,nc+3);  %扩展3行,3列
%角点、边缘点
NEW_image_slave(1,1) = image_slave(1,1);
NEW_image_slave(1,2:nc+1) = image_slave(1,1:nc);
NEW_image_slave(1,nc+2) = image_slave(1,nc);
NEW_image_slave(1,nc+3) = image_slave(1,nc);
NEW_image_slave(2:nr+1,1) = image_slave(1:nr,1);
NEW_image_slave(nr+2,1) = image_slave(nr,1);
NEW_image_slave(nr+3,1) = image_slave(nr,1);
NEW_image_slave(nr+2,2:nc+1) = image_slave(nr,1:nc);
NEW_image_slave(nr+2,nc+2) = image_slave(nr,nc);
NEW_image_slave(nr+2,nc+3) = image_slave(nr,nc);
NEW_image_slave(nr+3,:) = NEW_image_slave(nr+2,:);
NEW_image_slave(2:nr+1,nc+2) = image_slave(1:nr,nc);
NEW_image_slave(nr+3,nc+2) = image_slave(nr,nc);
NEW_image_slave(nr+3,nc+3) = image_slave(nr,nc);
NEW_image_slave(:,nc+3) = NEW_image_slave(:,nc+2);
%内点
NEW_image_slave(2:nr+1,2:nc+1) = image_slave(1:nr,1:nc);
%行权
Row_weight = zeros(4,1);
Col_weight = zeros(1,4);
%实虚分别插值
image_slave_regis = zeros(nr,nc);
for ii = 2:nr+1
    for jj = 2:nc+1
        
        [move_a,move_r] = every_subpixel_move( ii-1,jj-1,a,b,c,d,e,f );
%         fprintf('(%f,%f) ',move_a,move_r);
        %每个像素的压像素偏移拟合
        Row_weight(1) = Calculation_Weight(1+move_a);
        Row_weight(2) = Calculation_Weight(move_a);
        Row_weight(3) = Calculation_Weight(1-move_a);
        Row_weight(4) = Calculation_Weight(2-move_a);
        %列权
        Col_weight(1) = Calculation_Weight(1+move_r);
        Col_weight(2) = Calculation_Weight(move_r);
        Col_weight(3) = Calculation_Weight(1-move_r);
        Col_weight(4) = Calculation_Weight(2-move_r);
        %权矩阵
        Weight = Row_weight * Col_weight;
        
        I_real = real(NEW_image_slave(ii-1:ii+2,jj-1:jj+2)); %实部插值  4*4
        I_imag = imag(NEW_image_slave(ii-1:ii+2,jj-1:jj+2)); %虚部插值
        Real = sum(sum(I_real .* Weight));
        Imag = sum(sum(I_imag .* Weight));
        image_slave_regis(ii-1,jj-1) = Real + j * Imag;
    end
%     fprintf('\n');
end

function weight = Calculation_Weight(offset)  %% 计算权重函数
if abs(offset)<1
    weight = 1 - 2 * abs(offset)^2 +  abs(offset)^3;
elseif (abs(offset)<=2)&&(abs(offset)>=1)
    weight = 4 - 8 * abs(offset) + 5 * abs(offset)^2 - abs(offset)^3;
elseif abs(offset)>2
    weight = 0;
end
end
end
