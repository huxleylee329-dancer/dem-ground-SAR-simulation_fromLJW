%函数名称：AFilter
%函数功能：对输入的影像bArray进行融合滤波
%函数输入：Array是待滤波影像
%函数输出：bArrayF是滤波后影像
%对应论文第56页4.3.1自适应窗口预滤波式（4.6）
function Array_Output = Standard_Deviation(Array)

[Na,Nr]=size(Array);
Array_large=zeros(Na+2,Nr+2);
Array_large(2:Na+1,2:Nr+1)=Array;

Array_large(1,:) = Array_large(3,:); 
Array_large(end,:) = Array_large(end-2,:); 
Array_large(:,1) = Array_large(:,3);
Array_large(:,end) = Array_large(:,end-2);
Array_Output=zeros(Na,Nr);

for ia=2:Na+1
    for ir=2:Nr+1
        window=Array_large(ia-1:ia+1,ir-1:ir+1);%3*3的窗口
        uniform_phase=mean(mean(window));%对应论文第56页4.3.1自适应窗口预滤波式（4.6）
         Array_Output(ia-1,ir-1)=sqrt(sum(sum((window-uniform_phase).^2))./8);% Array_Output是矩阵每一个元素相较于均值的标准差
    end
end

end
