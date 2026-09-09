%函数名称：Mean_Filter_Correlation
%函数功能：对输入的影像Array进行融合滤波，这是没有改进的预滤波，窗口大小是固定的
%函数输入：Array是待滤波影像，Correlation是与Array对应大小的相干系数影像
%函数输出：Array_Output是滤波后影像
% function Array_Output = Mean_Filter_Correlation_1(Array, Correlation)
% %        
% %         if(Correlation > 0.5)
% % 
%             A=fspecial('average',[3 3]);    
%             Array_Output=imfilter(Array,A); 
% % 
% %         else 
% %             A=fspecial('average',[4 4]);    
% %             Array_Output=imfilter(Array,A);                
% %         end


% end
function Array_Output = Mean_Filter_Correlation_baseline(Array,win_pre)
            j=sqrt(-1);
            Array = exp(j*Array);%将输入的相位转换到复数域   
            A=fspecial('average',[win_pre win_pre]);    
            Array_Output=imfilter(Array,A);%对复数均值滤波
            Array_Output = angle(Array_Output);%将复数转换为相位
end