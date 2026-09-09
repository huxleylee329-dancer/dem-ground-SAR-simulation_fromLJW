%函数名称：Mean_Filter_Correlation
%函数功能：通过调整gamma和PSD_Window_mean改变window大小，与极限窗口estimate_window作比较...
%取较小的窗口后，对输入的影像Array进行均值滤波
%函数输入：Array是待滤波影像，Correlation是与Array对应大小的相干系数影像
%函数输出：Array_Output是滤波后影像
function [Array_Output,window] = Mean_Filter_Correlation(Array,gamma,PSD_Window,estimate_window)
%        
%         if(Correlation > 0.5)
%          
%             A=fspecial('average',[3 3]);    
%             Array_Output=imfilter(Array,A); 
% 
%         else 
%             A=fspecial('average',[4 4]);    
%             Array_Output=imfilter(Array,A);                
%         end
            [a,b]=size(PSD_Window);%PSD_Window表示估计窗口每一个像素的PSD
            PSD_Window_mean = sum(sum(PSD_Window))/(a*b);%估计窗口的PSD均值
            windowsize = round(1/gamma+PSD_Window_mean+1);%round是对矩阵元素四舍五入至最近的整数
            %Q:按照式（4.5），半径=（windowsize-1）/2=1/gamma+PSD_Window_mean
            %所以windowsize=2/gamma+2*PSD_Window_mean+1
%             windowsize = round(2/gamma+2*PSD_Window_mean+1);%round是对矩阵元素四舍五入至最近的整数
            
            if mod(windowsize,2)==0%mod为取余运算
            window = windowsize+1;
            else window = windowsize;%保证窗口大小为奇数
            end
            %根据论文57页（4.7）（4.8）可知，estimate_window为极限多视倍数，由系统参数和控制点先验信息估计得到 
            
            %取最小的窗口进行均值滤波;
            if (window<=estimate_window)
            A = fspecial('average',[window window]);    
            Array_Output = imfilter(Array,A); 
            else
            A = fspecial('average',[estimate_window estimate_window]);    
            Array_Output = imfilter(Array,A); 
            end 
end