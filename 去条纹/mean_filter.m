function mean_array=mean_filter(Array, Size_large)

[Na,Nr]=size(Array);
Array_large=zeros(Na+Size_large-1,Nr+Size_large-1);
radius=(Size_large-1)/2;
Array_large(radius+1:end-radius,radius+1:end-radius)=Array;
Array_large(1:radius,1:end) = flipud(Array_large(radius+1:2*radius,1:end)); 
Array_large(1:end,1:radius) = fliplr(Array_large(1:end,radius+1:2*radius));  
Array_large(end-radius+1:end,1:end) = flipud(Array_large(end-2*radius+1:end-radius,1:end));
Array_large(1:end,end-radius+1:end) = fliplr(Array_large(1:end,end-2*radius+1:end-radius));



window_sum=zeros(Na,Nr);

for m=1:Size_large
    for n=1:Size_large
        window=Array_large(m:end-(Size_large-m),n:end-(Size_large-n));      
        window_sum=window_sum+window;
    end
end

mean_array=window_sum./Size_large./Size_large;

end