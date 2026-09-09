function w =Weight(a)
if (abs(a)<=1)
    w =1 -2*(abs(a))^2 +(abs(a))^3;
else if(abs(a)>1 && abs(a)<=2)
        w =4 -8*abs(a) +5*(abs(a))^2 -(abs(a))^3;
    else
        w=0;
    end
end