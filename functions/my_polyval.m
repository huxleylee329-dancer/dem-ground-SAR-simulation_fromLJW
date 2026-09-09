function index = my_polyval(x,y,coef);
index = coef(1)*ones(length(x))+coef(2)*x+coef(3)*y+coef(4)*x.*x+coef(5)*x.*y+coef(6)*y.*y;
end