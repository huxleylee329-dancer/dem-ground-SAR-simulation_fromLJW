function [satetime] = getsatetime(Prf,numsate,path)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
Orbit=csvread([path]);%Matlab读取csv数据时，以0为开头第一行第一列
centertime=Orbit(round((length(Orbit(:,1))+1)/2));
one_secd=1/2/Prf;
satetime=(linspace(0,numsate-1,numsate)-(numsate-1)/2).*one_secd+centertime;
end

