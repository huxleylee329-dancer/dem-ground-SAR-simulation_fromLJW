function [sigmaCs] = reflectivity(AngleCs)
%由得到的入射角计算出RCS
sigma1=(double(AngleCs<=17)).*(2-0.7058823529412*AngleCs);
sigma2=(double((17<AngleCs)&(AngleCs<=30))).*(-10-0.3846153846154*(AngleCs-17));
sigma3=(double((30<AngleCs)&(AngleCs<=80))).*(-15-0.2*(AngleCs-30));
sigma4=(double((80<AngleCs)&(AngleCs<=90))).*(-25-7.5*(AngleCs-80));
sigma5=(double(AngleCs>90))*(-200);
sigmaCs=sigma1+sigma2+sigma3+sigma4+sigma5;
end

