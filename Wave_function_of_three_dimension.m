clear all
close all
clc

dx=pi/60;
col=0:dx:pi;
az=0:dx:2*pi;
[phi,theta]=meshgrid(az,col);
% 计算 l=3 的网格上的 P_{l}^{m}(\cos \theta)
l=3;
Plm=legendre(l,cos(theta));
% 由于 legendre 为 m 的所有值计算答案，因此 Plm 会包含一些额外的函数值。
% 提取 m=2 的值并丢弃其余值。
% 使用 reshape 函数将结果定向为与 phi 和 theta 具有相同大小的矩阵。
m=2;
if l~=0
    Plm=reshape(Plm(m+1,:,:),size(phi));
end
% 计算Y_3^2的球谐函数值
a=(2*l+1)*factorial(l-m);
b=4*pi*factorial(l+m);
C=sqrt(a/b);
Ylm=C.*Plm.*exp(1i*m*phi);
% 球面坐标转换为笛卡尔坐标并绘图
[Xm,Ym,Zm]=sph2cart(phi,pi/2-theta,abs(real(Ylm)));
surf(Xm,Ym,Zm)
title('$Y_3^2$ spherical harmonic','interpreter','latex')
