clear all
close all
clc

[X,Z]=meshgrid(linspace(-30,30,120));
Y=zeros(size(X));

psi=HydWave(4,2,0,X,Y,Z);
surf(X,Z,psi,'EdgeColor','none')
view(2);caxis([0,2])


function psi=HydWave(n,l,m,x,y,z)
% @author:slandarer

% 由坐标点计算向量模长及角度
r=vecnorm([x(:),y(:),z(:)]')';
theta=atan2(vecnorm([x(:),y(:)]')',z(:));
phi=atan2(y(:),x(:));

% 恢复矩阵型状
r=reshape(r,size(x));
theta=reshape(theta,size(x));
phi=reshape(phi,size(x));

% 利用MATLAB自带legendre函数计算球谐函数
Plm=legendre(l,cos(theta));
if l~= 0
    Plm=reshape(Plm(m+1,:,:),size(phi));
end
C=sqrt(((2*l+1)*factorial(l-m))/(4*pi*factorial(l+m)));
Ylm=C.*Plm.*exp(1i*m*phi);

% laguerreL函数部分计算
Lag=laguerreL(n-l-1,2*l+1,2.*r./n);

% 整合起来
psi=exp(-r./n).*(2.*r./n).^l.*Lag.*Ylm;
psi=real(conj(psi).*psi);
end

