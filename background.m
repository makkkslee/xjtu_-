clear all
close all
clc

function HydWaveDemo
% @author:slandarer

% 构造一个好看的颜色映射实属不易
map=[0.1294 0.0549 0.1725;0.2196 0.1608 0.2902;0.3882 0.1804 0.4941;
     0.4392 0.1922 0.4706;0.5333 0.2235 0.4392;0.6471 0.2588 0.3686;
     0.7137 0.2745 0.3294;0.7725 0.3059 0.2902;0.8510 0.3725 0.2275;
     0.9137 0.4196 0.1804;0.9608 0.5020 0.2000;0.9765 0.5529 0.2078;
     0.9804 0.6431 0.2549;0.9843 0.6627 0.2706;0.9765 0.7176 0.3412;
     0.9765 0.7686 0.4000;0.9765 0.8118 0.4902;0.9725 0.8510 0.5961;
     0.9882 0.9020 0.6667;1.0000 0.9451 0.8431;1.0000 0.9961 0.9804;
     1.0000 1.0000 1.0000];
Xi=1:size(map,1);Xq=linspace(1,size(map,1),800);
map=[interp1(Xi,map(:,1),Xq,'linear')',...
     interp1(Xi,map(:,2),Xq,'linear')',...
     interp1(Xi,map(:,3),Xq,'linear')'];
 
% 情况信息列表
condList=[1,2,0,0,-10,10,0.01;
          2,3,0,0,-20,20,0.01;
          6,2,1,0,-12,12,0.03;
          7,3,1,0,-20,20,0.08;
          8,3,1,1,-20,20,0.08;
          11,2,1,1,-12,12,0.03;
          12,3,2,0,-23,23,0.35;
          13,3,2,1,-23,23,0.35;
          14,3,2,2,-23,23,0.35;
          16,4,0,0,-36,36,0.008;
          17,4,1,0,-30,30,0.2;
          18,4,1,1,-30,30,0.1;
          19,4,2,0,-30,30,0.95;
          20,4,2,1,-30,30,0.95;
          21,4,2,2,-30,30,0.95;
          22,4,3,0,-35,35,6;
          23,4,3,1,-35,35,6;
          24,4,3,2,-35,35,6;
          25,4,3,3,-35,35,1.8];

fig=gcf;
fig.Color=[0,0,0];
set(fig,'InvertHardCopy','off');

% 文本信息
axh=subplot(5,5,[3,4,5]);
axh.XLim=[0,3];
axh.YLim=[0,1];
axh.XTick=[];
axh.YTick=[];
axh.XColor=[0,0,0];
axh.YColor=[0,0,0];
axh.Color=[0 0 0];
text(axh,.04,0.8,'Hydrogen   Electron   Orbita','Color',[1,1,1].*.98,...
    'FontName','cambria','FontWeight','bold','FontSize',15)
text(axh,.04,0.2,'$\psi_{n, l, m}=e^{-r / n}\left(\frac{2 r}{n}\right)^{l}\left[L_{n-l-1}^{2 l+1}\left(\frac{2 r}{n}\right)\right] Y_{l}^{m}(\theta, \phi)$',...
    'Color',[1,1,1].*.98,'FontSize',12,'Interpreter','latex')

% 文本信息2
axh2=subplot(5,5,[9,10]);
axh2.XLim=[0,2];
axh2.YLim=[0,1];
axh2.XTick=[];
axh2.YTick=[];
axh2.XColor=[0,0,0];
axh2.YColor=[0,0,0];
axh2.Color=[0 0 0];
text(axh2,.04,0.8,'Not normalized by:','Color',[1,1,1].*.98,...
    'FontName','cambria','FontSize',13)
text(axh2,.04,0.2,'$\sqrt{\left(\frac{2}{n a}\right)^{3} \frac{(n-l-1) !}{2 n(n+l) !}}$','Color',[1,1,1].*.98,...
    'FontName','cambria','FontSize',12,'Interpreter','latex')

% 循环绘图
for i=1:size(condList,1)
    ax=subplot(5,5,condList(i,1));
    
    % 绘制密度分布
    [X,Z]=meshgrid(linspace(condList(i,5),condList(i,6),120));
    Y=zeros(size(X));
    psi=HydWave(condList(i,2),condList(i,3),condList(i,4),X,Y,Z);
    surf(X,Z,psi,'EdgeColor','none')
    caxis(ax,[0,condList(i,7)]);colormap(map);view(2)
    axis tight
    
    % 坐标区域修饰
    ax.XLim=[condList(i,5),condList(i,6)];
    ax.YLim=[condList(i,5),condList(i,6)];
    ax.DataAspectRatio=[1,1,1];
    ax.XTick=[];
    ax.YTick=[];
    ax.XColor=[1,1,1].*.4;
    ax.YColor=[1,1,1].*.4;
    ax.LineWidth=1.5;
    ax.Box='on';
    ax.FontName='cambria';
    ax.XLabel.String=['(',num2str(condList(i,2)),',',...
                   num2str(condList(i,3)),',',...
                   num2str(condList(i,4)),')'];
    drawnow
end
saveas(fig,'result.png')
end
