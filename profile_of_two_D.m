% 绘制主量子数为1、2、3的氢原子波函数
close all;
clear all;
clc;

% 创建一个新的图形窗口
fig = figure;
fig.Position = [0, 0, 1600, 900];
fig.Color = [0, 0, 0];
set(fig, 'InvertHardCopy', 'off');

% 定义颜色映射
map = [0.1294 0.0549 0.1725; 0.2196 0.1608 0.2902; 0.3882 0.1804 0.4941;
       0.4392 0.1922 0.4706; 0.5333 0.2235 0.4392; 0.6471 0.2588 0.3686;
       0.7137 0.2745 0.3294; 0.7725 0.3059 0.2902; 0.8510 0.3725 0.2275;
       0.9137 0.4196 0.1804; 0.9608 0.5020 0.2000; 0.9765 0.5529 0.2078;
       0.9804 0.6431 0.2549; 0.9843 0.6627 0.2706; 0.9765 0.7176 0.3412;
       0.9765 0.7686 0.4000; 0.9765 0.8118 0.4902; 0.9725 0.8510 0.5961;
       0.9882 0.9020 0.6667; 1.0000 0.9451 0.8431; 1.0000 0.9961 0.9804;
       1.0000 1.0000 1.0000];
Xi = 1:size(map, 1);
Xq = linspace(1, size(map, 1), 800);
map = [interp1(Xi, map(:, 1), Xq, 'linear')', ...
       interp1(Xi, map(:, 2), Xq, 'linear')', ...
       interp1(Xi, map(:, 3), Xq, 'linear')'];

% 绘制公式部分
formula_ax = axes('Position', [0.50, 0.7, 0.2, 0.2]);
formula_ax.XLim = [0, 3];
formula_ax.YLim = [0, 1];
formula_ax.XTick = [];
formula_ax.YTick = [];
formula_ax.XColor = [0, 0, 0];
formula_ax.YColor = [0, 0, 0];
formula_ax.Color = [0 0 0];
text(formula_ax, .04, 0.8, 'Hydrogen   Electron   Orbitals', 'Color', [1, 1, 1] * .98, ...
    'FontName', 'cambria', 'FontWeight', 'bold', 'FontSize', 40)
text(formula_ax, .04, 0.2, '$\psi_{n, l, m}=e^{-r / n}\left(\frac{2 r}{n}\right)^{l}\left[L_{n-l-1}^{2 l+1}\left(\frac{2 r}{n}\right)\right] Y_{l}^{m}(\theta, \phi)$', ...
    'Color', [1, 1, 1] * .98, 'FontSize', 30, 'Interpreter', 'latex')

% 定义量子数组合
quantum_numbers = {
    1, 0, 0;  % n=1, l=0, m=0
    2, 0, 0;  % n=2, l=0, m=0
    2, 1, 0;  % n=2, l=1, m=0
    2, 1, 1;  % n=2, l=1, m=1
    3, 0, 0;  % n=3, l=0, m=0
    3, 1, 0;  % n=3, l=1, m=0
    3, 1, 1;  % n=3, l=1, m=1
    3, 2, 0;  % n=3, l=2, m=0
    3, 2, 1;  % n=3, l=2, m=1
    3, 2, 2   % n=3, l=2, m=2
};

% 绘制每个量子数对应的波函数
row_positions = [0.775, 0.525, 0.275, 0.025];  % 每行的垂直位置
col_positions = [0.05, 0.3, 0.55, 0.8];  % 每列的水平位置
row_height = 0.2;  % 每行的高度
col_width = 0.2;  % 每列的宽度

for i = 1:size(quantum_numbers, 1)
    n = quantum_numbers{i, 1};
    l = quantum_numbers{i, 2};
    m = quantum_numbers{i, 3};
    
    % 计算子图位置
    if i == 1
        row = 1;
        col = 1;
    elseif i <= 3
        row = 2;
        col = i - 1;
    elseif i <= 6
        row = 3;
        col = i - 3;
    else
        row = 4;
        col = i - 6;
    end
    
    ax = axes('Position', [col_positions(col), row_positions(row), col_width, row_height]);
    
    % 生成网格
    [X, Z] = meshgrid(linspace(-30, 30, 120));
    Y = zeros(size(X));
    
    % 计算波函数
    psi = HydWave(n, l, m, X, Y, Z);
    
    % 绘制波函数
    surf(X, Z, psi, 'EdgeColor', 'none');
    colormap(map);
    view(2);
    axis tight;
    caxis([0, max(max(psi))]);
    
    % 设置坐标轴属性
    ax.XLim = [-30, 30];
    ax.YLim = [-30, 30];
    ax.DataAspectRatio = [1, 1, 1];
    ax.XTick = [];
    ax.YTick = [];
    ax.XColor = [1, 1, 1] * 0.4;
    ax.YColor = [1, 1, 1] * 0.4;
    ax.LineWidth = 1.5;
    ax.Box = 'on';
    ax.FontName = 'cambria';
    
    % 在子图下方添加量子数标注
    title_text = ['(', num2str(n), ',', num2str(l), ',', num2str(m), ')'];
    ax.Title.String = title_text;
    ax.Title.FontSize = 10;
    ax.Title.Color = [1, 1, 1] * 0.9;
end

% 保存图像
if ~exist('picture/Electron_contour_of_the_wave_function_of_hydrogen_atom', 'dir')
    mkdir picture/Electron_contour_of_the_wave_function_of_hydrogen_atom;  % 如果文件夹不存在，则创建文件夹
end
saveas(fig, 'picture/Electron_contour_of_the_wave_function_of_hydrogen_atom/Hydrogen_Orbitals.png');  % 保存到 picture 文件夹

% 氢原子波函数计算函数
function psi = HydWave(n, l, m, x, y, z)
    % 由坐标点计算向量模长及角度
    r = vecnorm([x(:), y(:), z(:)]')';
    theta = atan2(vecnorm([x(:), y(:)]')', z(:));
    phi = atan2(y(:), x(:));
    
    % 恢复矩阵形状
    r = reshape(r, size(x));
    theta = reshape(theta, size(x));
    phi = reshape(phi, size(x));
    
    % 利用MATLAB自带legendre函数计算球谐函数
    Plm = legendre(l, cos(theta));
    if l ~= 0
        Plm = reshape(Plm(m + 1, :, :), size(phi));
    end
    C = sqrt(((2 * l + 1) * factorial(l - m)) / (4 * pi * factorial(l + m)));
    Ylm = C .* Plm .* exp(1i * m * phi);
    
    % laguerreL函数部分计算
    if n == l + 1
        Lag = ones(size(r));
    else
        Lag = laguerreL(n - l - 1, 2 * l + 1, 2 .* r ./ n);
    end
    
    % 整合起来
    psi = exp(-r ./ n) .* (2 .* r ./ n) .^ l .* Lag .* Ylm;
    psi = real(conj(psi) .* psi);
end

% Laguerre多项式计算函数
function L = laguerreL(n, alpha, x)
    if n == 0
        L = ones(size(x));
    elseif n == 1
        L = 1 + alpha - x;
    else
        L0 = ones(size(x));
        L1 = 1 + alpha - x;
        for k = 2:n
            L = ((2 * k + alpha - 1 - x) .* L1 - (k + alpha - 1) .* L0) / k;
            L0 = L1;
            L1 = L;
        end
    end
end