clear all
close all
clc

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

% 确保 picture 文件夹存在
if ~exist('picture/Electron_contour_of_the_wave_function_of_hydrogen_atom', 'dir')
    mkdir picture/Electron_contour_of_the_wave_function_of_hydrogen_atom;
end

% 绘制每个量子数对应的波函数并单独保存
for i = 1:size(quantum_numbers, 1)
    n = quantum_numbers{i, 1};
    l = quantum_numbers{i, 2};
    m = quantum_numbers{i, 3};
    
    % 创建新图形窗口
    fig = figure;
    fig.Position = [100, 100, 800, 600];
    fig.Color = [0, 0, 0];
    
    % 绘制球谐函数部分
    dx = pi/60;
    col = 0:dx:pi;
    az = 0:dx:2*pi;
    [phi, theta] = meshgrid(az, col);
    
    % 计算球谐函数
    Plm = legendre(l, cos(theta));
    if l ~= 0
        Plm = reshape(Plm(m + 1, :, :), size(phi));
    end
    
    a = (2*l + 1) * factorial(l - m);
    b = 4 * pi * factorial(l + m);
    C = sqrt(a / b);
    Ylm = C .* Plm .* exp(1i * m * phi);
    
    % 球面坐标转换为笛卡尔坐标并绘图
    [Xm, Ym, Zm] = sph2cart(phi, pi/2 - theta, abs(real(Ylm)));
    surf(Xm, Ym, Zm);
    colormap(map);
    shading interp;
    
    % 设置图形属性
    title(['$Y_{', num2str(l), '}^{', num2str(m), '}$ spherical harmonic for n = ', num2str(n)], ...
          'Interpreter', 'latex', 'FontSize', 14, 'Color', [1, 1, 1] * 0.9);
    axis equal;
    axis tight;
    view(3);
    grid off;
    
    % 保存图像到 picture 文件夹
    filename = sprintf('picture/Electron_contour_of_the_wave_function_of_hydrogen_atom/Hydrogen_Orbital_n%d_l%d_m%d.png', n, l, m);
    saveas(fig, filename);
    
end

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