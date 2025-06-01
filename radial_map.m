clear all
close all
clc

% 定义主量子数范围
n_values = [1, 2, 3];

% 第一玻尔半径（单位：皮米）
a0 = 52.9;

% 确保 picture/radial_map 文件夹存在
if ~exist('picture', 'dir')
    mkdir picture;
end
if ~exist('picture/radial_map', 'dir')
    mkdir('picture/radial_map');
end

% 绘制汇总图
fig = figure;
fig.Position = [100, 100, 800, 600];
fig.Color = [1, 1, 1];

% 定义归一化半径范围 r/a0，确保覆盖所有主量子数的范围
max_r_normalized = max(n_values)^2 * 5;
r_normalized = linspace(0, max_r_normalized, 400);
r = r_normalized * a0;  % 对应的实际半径（单位：皮米）

hold on;
all_max_P = [];

for n = n_values
    l_values = 0:n-1;
    for l = l_values
        % 计算径向概率分布
        P = zeros(size(r));
        for i = 1:length(r)
            P(i) = RadialProbability(n, l, r(i), r);
        end
        
        % 绘制当前角量子数的径向概率分布曲线
        plot(r_normalized, P, 'LineWidth', 1.5, 'DisplayName', ['n = ', num2str(n), ', l = ', num2str(l)]);
        
        % 记录当前曲线的最大值
        all_max_P = [all_max_P, max(P)];
    end
end
hold off;

% 设置图形属性
xlabel('Radial Distance (r/a_0)', 'FontSize', 12, 'Color', [0, 0, 0]);
ylabel('Radial Probability Density', 'FontSize', 12, 'Color', [0, 0, 0]);
title('Radial Probability Distribution for n = 1, 2, 3', ...
    'Interpreter', 'latex', 'FontSize', 14, 'Color', [0, 0, 0]);
grid on;
legend('show');
axis tight;

% 调整坐标轴范围以确保图像美观
ylim([0 max(all_max_P) * 1.1]);

% 添加注释说明单位转换
text(0.5, 0.95, ['First Bohr Radius (a_0) = ', num2str(a0), ' pm'], 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0, 0, 0]);

% 保存图像到 picture/radial_map 文件夹
filename = 'picture/radial_map/Radial_Probability_Summary.png';
saveas(fig, filename);

% 循环绘制每个主量子数的径向概率分布图
for n = n_values
    max_r_normalized = n^2 * 5;  % 归一化后的最大半径范围
    r_normalized = linspace(0, max_r_normalized, 400);
    r = r_normalized * a0;  % 实际半径（单位：皮米）

    fig = figure;
    fig.Position = [100, 100, 800, 600];
    fig.Color = [1, 1, 1];

    l_values = 0:n-1;

    hold on;
    for l = l_values
        % 计算径向概率分布
        P = zeros(size(r));
        for i = 1:length(r)
            P(i) = RadialProbability(n, l, r(i), r);
        end

        % 绘制当前角量子数的径向概率分布曲线
        plot(r_normalized, P, 'LineWidth', 1.5, 'DisplayName', ['l = ', num2str(l)]);
    end
    hold off;

    % 设置图形属性
    xlabel('Radial Distance (r/a_0)', 'FontSize', 12, 'Color', [0, 0, 0]);
    ylabel('Radial Probability Density', 'FontSize', 12, 'Color', [0, 0, 0]);
    title(['Radial Probability Distribution for n = ', num2str(n)], ...
        'Interpreter', 'latex', 'FontSize', 14, 'Color', [0, 0, 0]);
    grid on;
    legend('show');
    axis tight;

    % 添加注释说明单位转换
    text(0.5, 0.95, ['First Bohr Radius (a_0) = ', num2str(a0), ' pm'], 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0, 0, 0]);

    % 调整坐标轴范围以确保图像美观
    ylim([0 max(P)*1.1]);

    % 保存图像到 picture/radial_map 文件夹
    filename = sprintf('picture/radial_map/Radial_Probability_n%d.png', n);
    saveas(fig, filename);
end

% ===================== 函数部分 =========================

% 计算归一化的径向概率分布
function P = RadialProbability(n, l, r_single, r_array)
    R = RadialWavefunction(n, l, r_single);
    P_unnorm = r_single^2 * abs(R)^2;
    
    % 为整个 r_array 计算 R 和未归一化 P，用于积分归一化
    R_full = arrayfun(@(r) RadialWavefunction(n, l, r), r_array);
    P_full = r_array.^2 .* abs(R_full).^2;
    
    normalization_factor = trapz(r_array, P_full);
    P = P_unnorm / normalization_factor;
end

% 计算径向波函数
function R = RadialWavefunction(n, l, r)
    a0 = 52.9;
    rho = (2 * r) / (n * a0);
    normalization = sqrt((2 / n)^3 * factorial(n - l - 1) / (2 * n * factorial(n + l)));
    Lag = LaguerreL(n - l - 1, 2 * l + 1, rho);
    R = normalization * exp(-rho / 2) * (rho^l) * Lag;
end

% Laguerre 多项式计算
function L = LaguerreL(n, alpha, x)
    if n == 0
        L = 1;
    elseif n == 1
        L = 1 + alpha - x;
    else
        L0 = 1;
        L1 = 1 + alpha - x;
        for k = 2:n
            L = ((2*k + alpha - 1 - x) * L1 - (k + alpha - 1) * L0) / k;
            L0 = L1;
            L1 = L;
        end
    end
end