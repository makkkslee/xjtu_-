clear all
close all
clc

% 定义主量子数范围
n_values = [1, 2, 3];

% 定义角度范围（极角 theta 从 0 到 pi）
theta = linspace(0, pi, 180);

% 确保 picture/angular_map 文件夹存在
if ~exist('picture', 'dir')
    mkdir picture;
end
if ~exist('picture/angular_map', 'dir')
    mkdir('picture/angular_map');
end

% 循环绘制每个主量子数和角量子数的角向概率分布图
for n = n_values
    % 定义当前主量子数可能的角量子数范围
    l_values = 0:n-1;
    
    for l = l_values
        % 创建新图形窗口
        fig = figure;
        fig.Position = [100, 100, 800, 600];
        fig.Color = [1, 1, 1];
        
        % 定义当前角量子数可能的磁量子数范围
        m_values = -l:l;
        
        % 绘制每个磁量子数的角向概率分布曲线
        hold on;
        for m = m_values
            % 计算角向概率分布
            P = zeros(size(theta));
            for i = 1:length(theta)
                P(i) = AngularProbability(l, m, theta(i));
            end
            
            % 绘制当前磁量子数的角向概率分布曲线
            plot(theta, P, 'LineWidth', 1.5, 'DisplayName', ['m = ', num2str(m)]);
        end
        hold off;
        
        % 设置图形属性
        xlabel('Polar Angle (\theta)', 'FontSize', 12, 'Color', [0, 0, 0]);
        ylabel('Angular Probability Density', 'FontSize', 12, 'Color', [0, 0, 0]);
        title(['Angular Probability Distribution for n = ', num2str(n), ', l = ', num2str(l)], ...
            'Interpreter', 'latex', 'FontSize', 14, 'Color', [0, 0, 0]);
        grid on;
        legend('show');
        axis tight;
        
        % 调整坐标轴范围以确保图像美观
        xlim([0, pi]);
        xticks(linspace(0, pi, 7));
        xticklabels({'0', '\pi/6', '\pi/3', '\pi/2', '2\pi/3', '5\pi/6', '\pi'});
        
        % 保存图像到 picture/angular_map 文件夹
        filename = sprintf('picture/angular_map/Angular_Probability_n%d_l%d.png', n, l);
        saveas(fig, filename);
        
    end
end

% 角向概率分布计算函数
function P = AngularProbability(l, m, theta)
    % 计算球谐函数 Y_l^m(theta, phi)
    Ylm = SphericalHarmonic(l, m, theta);
    
    % 计算角向概率分布 P(theta) = |Y_l^m(theta, phi)|^2
    P = abs(Ylm)^2;
end

% 球谐函数计算函数
function Ylm = SphericalHarmonic(l, m, theta)
    % 检查 l 和 m 的有效性
    if abs(m) > l
        error('m 的绝对值不能大于 l');
    end
    
    % 计算归一化因子
    normalization = sqrt((2*l + 1) * factorial(l - abs(m)) / (4 * pi * factorial(l + abs(m))));
    
    % 计算连带Legendre多项式
    Plm = Legendre(l, abs(m), cos(theta));
    
    % 计算球谐函数
    Ylm = normalization * Plm * exp(1i * m * 0);  % phi=0
end

% 连带Legendre多项式计算函数
function Plm = Legendre(l, m, x)
    % 使用MATLAB内置的legendre函数
    Plm = legendre(l, x);
    
    % 确保索引有效
    if l == 0
        Plm = Plm(1);
    else
        % MATLAB的legendre函数返回一个矩阵，其中每一行对应不同的m值（从0到l）
        % 为了得到m的值从-l到l，需要进行索引转换
        % 这里我们只关心m的绝对值，因为Y_l^m和Y_l^-m的角向分布相同
        % 所以我们取绝对值m，并相应地调整索引
        Plm = Plm(m + 1, :);  % m从0到l，对应索引从1到l+1
    end
end