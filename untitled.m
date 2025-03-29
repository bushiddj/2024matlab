clear;
clc;

% 定义材料的结构体
materials(1).name = '纳米气凝胶层';
materials(1).thermal_conductivity = 0.018; % 热导率K (W/m·K)
materials(1).density = 180;                % 密度ρ (kg/m^3)
materials(1).specific_heat = 1000;         % 比热容C (J/kg·K)
materials(1).thickness = 0.04;             % 厚度(m)
materials(1).thickness_range = [0.02, 0.05]; % 允许厚度范围 (m)
materials(1).unit_cost = 8500;             % 单位体积成本 (元/m^3)

materials(2).name = '岩棉层';
materials(2).thermal_conductivity = 0.042; % 热导率K (W/m·K)
materials(2).density = 120;                % 密度ρ (kg/m^3)
materials(2).specific_heat = 800;          % 比热容C (J/kg·K)
materials(2).thickness = 0.08;             % 厚度(m)
materials(2).thickness_range = [0.04, 0.1];% 允许厚度范围 (m)
materials(2).unit_cost = 2800;             % 单位体积成本 (元/m^3)

materials(3).name = '铝箔反射层';
materials(3).thermal_conductivity = 0.2;   % 热导率K (W/m·K)
materials(3).density = 2700;               % 密度ρ (kg/m^3)
materials(3).specific_heat = 900;          % 比热容C (J/kg·K)
materials(3).thickness = 0.01;             % 厚度(m)
materials(3).thickness_range = [0.005, 0.02]; % 允许厚度范围 (m)
materials(3).unit_cost = 1500;             % 单位体积成本 (元/m^3)

% 定义熔盐的结构体
molten_salt.name = '熔盐';
molten_salt.density = 1800;                % 密度ρ (kg/m^3)
molten_salt.specific_heat = 1500;          % 比热容C (J/kg·K)
molten_salt.thermal_conductivity = 0.55;   % 热导率K (W/m·K)
molten_salt.diameter = 0.3;                % 直径(m)
molten_salt.velocity = 0.5;                % 速度(m/s)

% 参数设置
z_max = 600;  % z 的最大值
r_max = 0.15 + materials(1).thickness + materials(2).thickness + materials(3).thickness; % r 的最大值
Nz = 281*2;     % z 方向的离散点数
Nr = 281*2;     % r 方向的离散点数
dz = z_max / (Nz - 1);
dr = r_max / (Nr - 1);

a1 = molten_salt.thermal_conductivity / (molten_salt.density * molten_salt.specific_heat);
a2 = materials(1).thermal_conductivity / (materials(1).density * materials(1).specific_heat);
a3 = materials(2).thermal_conductivity / (materials(2).density * materials(2).specific_heat);
a4 = materials(3).thermal_conductivity / (materials(3).density * materials(3).specific_heat);

p1 = molten_salt.velocity / (molten_salt.density * molten_salt.specific_heat);
max_dt = 0.5 * min([dr^2 / a1, dr^2 / a2, dr^2 / a3, dr^2 / a4]);
dt = 1.0;    % 时间步长和总时间
hour = 8;
minute = 40;
k0 = molten_salt.thermal_conductivity;
k1 = materials(1).thermal_conductivity;
k2 = materials(2).thermal_conductivity;
k3 = materials(3).thermal_conductivity;

p2 = 10.5;
p3 = temperature_at_time(60*hour+minute);

% 初始化温度场
T = zeros(Nr, Nz);

% 初始条件
r_threshold = 0.15;
for i = 1:Nr
    if (i-1) * dr <= r_threshold
        T(i, :) = 600;
    else
        T(i, :) = temperature_at_time(60*hour+minute);
    end
end

% 这里需要根据实际 r 坐标确定界面：
r_interface1 = 0.15; % 第一区域与第二区域的界面
r_interface2 = r_interface1+materials(1).thickness;  % 第二区域与第三区域界面
r_interface3 = r_interface2+materials(2).thickness;  % 第三区域与第四区域界面
%
% 计算对应的网格索引（假设 r 从0开始，网格步长 dr）
i_interface1 = round(r_interface1/dr) + 1;
i_interface2 = round(r_interface2/dr) + 1;
i_interface3 = round(r_interface3/dr) + 1;

% 计算移动步长
move_step  = round(molten_salt.velocity * dt / dz);

% 初始化存储数据的矩阵
total_hours = 24*3;  % 记录小时的数据
T_history = zeros(total_hours, Nr);  % 存储每小时的温度分布
time_hours = zeros(1, total_hours);  % 存储时间点
hour_count = 1;
last_recorded_hour = -1;

% 时间步进
for t0 = 60*hour+minute:dt:60*hour+minute+60*total_hours
    t=mod(t0,60*24);
    p3 = temperature_at_time(t);
    T_old = T;
    
    % ADI 方法的第一步：隐式处理 r 方向扩散
    for j = 2:Nz-1
        A_r = zeros(Nr, Nr);
        b_r = zeros(Nr, 1);
        for i = 2:Nr-1
            r_pos = (i-1) * dr;
            if r_pos <= r_threshold
                alpha = a1;
            elseif r_pos <= r_threshold + materials(1).thickness
                alpha = a2;
            elseif r_pos <= r_threshold + materials(1).thickness + materials(2).thickness
                alpha = a3;
            else
                alpha = a4;
            end
            
            A_r(i, i-1) = -alpha / dr^2 - alpha / (2 * dr * r_pos);
            A_r(i, i) = 1 / dt + 2 * alpha / dr^2;
            A_r(i, i+1) = -alpha / dr^2 + alpha / (2 * dr * r_pos);
            b_r(i) = T_old(i, j) / dt;
        end
           % --- 修改内部界面处的离散方程 ---
    % 在 r = r_interface1 处（界面1）：两侧分别使用 k0 (左侧) 和 k1 (右侧)
 % 在 r = r_interface1 处（界面1）：两侧分别使用 k0 (左侧) 和 k1 (右侧)
    if i_interface1 > 1 && i_interface1 < Nr
        A_r(i_interface1, :) = 0;
        A_r(i_interface1, i_interface1-1) = -k0 / (dr/2);
        A_r(i_interface1, i_interface1)   = (k0 / (dr/2) + k1 / (dr/2));
        A_r(i_interface1, i_interface1+1) = -k1 / (dr/2);
        b_r(i_interface1) = 0;

        b_r(i_interface1) = 0;
    end
    % 在 r = r_interface2 处（界面2）：两侧分别使用 k1 (左侧) 和 k2 (右侧)
    if i_interface2 > 1 && i_interface2 < Nr
        A_r(i_interface2, :) = 0;
        A_r(i_interface2, i_interface2-1) = -k1 / (dr/2);
        A_r(i_interface2, i_interface2)   = (k1 / (dr/2) + k2 / (dr/2));
        A_r(i_interface2, i_interface2+1) = -k2 / (dr/2);
        b_r(i_interface2) = 0;
    end
    % 在 r = r_interface3 处（界面3）：两侧分别使用 k2 (左侧) 和 k3 (右侧)
    if i_interface3 > 1 && i_interface3 < Nr
        A_r(i_interface3, :) = 0;
        A_r(i_interface3, i_interface3-1) = -k2 / (dr/2);
        A_r(i_interface3, i_interface3)   = (k2 / (dr/2) + k3 / (dr/2));
        A_r(i_interface3, i_interface3+1) = -k3 / (dr/2);
        b_r(i_interface3) = 0;
    end
    % -----------------------------
    % -----------------------------
    
    % r=0处的边界条件（轴对称性）
    A_r(1, 1) = 1;
    A_r(1, 2) = -1;
    b_r(1) = 0;
    % r=r_max 处的边界条件（例如，对流或指定温度）
    A_r(Nr, Nr-1) = -k3 / dr;
    A_r(Nr, Nr) = k3 / dr + p2;
    b_r(Nr) = p2 * p3;
    T(:, j) = A_r \ b_r;
    end
    
    % ADI 方法的第二步：隐式处理 z 方向扩散
    for i = 2:Nr-1
        A_z = zeros(Nz, Nz);
        b_z = zeros(Nz, 1);
        for j = 2:Nz-1
            r_pos = (i-1) * dr;
            if r_pos <= r_threshold
                alpha = a1; conv_term = p1;
            elseif r_pos <= r_threshold + materials(1).thickness
                alpha = a2; conv_term = 0;
            elseif r_pos <= r_threshold + materials(1).thickness + materials(2).thickness
                alpha = a3; conv_term = 0;
            else
                alpha = a4; conv_term = 0;
            end
            
           A_z(j, j-1) = -alpha / dz^2 - max(conv_term, 0) / dz;
           A_z(j, j) = 1 / dt + 2 * alpha / dz^2 + abs(conv_term) / dz;
           A_z(j, j+1) = -alpha / dz^2 - min(conv_term, 0) / dz;
            b_z(j) = T(i, j) / dt;
        end
        
        A_z(1, 1) = 1; b_z(1) = 600;
        A_z(Nz, Nz-1) = -p2;
        A_z(Nz, Nz) = p2 + 1/dt;
        b_z(Nz) = p2 * p3;
  
        
        T(i, :) = (A_z \ b_z)';
    end
    
    % 检查收敛性
    if max(max(abs(T - T_old))) < 1e-6
        fprintf('Solution converged at time t = %.2f\n', t);
        break;
    end
    
    % 更新温度场，模拟x方向移动
    %for j = Nz:-1:2
    %T(:, j) = T(:, j) - (dt / dz) * molten_salt.velocity .* (T(:, j) - T(:, j-1));
    %end
    T(:,1)=T(:,2);
    % 绘制中间结果，每隔一段时间绘制一次
    if mod(t0, 100) < dt % 每隔10个时间单位绘制一次
        surf(linspace(0, z_max, Nz), linspace(0, r_max, Nr), T, 'EdgeColor', 'none');
        colorbar;
        xlabel('z');
        ylabel('r');
        zlabel('Temperature');
        title(['Temperature Distribution at time t = ', num2str(t0)]);
        pause(0.01); % 暂停以便可视化更新
    end
    if floor(t0/60) > last_recorded_hour && hour_count <= total_hours
    % 记录当前小时的温度分布
    T_history(hour_count, :) = T(:,end);
    time_hours(hour_count) = t0/60;  % 转换为小时
    hour_count = hour_count + 1;
    last_recorded_hour = floor(t0/60);
    end

end

% 在主循环结束后，绘制三维图
figure(3);
[R_grid, Time_grid] = meshgrid(linspace(0, r_max, Nr), time_hours);  % 生成温度和时间的网格
surf(R_grid, Time_grid, T_history, 'EdgeColor', 'none');  
xlabel('半径 r (m)');
ylabel('时间 t (小时)');
zlabel('温度分布 (K)');
title(['240小时内温度分布随时间和温度的变化 (z = ', num2str(Nz*dz), ' m)']);
colorbar;
view(45, 30);  % 设置视角


% 在三维图中添加材料界面的平面
hold on;
zLimits = zlim;
tLimits = [min(time_hours) max(time_hours)];
% 添加材料界面的平面
patch([r_interface1 r_interface1 r_interface1 r_interface1], ...
      [tLimits(1) tLimits(1) tLimits(2) tLimits(2)], ...
      [zLimits(1) zLimits(2) zLimits(2) zLimits(1)], ...
      'r', 'FaceAlpha', 0.1);
patch([r_interface2 r_interface2 r_interface2 r_interface2], ...
      [tLimits(1) tLimits(1) tLimits(2) tLimits(2)], ...
      [zLimits(1) zLimits(2) zLimits(2) zLimits(1)], ...
      'r', 'FaceAlpha', 0.1);
patch([r_interface3 r_interface3 r_interface3 r_interface3], ...
      [tLimits(1) tLimits(1) tLimits(2) tLimits(2)], ...
      [zLimits(1) zLimits(2) zLimits(2) zLimits(1)], ...
      'r', 'FaceAlpha', 0.1);

% 添加图例说明
legend('温度分布', '材料界面');


% 取 T_history 的第 i_interface1 列
temperature_data = T_history(:, i_interface1);

% 绘制二维图像
figure;
plot(Time_grid, temperature_data, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
xlabel('时间 (小时)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('温度 (K)', 'FontSize', 12, 'FontWeight', 'bold');
title('温度随时间变化图', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');