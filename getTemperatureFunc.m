function [Tmax, Tmin] = getTemperatureFunc(materials)
% 定义熔盐的结构体
molten_salt.name = '熔盐';
molten_salt.density = 1800;                % 密度ρ (kg/m^3)
molten_salt.specific_heat = 1500;          % 比热容C (J/kg·K)
molten_salt.thermal_conductivity = 0.55;   % 热导率K (W/m·K)
molten_salt.diameter = 0.3;                % 直径(m)
molten_salt.velocity = 0.5;                % 速度(m/s)

% 参数设置
z_max = 600;  % z 的最大值
r_max = 0.15 + materials(1).thickness + materials(2).thickness + materials(3).thickness+0.02; % r 的最大值
Nz = 100;     % z 方向的离散点数
Nr = 100;     % r 方向的离散点数
dz = z_max / (Nz - 1);
dr = r_max / (Nr - 1);

a1 = molten_salt.thermal_conductivity / (molten_salt.density * molten_salt.specific_heat);
a2 = materials(1).thermal_conductivity / (materials(1).density * materials(1).specific_heat);
a3 = materials(2).thermal_conductivity / (materials(2).density * materials(2).specific_heat);
a4 = materials(3).thermal_conductivity / (materials(3).density * materials(3).specific_heat);
a5 = materials(4).thermal_conductivity / (materials(4).density * materials(4).specific_heat);

p1 = molten_salt.velocity / (molten_salt.density * molten_salt.specific_heat);
max_dt = 0.5 * min([dr^2 / a1, dr^2 / a2, dr^2 / a3, dr^2 / a4]);
dt = 1.2;    % 时间步长和总时间
hour = 5;
minute = 30;
record_interval=3;
k0 = molten_salt.thermal_conductivity;
k1 = materials(1).thermal_conductivity;
k2 = materials(2).thermal_conductivity;
k3 = materials(3).thermal_conductivity;
k4 = materials(4).thermal_conductivity;

p2 = 10.5;
p3 = temperature_at_time(60*hour+minute);
r_salt_max=0.02;
v_threshold=8;
isnight=0;
wind_speeds = [20, 15, 4.0, 4.0, 4.0, 4.0, 0, 0, 0, 0, 0, 0, ...
               0, 0, 0,0, 6.0, 6.0, 6.0, 6.0,6.0, 6.0, 8.0, 8.0];
m3_absorb=0.6;
m3_ref=0.9;
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
r_interface4 = r_interface3+materials(3).thickness;  % 第四区域与第五区域界面
r_interface5 = r_interface4+materials(4).thickness;
%
% 计算对应的网格索引（假设 r 从0开始，网格步长 dr）
i_interface1 = round(r_interface1/dr) + 1;
i_interface2 = round(r_interface2/dr) + 1;
i_interface3 = round(r_interface3/dr) + 1;
i_interface4 = round(r_interface4/dr) + 1;
i_interface5 = round(r_interface5/dr) + 1;
% 计算移动步长
move_step  = round(molten_salt.velocity * dt / dz);

% 初始化存储数据的矩阵
total_hours = 0.5;  % 记录小时的数据
T_history = zeros(ceil(total_hours * 60 / record_interval), Nr);  
time_records = zeros(1, ceil(total_hours * 60 / record_interval));  % 存储时间点
record_count = 1;
last_recorded_time = -record_interval;  % 初始化为负的记录间隔
% 时间步进
for t0 = 60*hour+minute:dt:60*hour+minute+60*total_hours
    t=mod(t0,60*24);
    minute =mod(t,60);
    hour = floor(t/60);
    p3 = temperature_at_time(t);
    T_old = T;
    wind_speed=wind_speeds(hour + 1);

    if (t>9*60+13)&&(t<12*60+47)
    materials(1).thermal_conductivity = 0.018+0.004*(G_t(t)-800)/(950-800);
    else
        materials(1).thermal_conductivity =0.018;
    end

    materials(4).thickness=r_salt_max*(1-exp(-wind_speed/v_threshold));
    r_interface5 = r_interface4 +materials(4).thickness;
    i_interface5 = round(r_interface5 / dr) + 1;
    % 动态调整 Nz 并插值温度数据
    Nz_new = round(z_max / dz) + 1;
    if Nz_new ~= Nz
        T_interp = zeros(Nr, Nz_new);
        z_old = linspace(0, z_max, Nz);
        z_new = linspace(0, z_max, Nz_new);
        for i = 1:Nr
            T_interp(i, :) = interp1(z_old, T(i, :), z_new, 'linear');
        end
        T = T_interp;
        Nz = Nz_new;
    end

    for i = 1:Nr
        r_pos = (i-1) * dr;
        if r_pos > r_interface5
            % 如果网格点位于沙土层之外，将温度设置为环境温度
            T(i, :) = temperature_at_time(t0);
        end
    end

    if (t < 6*60) || (t > 16)
    materials(2).thermal_conductivity = 0.042+0.008;
    isnight=1;
    else
    materials(2).thermal_conductivity = 0.042;
    isnight=0;
    end
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
            elseif r_pos <= r_threshold + materials(1).thickness + materials(2).thickness+ materials(3).thickness
                alpha = a4;
            else
                alpha = a5;
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
    if i_interface4 > 1 && i_interface4 < Nr
        A_r(i_interface4, :) = 0;
        A_r(i_interface4, i_interface4-1) = -k3 / (dr/2);
        A_r(i_interface4, i_interface4)   = (k3 / (dr/2) + k4 / (dr/2));
        A_r(i_interface4, i_interface4+1) = -k4 / (dr/2);
        b_r(i_interface4) = 0;
    end
    % -----------------------------
    % -----------------------------
    
    % r=0处的边界条件（轴对称性）
    A_r(1, 1) = 1;
    A_r(1, 2) = -1;
    b_r(1) = 0;
    % r=r_max 处的边界条件（例如，对流或指定温度）
        if(isnight==1)
            A_r(Nr, Nr-1) = -k4 / dr;
            A_r(Nr, Nr) = (k4 / dr + p2);
            b_r(Nr) = p2 * p3;
            T(:, j) = A_r \ b_r;
        else
            A_r(Nr, Nr-1) = -k3 / dr;
            A_r(Nr, Nr) = k3 / dr + p2;
            b_r(Nr) = p2 * p3+(m3_ref-m3_absorb)*G_t(t);
            T(:, j) = A_r \ b_r;
        end
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
            elseif r_pos <= r_threshold + materials(1).thickness + materials(2).thickness+ materials(3).thickness
                alpha = a4; conv_term = 0;
            else
                alpha = a5; conv_term = 0;
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
    if(hour>=5)&&(hour<=6)&&(minute>=45)||(hour>=6)&&(hour<=7)&&(minute<=15)
        if t0 - last_recorded_time >= record_interval
            T_history(record_count, :) = T(:, end);
            time_records(record_count) = t0;  % 保留时间为分钟
            record_count = record_count + 1;
            last_recorded_time = t0;
        end
    end
end
T1=T_history(1 :5,i_interface1)
Tmax=max(T1)
Tmin=min(T1)
end