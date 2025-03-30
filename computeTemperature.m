
function T_history = computeTemperature(materials)
    % 初始化参数（从原代码提取关键参数）
    Nz = 100;
    Nr = 100;
    z_max = 600;
    r_max = 0.15 + materials(1).thickness + materials(2).thickness + materials(3).thickness + 0.02;
    dz = z_max / (Nz - 1);
    dr = r_max / (Nr - 1);
    dt = 1.2;
    total_hours = 0.5;
    record_interval = 3;
    
    % 计算所需时间步数
    num_records = ceil(total_hours * 60 / record_interval);
    T_history = zeros(num_records, Nr);
    
    % 初始化温度场
    T = zeros(Nr, Nz);
    T(:, :) = 600;  % 假设初始温度
    
    % 迭代计算温度分布（简化计算逻辑）
    for t_step = 1:num_records
        % 这里可以插入完整的 ADI 方法或其他计算方法
        % 这里只是简单模拟温度变化
        T = T - 0.1 * rand(Nr, Nz);  % 假设温度随时间降低
        T_history(t_step, :) = T(:, end);  % 记录末端温度
    end
end
