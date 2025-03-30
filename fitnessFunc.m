function fitness = fitnessFunc(x, materials)
    % x: 决策变量，各层材料的厚度 [气凝胶, 岩棉, 铝箔]
    
    % 更新材料厚度
    for i = 1:length(materials)-1
        materials(i).thickness = x(i);
    end

    % 计算温度
    [Tmax, Tmin] = getTemperatureFunc(materials);

    % 管道参数
    L = 600; % 管道长度 (m)
    R0 = 0.15; % 初始管道半径（内半径）
    total_thickness = sum(x);
    % 计算总成本
    total_cost = 0;
    R_inner = R0; % 当前层的内半径
    for i = 1:length(materials)
        R_outer = R_inner + materials(i).thickness; % 计算外半径
        volume = pi * L * (R_outer^2 - R_inner^2); % 圆环体积
        total_cost = total_cost + volume * materials(i).unit_cost;
        R_inner = R_outer; % 更新下一层的内半径
    end

    % 设定平滑惩罚值
    penalty = 0;
    if Tmax < 539
        penalty = penalty + (540 - Tmax)^2 * 1000; % 平滑二次惩罚
    end
    if Tmin < 399
        penalty = penalty + (400 - Tmin)^2 * 1000; % 平滑二次惩罚
    end
     % 厚度约束惩罚
    if total_thickness > 0.15
        penalty = penalty + (total_thickness - 0.15)^2 * 1e6; % 强烈惩罚
    end
    % 适应度函数（成本 + 约束惩罚）
    fitness = total_cost + penalty;
end
