function fitness = fitnessFunc(thicknesses, materials, getTemperatureFunc, design_targets)

    % 更新材料厚度
    for i = 1:length(materials)
        materials(i).thickness = thicknesses(i);
    end

    % 计算材料成本
    total_cost = 0;
    for i = 1:length(materials)
        volume = materials(i).thickness * pi * (0.15 + sum(thicknesses(1:i-1)))^2;
        total_cost = total_cost + volume * materials(i).unit_cost;
    end

    % 获取白天和夜间的管道末端熔盐温度
    [day_temp, night_temp] = getTemperatureFunc(materials);

    % 计算温度惩罚
    day_penalty = max(0, design_targets.day - day_temp);
    night_penalty = max(0, design_targets.night - night_temp);
    temperature_penalty = day_penalty + night_penalty;

    % 总惩罚
    fitness = total_cost + 1e6 * temperature_penalty; % 1e6是一个大权重，确保温度目标优先
end