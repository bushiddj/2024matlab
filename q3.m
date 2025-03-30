%问题2
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

materials(4).name = '沙土';
materials(4).thermal_conductivity = 0.3;   % 热导率K (W/m·K)
materials(4).density = 1800;               % 密度ρ (kg/m^3)
materials(4).specific_heat = 1000;          % 比热容C (J/kg·K)
materials(4).thickness = 0;             % 厚度(m)
materials(4).thickness_range = [0.00, 0.02]; % 允许厚度范围 (m)
% 定义熔盐的结构体
molten_salt.name = '熔盐';
molten_salt.density = 1800;                % 密度ρ (kg/m^3)
molten_salt.specific_heat = 1500;          % 比热容C (J/kg·K)
molten_salt.thermal_conductivity = 0.55;   % 热导率K (W/m·K)
molten_salt.diameter = 0.3;                % 直径(m)
molten_salt.velocity = 0.5;                % 速度(m/s)
[Tmax2, Tmin2] = getTemperatureFunc(materials)