function T = temperature_at_time(t)
    
    % 判断当前时间是白天还是其他时间
    if t >= 6*60 && t <= 16*60
        % 白天时间：6:00 到 16:00
        beta = 6; % 题目中给出的beta值
        T = 48 - beta * abs(t/60 - 11);
    else
        % 其他时间
        T = -15;
    end
end