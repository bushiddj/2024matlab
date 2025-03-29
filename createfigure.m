function createfigure(xdata1, ydata1, zdata1, ZData2, YData2, XData2, XData3, XData4)
%CREATEFIGURE(xdata1, ydata1, zdata1, ZData2, YData2, XData2, XData3, XData4)
%  XDATA1:  surface xdata
%  YDATA1:  surface ydata
%  ZDATA1:  surface zdata
%  ZDATA2:  patch zdata
%  YDATA2:  patch ydata
%  XDATA2:  patch xdata
%  XDATA3:  patch xdata
%  XDATA4:  patch xdata

%  由 MATLAB 于 29-Mar-2025 21:41:11 自动生成question1

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 创建 surf
surf(xdata1,ydata1,zdata1,'DisplayName','温度分布','Parent',axes1,...
    'EdgeColor','none');

% 创建 patch
patch('DisplayName','材料界面','Parent',axes1,'ZData',ZData2,'YData',YData2,...
    'XData',XData2,...
    'FaceAlpha',0.1,...
    'FaceColor',[1 0 0]);

% 创建 patch
patch('Parent',axes1,'ZData',ZData2,'YData',YData2,'XData',XData3,...
    'FaceAlpha',0.1,...
    'FaceColor',[1 0 0]);

% 创建 patch
patch('Parent',axes1,'ZData',ZData2,'YData',YData2,'XData',XData4,...
    'FaceAlpha',0.1,...
    'FaceColor',[1 0 0]);

% 创建 zlabel
zlabel('温度分布 (K)');

% 创建 ylabel
ylabel('时间 t (小时)');

% 创建 xlabel
xlabel('半径 r (m)');

% 创建 title
title('168小时内温度分布随时间和温度的变化');

view(axes1,[112.18775376714 19.4214954222596]);
grid(axes1,'on');
hold(axes1,'off');
% 创建 colorbar
colorbar(axes1);

% 创建 legend
legend(axes1,'show');

