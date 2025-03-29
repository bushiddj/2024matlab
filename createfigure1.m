function createfigure1(X1, YMatrix1)
%CREATEFIGURE1(X1, YMatrix1)
%  X1:  x 数据的向量
%  YMATRIX1:  y 数据的矩阵

%  由 MATLAB 于 29-Mar-2025 21:41:33 自动生成 question1

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot(X1,YMatrix1,'MarkerFaceColor',[0 0 1],'Marker','o','LineWidth',1.5,...
    'Parent',axes1);

% 创建 ylabel
ylabel('温度 (K)','FontWeight','bold');

% 创建 xlabel
xlabel('时间 (小时)','FontWeight','bold');

% 创建 title
title('温度随时间变化图','FontWeight','bold');

box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% 设置其余坐标区属性
set(axes1,'FontSize',12,'FontWeight','bold');
