function createfigure(X1, Y1, S1, C1, Y2, C2)
%CREATEFIGURE(X1, Y1, S1, C1, Y2, C2)
%  X1:  scatter x
%  Y1:  scatter y
%  S1:  scatter s
%  C1:  scatter c
%  Y2:  scatter y
%  C2:  scatter c

%  Auto-generated by MATLAB on 09-Aug-2019 01:31:45

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create scatter
scatter(X1,Y1,S1,C1,'DisplayName','umax','MarkerFaceColor',[0 0 0],...
    'MarkerEdgeColor',[0.635294139385223 0.0784313753247261 0.184313729405403],...
    'Marker','square');

% Create scatter
scatter(X1,Y2,S1,C2,'DisplayName','target','MarkerEdgeColor',[0 0 0]);

% Create xlabel
xlabel('C:N ratios');

% Create ylabel
ylabel('umax (/h) and target production (mmol/gDW/h)');

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0 50]);
%box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontName','Arial','FontSize',12,'XColor',[0 0 0],'XMinorTick',...
    'on','YColor',[0 0 0],'YMinorTick','on','ZColor',[0 0 0]);
% Create legend
legend(axes1,'show');

