function [] = plotAbnormalitySignal(AbnormalSignal,AbnormalSignalSmooth,h)
hold on
xlab1 = xlabel('Time instants ($k$)','interpreter','latex');
xlab1.FontSize = 13;
ylab1 = ylabel('Abnormality measurments','interpreter','latex');
ylab1.FontSize = 14;
plot1 = plot(AbnormalSignal ,'--b','LineWidth',1);
plot2 = plot(AbnormalSignalSmooth ,'-r','LineWidth',2.2);
plot3 = plot([1:size(AbnormalSignal,2)],ones(size(AbnormalSignal,2),1),'--k','LineWidth',2);

h.Position = [1693         571        1663         365];
axis([0  size(AbnormalSignal,2) 0 max(AbnormalSignal)+0.5])
grid on
grid minor

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
rec1 = rectangle('Position',[66 0 (112-66) max(AbnormalSignal)+5],'Curvature',0);
rec1.EdgeColor = 'none';
rec1.FaceColor = [1 0 0 0.099];
rec2 = rectangle('Position',[395 0  (457-395) max(AbnormalSignal)+5],'Curvature',0);
rec2.EdgeColor = 'none';
rec2.FaceColor = [1 0 0 0.099];
hline = line(NaN,NaN,'LineWidth',10,'LineStyle','-','Color',[1 0 0 0.099]);
lg = legend([plot1 plot2 plot3,hline],{'Abnormality Signal',...
    'Smooth Abnormality Signal','Threshold','Abnormal area'},'interpreter','latex');
lg.FontSize = 12;
lg.Position = [0.7888    0.7694    0.1341    0.1274];
end

