function [fig_hdl,sb1_hdl,sb2_hdl] = topostats_PlotStats(stats)
% topostats_PlotStats(stats) plots a figure with the results provided by
%  any of the topostats_* analyses (TANOVA, GFP, TCT)

fig_hdl = figure;

sb1_hdl = subplot(3,1,[1 2]);
imagesc(stats.time,[0 max(stats.stat)*1.2],~stats.mask,[-1 1]); hold on;
colormap(gca,gray(4));
set(gca,'YDir','normal');
hold on;
plot(stats.time,stats.stat,'k','LineWidth',1);
grid on;
ylim([0 max(stats.stat)*1.2]);
xlabel('Time');
ylabel(strrep(stats.label{1},'topostats_',''));


sb2_hdl = subplot(3,1,3);
imagesc(stats.time,[0 1],~stats.mask,[-1 1]); hold on;
set(gca,'YDir','normal');
ylim([0 1]);
colormap(gca,gray(4));
area(stats.time,stats.prob,'FaceColor',[0 0 0]);
yline(0.05,'r');
xlabel('Time');
ylabel('p-value');
grid on;

end