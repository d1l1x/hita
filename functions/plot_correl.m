function plot_correl(r,R11,R22)
close all
h=plot(r,R11,'r-s');hold on;
set(h,'LineWidth',1);
h=plot(r,R22,'*-g');
set(h,'LineWidth',1);

saveas(gcf,'./doc/correlation.eps','psc2')
end