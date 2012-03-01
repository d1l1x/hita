function h=plot_spec(k,spectrum)
display('Save results...')
close all
h = figure('visible','on');% u=u-mean2(u);
% v=v-mean2(v);
% w=w-mean2(w);
data=importdata('spectrumSpektralCode.txt');
kC=data(:,1);
EC=data(:,2);
% vkp=importdata('data/3D/CTRL_TURB_ENERGY');
% h=loglog(vkp(:,1),vkp(:,3),'*-b');hold on
% set(h,'LineWidth',1);
h=axes;
set(h,'Fontsize',14)
set(h,'Xscale','log');
set(h,'Yscale','log');
set(h,'Linewidth',2);
set(h,'box','on');
hold on
h=loglog(k,spectrum,'r-s');
set(h,'LineWidth',1);
% h=loglog(kC,EC,'*-g');
% set(h,'LineWidth',1);
% legend('VKP','Dietzsch','Comte-Bellot')


saveas(gcf,'./doc/spectrum.eps','psc2')
end