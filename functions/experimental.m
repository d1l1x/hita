function kin=experimental(k,spectrum)
Imax = size(k,2);
kin=0;
dk=1;
%for i=1:Imax-1
%     kin=kin+(spectrum(i)+spectrum(i+1))/2*dk;
% end
Imax=63;
for i=1:Imax-1
    kin=kin+spectrum(1+i*(size(spectrum,1)-1)/Imax);
end
kin=(kin+(spectrum(1)+spectrum(end))/2)*(k(end)-k(1))/Imax;