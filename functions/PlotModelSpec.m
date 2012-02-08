function PlotModelSpec(k,spectrum,Dissipation,up,Lx,eta,nu)
    % Von Karman-Pao Spektren
    close all
    kd = k(end);
    ke = pi/Lx/2;
    A = 1.5;
    up = mean2(up);

    VKP1 = A*up^5/Dissipation.*(k./ke).^4./(1+(k./ke).^2).^(17/6).*exp(-3/2*A.*(k./kd).^(4/3));

    kd = 1./eta;
    VKP2 = 1.5*(k./kd).^(-5/3)./(Dissipation*nu^5)^(-1/4).*exp(-1.5*1.5.*(k./kd).^(4/3));

    % Kolmogorov Spektrum
    Kolmo=1.5*Dissipation^(2/3)*(k.^(-5/3));

    % Plot spectra
    h=loglog(k,Kolmo,k,VKP1,k,VKP2,k,spectrum);
    set(h,'LineWidth',2);

    h=legend('Kolmogorov','VKP1','VKP2','Computed');
    set(h,'Location','SouthWest')
end
