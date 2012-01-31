function [Dissipation,kinetic_E] = spec_prop(E,k,nu)
    kinetic_E = trapz(k,E);
    Dissipation = trapz(k,2*nu.*k.^2.*E');
end