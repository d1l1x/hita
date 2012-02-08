function [eta,u_eta,tau]=KolmoScale(nu,Dissipation)
    eta = (nu^3/Dissipation)^(1/4);
    u_eta = (nu*Dissipation)^(1/4);
    tau = (nu/Dissipation)^(1/2);
end