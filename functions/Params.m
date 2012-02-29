function [u,v,w,Lx,nu]=Params(uvel,vvel,wvel,dim)
    Lx=2*pi; %edge length
    nu=1.7e-5; % viscosity
    u=reshape(uvel,dim,dim,dim); % reshape to 3D
    v=reshape(vvel,dim,dim,dim);
    w=reshape(wvel,dim,dim,dim);
    clear uvel vvel wvel % save memory
end
