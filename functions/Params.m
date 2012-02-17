function [u,v,w,dim,Lx,dx,nu]=Params(uvel,vvel,wvel)
    dim=257; % number of points in one dimension
    Lx=0.4;
%     dim=33;
%     Lx=3.2e-2; % domain size
    Ly=Lx;
    Lz=Lx;
    dx=Lx/(dim-1); % grid spacing
    dy=dx;
    dz=dx;
    nu=1.7e-5; % viscosity
    u=reshape(uvel,dim,dim,dim); % reshape arrays to have them in 3D
    v=reshape(vvel,dim,dim,dim);
    w=reshape(wvel,dim,dim,dim);
    clear uvel vvel wvel
end