close all
clear all
clc

flag='3D';

if (strcmp('3D',flag))
    uvel=importdata('OUTPUT/3D/uvel');
    vvel=importdata('OUTPUT/3D/vvel');
    wvel=importdata('OUTPUT/3D/wvel');
else
    uvel=importdata('OUTPUT/2D/uvel');
    vvel=importdata('OUTPUT/2D/vvel');
end
%% compute transformations
if (strcmp('3D',flag))
    dim=256;
    Lx=5e-3;
    Ly=Lx;
    Lz=Lx;
    dx=Lx/dim;
    dy=dx;
    dz=dx;
    u=reshape(uvel,dim,dim,dim);
    v=reshape(vvel,dim,dim,dim);
    w=reshape(wvel,dim,dim,dim);
else
    dim=1024;
    Lx=1E-2;
    Ly=Lx;
    dx=Lx/dim;
    dy=dx;
    u=reshape(uvel,dim,dim);
    v=reshape(vvel,dim,dim);
end
%% compute fft
if (strcmp('3D',flag))
    NFFT = 2.^nextpow2(size(u));
    u_fft=fftn(u,NFFT)/dim^3;
    u_fft = u_fft(1:NFFT(1)/2+1,1:NFFT(2)/2+1,1:NFFT(2)/2+1);
    
    NFFT = 2.^nextpow2(size(v));
    v_fft=fftn(v,NFFT)/dim^3;
    v_fft = v_fft(1:NFFT(1)/2+1,1:NFFT(2)/2+1,1:NFFT(2)/2+1);
    
    NFFT = 2.^nextpow2(size(w));
    w_fft=fftn(w,NFFT)/dim^3;
    w_fft = w_fft(1:NFFT(1)/2+1,1:NFFT(2)/2+1,1:NFFT(2)/2+1);
    
    phi_x=u_fft.*conj(u_fft);
    phi_y=v_fft.*conj(v_fft);
    phi_z=w_fft.*conj(w_fft);
else
    NFFT = 2.^nextpow2(size(u));
    u_fft=fftn(u,NFFT)/dim^2;
    u_fft = u_fft(1:NFFT(1)/2+1,1:NFFT(2)/2+1);
    NFFT = 2.^nextpow2(size(u));
    v_fft=fftn(v,NFFT)/dim^2;
    v_fft = v_fft(1:NFFT(1)/2+1,1:NFFT(2)/2+1);
       
    phi_x=u_fft.*conj(u_fft);
    phi_y=v_fft.*conj(v_fft);
end

%% compute correlations

if (strcmp('3D',flag))
    R11=ifftn(fftn(u).*conj(fftn(u))/dim^3/std2(u)^2);
    R22=ifftn(fftn(v).*conj(fftn(v))/dim^3/std2(v)^2);
    R33=ifftn(fftn(w).*conj(fftn(w))/dim^3/std2(w)^2);
    R11=R11(1:round(size(R11,1)/2),1,1);
    R22=R22(1:round(size(R22,1)/2),1,1);
    R33=R33(1:round(size(R33,1)/2),1,1);
    r = linspace(0,Lx/2,dim/2)/(Lx/2);
else
    R11= ifftn(fftn(u).*conj(fftn(u))/dim^2/std2(u)^2);
    R22=ifftn(fftn(v).*conj(fftn(v))/dim^2/std2(v)^2);
    R11=R11(1:round(size(R11,1)/2),1);
    R22=R22(1:round(size(R22,1)/2),1);
    r = linspace(0,Lx/2,dim/2)/(Lx/2);
    test = R11 + r'/2.*gradient(R11,1/256);
    plot(r,R11,r,R22,r,test)
    h=line([0 1],[0 0],'Color',[0 0 0],'LineWidth',1.0);
end
%% compute length scales
L11=trapz(r,R11);
L22=trapz(r,r'.*R22);
hold on
rectangle('Position',[0,0,L11,1],'LineWidth',2,'LineStyle','--')

%% compute 1D spectrum
L=length(R11);
NFFT=2^nextpow2(L);
spec=fft(R11,NFFT)/L;
f = L/2*linspace(0,1,NFFT/2+1);

slope=1.5*664092^(2/3)*(f.^(-5/3));
loglog(f,2*abs(spec(1:NFFT/2+1)));
% hold on
loglog(f,slope);
%% compute spectrum
% spec = zeros(round(dim*dim*dim/8),1);
if (strcmp('3D',flag))
%     phi = u_fft;
%     phi(:,:,:)=0.0;
%     for k=1:dim
%         for j=1:dim
%             for i=1:dim
    %             kappa = sqrt(i*i+j*j+k*k);
    %             kappa_pos=int16(kappa);
    %             if (kappa_pos <= size(spec,1))
    %                 spec(kappa_pos) = spec(kappa_pos)+kappa*kappa*(...
    %                 + real(u_fft(i,j,k))*real(u_fft(i,j,k))+imag(u_fft(i,j,k))*imag(u_fft(i,j,k)) ...
    %                 + real(v_fft(i,j,k))*real(v_fft(i,j,k))+imag(v_fft(i,j,k))*imag(v_fft(i,j,k)) ...
    %                 + real(w_fft(i,j,k))*real(w_fft(i,j,k))+imag(w_fft(i,j,k))*imag(w_fft(i,j,k)));

    %             end
    %             spec(kappa_pos) = spec(kappa_pos) + kappa*kappa*0.5*(phi_x(i,j,k).^+phi_y(i,j,k).^2+phi_z(i,j,k).^2);
                  phi = 0.5*(phi_x+phi_y+phi_z);
%             end
%         end
%     end
else
%     phi = u_fft;
%     phi(:,:)=0.0;
%     for j=1:dim
%         for i=1:dim
%             phi(i,j) = phi(i,j) +(phi_x(i,j)+phi_y(i,j));
%         end
%     end
    phi = 0.5*phi_x+phi_y;
end
%% compute k vector
if (strcmp('3D',flag))
maxdim = sqrt(dim^2*(2*pi/Lx)^2+dim^2*(2*pi/Ly)^2+dim^2*(2*pi/Lz)^2);
    E=zeros(uint64(sqrt(3*dim^2)),1);
    kk=zeros(uint64(sqrt(3*dim^2)),1);
    dim = size(phi,1);
    for k=1:dim
        for j=1:dim
            for i=1:dim
                kappa=sqrt(i*i*(2*pi/Lx)^2+j*j*(2*pi/Ly)^2+k*k*(2*pi/Lz)^2);
                kappa_pos=uint64(sqrt(i*i+j*j+k*k));    
                E(kappa_pos) = E(kappa_pos) + phi(i,j,k);
                kk(kappa_pos) = kappa;
            end
        end
    end
else
    maxdim = sqrt(dim^2*(2*pi/Lx)^2+dim^2*(2*pi/Ly)^2);
    E=zeros(uint64(sqrt(dim^2+dim^2)),1);
    kk=zeros(uint64(sqrt(dim^2+dim^2)),1);
    dim = size(phi,1);
    for j=1:dim
        for i=1:dim
            kappa=sqrt(i*i*(2*pi/Lx)^2+j*j*(2*pi/Ly)^2);
            kappa_pos=uint64(sqrt(i*i+j*j));
            E(kappa_pos) = E(kappa_pos) + phi(i,j);
            kk(kappa_pos) = kappa;
        end
    end
end
%% compute 1D spectrumcl
close all
% E=E*4*pi;
% E=E.*kk.^2;
slope=1.5*664092^(2/3)*(kk.^(-5/3));
test=importdata('INPUT/2D/CTRL_TURB_ENERGY');
%
% dissip=664092;
dissip=693207;
up=17;
L=Lx;
kkke=kk./(2*pi)*L;
kkkd=kk./(2*pi*100)*L;
VKP = 1.5*17^5/dissip.*(kkke).^4./(1+kkke.^2).^(17/6).*exp(-3/2*1.5.*(kkkd).^(4/3));
%
loglog(test(:,1),test(:,3),kk,slope,kk,VKP,kk(2:end),E(2:end))
h=legend('Initial','Kolmogorov','VKP','Computed');
set(h,'Location','SouthWest')
