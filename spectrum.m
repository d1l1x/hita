%% Small DNS postprocessor
% INTRODUCTORY TEXT

%% Clear complete workspace
% Its always a good idea to clear the complete workspace and the command
% window also closing all figures might be helpful. You may also use the
% header defin some neccessary flags distinguishing bewteen different data
% sets.
close all
clear all
clc

flag='3D';
datadir='data';

%% Read data files
% Read in the data files and measure the time for reading. The output of
% the tic/toc block is in seconds. What you should get from the tic/toc
% block is that most of the time is spend during data I/O. The actual
% computation needs only ??? of the time of the I/O operations.

%%% 3D
if (strcmp('3D',flag))
    tic; % enable timer
    uvel=importdata([datadir,'/',flag,'/uvel']);
    vvel=importdata([datadir,'/',flag,'/vvel']);
    wvel=importdata([datadir,'/',flag,'/wvel']);
    time_reading = toc; % end timer
%     test=importdata('data/test/CFX_velocity_field.dat');
%     uvel=test(:,1);
%     vvel=test(:,2);
%     wvel=test(:,3);
end
%%% 2D
if (strcmp('2D',flag))
    tic;
    uvel=importdata([datadir,'/',flag,'/uvel']);
    vvel=importdata([datadir,'/',flag,'/vvel']);
    time_reading = toc;
end
%% Set some neccessary parameters
% For further computations it is important to define some parmeters of the
% DNS simulation such as domain size, grid spacing, and the number of grid
% points.
%%% 3D
if (strcmp('3D',flag))
    dim=256; % number of points in one dimension
    Lx=5e-3; % domain size
%     dim=33;
%     Lx=3.2e-2; % domain size
    Ly=Lx;
    Lz=Lx;
    dx=Lx/dim; % grid spacing
    dy=dx;
    dz=dx;
    nu=1.7e-5; % viscosity
    u=reshape(uvel,dim,dim,dim); % reshape arrays to have them in 3D
    v=reshape(vvel,dim,dim,dim);
    w=reshape(wvel,dim,dim,dim);
end
%%% 2D
if (strcmp('2D',flag))
    dim=1024; % number of points in one dimension
    Lx=1E-2;  % domain size
    Ly=Lx;
    dx=Lx/dim; % grid spacing
    dy=dx;
    u=reshape(uvel,dim,dim); % reshape arrays to have them in 2D
    v=reshape(vvel,dim,dim);
end
%% Compute FFT
% This is the most important part of the script. Since the performance of
% an actual DFT is rather bad the preferred choice is a FFT. The FFT
% approach is fastest if the data set to be transformed has a size that is
% a multiple of two. Thats why the function *nextpow2* is used to get the
% next powert of two approximating the dimension _dim_ of the data set. As
% a consequence the data set is zero padded or truncated. _Since the output
% of an FFT operation is symmetric we only need to save half the transform_.
%%
% 
% <latex>
%   \begin{equation}
%      \Phi_{ij}(\kappa)=\frac{1}{(2\,\pi)^3}\iiint\limits^{\infty}_{-\infty}
%                           R_{ij}(\mathbf{r})\,\mathrm{e}^{-i\mathrm{\kappa}r}
%                           \,\mathrm{d}\mathbf{r}
%   \end{equation}
% After the transformation of all velocity components we have to compute
% the velocity correlation tensor $\Phi$ . From theory we know
%   \begin{equation}
%   (u_i*u_j)=\int\limits_{-\infty}^{\infty}u_i^{*}(\mathbf{x})\,
%                   u_j(\mathbf{x}+\mathbf{r})\,\mathrm{d}\mathbf{r}.
%   \end{equation}
% Since all our data sets are transformed (and we are in the Fourier space)
% the last expression can be simply computed by multiplying
%   \begin{equation}
%       \mathfrak{F}\left\{u_i*u_j\right\} = \alpha\cdot 
%           \left\{\mathfrak{F}\left\{u_i\right\}\right\}^{*}\cdot
%                  \mathfrak{F}\left\{u_j\right\},
%   \end{equation}
% where $\alpha$ is a normalization factor.
% </latex>
scaling = 1;
extension = 0;
if (strcmp('3D',flag))
    tic; % start timer
    NFFT = 2.^nextpow2(size(u)+extension); % next power of 2 fitting the length of u
    u_fft=fftn(u,NFFT)./scaling; %2 pi --> definition of FFT 
    %
    NFFT = 2.^nextpow2(size(v)+extension);
    v_fft=fftn(v,NFFT)./scaling;
    %
    NFFT = 2.^nextpow2(size(w)+extension);
    w_fft=fftn(w,NFFT)./scaling;
    time_fft=toc; % get final time for all transformations
    
    Rij_x=(u_fft.*conj(u_fft)); % compute velo. correlation tensor
    Rij_y=(v_fft.*conj(v_fft));
    Rij_z=(w_fft.*conj(w_fft));
    
end
if (strcmp('2D',flag))
    scaling = 1;
    extension = 0;
    tic; %start timer
    NFFT = 2.^(nextpow2(size(u))+extension);
    u_fft=fft2(u,NFFT(1),NFFT(2))./scaling; %2 pi --> definition of FFT 
    %
    NFFT = 2.^(nextpow2(size(v))+extension);
    v_fft=fft2(v,NFFT(1),NFFT(2))./scaling;
    %
    Rij_x=(u_fft.*conj(u_fft));%/size(u,1).^2/size(u,2).^2;
    Rij_y=(v_fft.*conj(v_fft));%/size(v,1).^2/size(v,2).^2;
    
end

%% Compute correlations
% Computing a correlation can be a tedious work (requireing tremendeous
% effort) especially if you have large data sets. From theory it is well
% known that the multiplication of the transform of a data set and its
% complex conjugate are an accurate representation of the correlation
% function. Using the FFT approach this gives an enormeous speed advantage.
% Since we already computed the veloity correlation tensor we may use this
% result in order to compute the correlation tensor.
%%
% 
% <latex>
%   \begin{equation}
%       R_{ij} = \frac{cov(U_i,U_j)}{\sqrt{\sigma_i^2\,\sigma_j^2}}
%              = \frac{(u_i'-\mu_i)\,(u_j-\mu_j)}{\sqrt{\sigma_i^2\,\sigma_j^2}}
%   \end{equation}
% </latex>
% 
if (strcmp('3D',flag))
    NFFT = 2.^nextpow2(size(u_fft));
    R1=ifftn(Rij_x,NFFT)/NFFT(1)/NFFT(2)/NFFT(3)/std2(u)^2;
    NFFT = 2.^nextpow2(size(v_fft));
    R2=ifftn(Rij_y,NFFT)/NFFT(1)/NFFT(2)/NFFT(3)/std2(v)^2;
    NFFT = 2.^nextpow2(size(w_fft));
    R3=ifftn(Rij_z,NFFT)/NFFT(1)/NFFT(2)/NFFT(3)/std2(w)^2;
    %
%     R33 = reshape(R3(1,1,:),NFFT(1),1);
%     R22 = reshape(R2(1,:,1),NFFT(1),1);
%     R11 = reshape(R1(:,1,1),NFFT(1),1);
    
    R11 = (reshape(R3(1,1,:),NFFT(1),1)+R2(1,:,1)'+R1(:,1,1))/3;
    R11 = R11(1:size(u_fft)/2+1);
    %
    R1_22 = (R1(1,:,1)+R3(1,:,1))/2;
    R2_22 = (R2(:,1,1)+R3(:,1,1))/2;
    R3_22 = (reshape(R1(1,1,:),size(u_fft,1),1)+reshape(R2(1,1,:),size(u_fft,1),1))/2;
    
    R22 = (R1_22'+R2_22+R3_22)/3;
    R22 = R22(1:size(u_fft)/2+1);
%     R22 = (R11(1,:,1) + )/3;
    
    
%     R111 = (R11+R22+R33)/3;
%     R11 = R111(1:NFFT(1)/2+1);
%       
      
%     R22=R2(1:round(size(R2,1)/2),1,1);
%     R33=R3(1:round(size(R3,1)/2),1,1);
    %
    r = linspace(0,Lx/2,size(u_fft,1)/2+1)/(Lx/2);
end
if (strcmp('2D',flag))
    NFFT = 2.^nextpow2(size(u_fft));
    R1 = ifft2(Rij_x,NFFT(1),NFFT(2))...
                ./NFFT(1)./NFFT(2)./std2(u)^2; %...
                %.*(2*pi)^4; % scaling due to division by 2*pi
    %
    NFFT = 2.^nextpow2(size(v_fft));
    R2 = ifft2(Rij_y,NFFT(1),NFFT(2))...
                ./NFFT(1)./NFFT(2)./std2(v)^2; %...
                %.*(2*pi)^4; % scaling due to division by 2*pi                
    R11 = (R1(:,1)+R2(1,:)')/2;
    R11 = R11(1:size(u_fft)/2+1);
    R22 = (R1(1,:)+R2(:,1)')/2;
    R22 = R22(1:size(u_fft)/2+1);
    %
    r = linspace(0,Lx/2,size(u_fft,1)/2+1)/(Lx/2); % get the radius
    %%
    % 
    % <latex>
    %   From theory we know that the transverse correlation could also be
    %   computed from the longitudinal correlation by
    %   \begin{equation}
    %       g(r) = f + \frac{r}{2}\frac{\partial f}{\partial r}
    %   \end{equation}
    % </latex>
    % 
    g_r = R11 + r'/2.*gradient(R11,max(diff(r)));
end
plot(r,R11,r,R22,r,g_r)
legend('R11','R22','g_r');
h=line([0 1],[0 0],'Color',[0 0 0],'LineWidth',1.0);
% 2D graphs of correlation function
pcolor(fftshift(R1));shading interp;title('R11');
figure
pcolor(fftshift(R2));shading interp;title('R22');
%% Compute length scales
% Computing the length scales is rather easy. The longitudinal and
% transverse length scale are defined through
%%
%
% <latex>
%   \begin{eqnarray}
%       L_{11} &= \int\limits_0^{\infty}R_{11}\,\mathrm{d}r\\
%       L_{22} &= \int\limits_0^{\infty}R_{22}\,\mathrm{d}r
%   \end{eqnarray}
%   Since our data is not represented in an analytical manner we may use a
%   numerical integration routine. Matlab supporty only one numerical
%   integration scheme, namely the Trapezoidal numerical integration. For
%   more information about integration routines you can visit the 
%   \href{http://www.mathworks.de/support/solutions/en/data/1-1679J/index.html}
%   {Mathworks Matlab}
%   support page. 
% </latex>
L11=trapz(r,R11);
L22=trapz(r,R22);
hold on
rectangle('Position',[0,0,L11,1],'LineWidth',2,'LineStyle','--')
%% Spectrum computation
% In general the spectrum of a phyiscal quantity has three dimensions
% whe./1024^4reas the direction in wavenumber space is indicated by $\kappa_1$,
% $\kappa_2$ and $\kappa_3$. Opposed to this relatively extensive
% computation one also might get an idea of the spectral distribution
% calculating the one dimensional spectra. This is achieved by Fourier
% transforming the previously computed correlation functions.
%%
%
% <latex>
%   \begin{equation}
%       E_{ij}(\kappa_1) = \frac{1}{\pi} \int\limits_{-\infty}^{\infty}
%                                \mathbf{R}_{ij}(e_1r_1)\,\mathrm{e}^{-i\kappa_1 r_1}
%                                 \mathrm{d}r_1
%   \end{equation}
% </latex>
%
%%% Compute 1D spectrum
L=length(R11);
NFFT=2^nextpow2(L);
E11=fft(R11,NFFT)/L.*2/pi.*std2(u)^2;
%
L=length(R22);
NFFT=2^nextpow2(L);
E22=fft(R22,NFFT)/L.*2/pi;

f = linspace(0,1,NFFT)*2*pi/dx;

slope=1.5*664092^(2/3)*(f.^(-5/3));
% loglog(f,2*abs(spec_(1:NFFT/2+1)));
% hold on
% loglog(f,slope);
%% Compute 3D spectrum
% In order to avoid aliasing effects usually connected with a one
% dimensional spectrum it is also possible to produce correlations that
% involve all possible directions. The three dimensional Fourier
% transformation of such a correlation produces a spectrum that not only
% depends on a single wavenumber but on the wavenumber vector $\kappa_i$.
% Though the directional information contained in $\kappa_i$ eliminates the
% aliasing problem the complexity makes a physical reasoning impossible.
% For homogeneous isotropic turbulence the situation can be simplified by
% integrating the three dimensional spectrum over spherical shells.
%%
%
% <latex>
%   \begin{equation}
%       E(\kappa) = \oiint E(\boldsymbol\kappa)\mathrm{d}S(\kappa)
%                 = \oiint \frac{1}{2}\,Phi_{ii}(\boldsymbol\kappa)\mathrm{d}S(\kappa)
%   \end{equation}
%   Since the surface of a sphereis completly determined by its radius the
%   surface integral can be solved analytically.
%   \begin{equation}
%       \oiint(\,)\mathrm{d}S(\kappa) = 4\pi\kappa^2\cdot(\,)
%   \end{equation}
% This leads to
%   \begin{equation}
%       E(|\kappa|) = \frac{1}{2}\,\Phi_{ii}(|\boldsymbol\kappa|)
%   \end{equation}
% </latex>
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
                  phi = 0.5*(Rij_x+Rij_y+Rij_z);
%                   phi = fftshift(phi);
%                   phi = phi(1:round(size(Rij_x,1)/2),...
%                             1:round(size(Rij_y,1)/2),...
%                             1:round(size(Rij_z,1)/2));
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
    phi = 0.5*(Rij_x./(size(Rij_x,1)*size(Rij_x,2))...
              +Rij_y./(size(Rij_y,1)*size(Rij_y,2)));
%     phi = phi(round(size(Rij_x,1)/2+1:end),...
%               round(size(Rij_y,1)/2+1:end));
%     phi = phi(1:round(size(phi,1)));
end
%% Compute $\kappa$ vector
% From the previous section we know that the only independent we have in the
% ``system'' is the magnitude of the wave number vector, i.e. 
% $\kappa=|\boldsymbol\kappa|=\sqrt{\kappa_1+\kappa_2+\kappa_3}$. Secondly
% we have to compute the sum $\Phi_{ii}=\Phi_{11}+\Phi_{22}+\Phi_{33}$ and
% take into account its dependence on $|\boldsymbol\kappa|$.
if (strcmp('3D',flag))
    dim = size(phi,1)/2+1;
%     maxdim = sqrt(3*dim^2*(2*pi/Lx)^2);
    E=zeros(round(sqrt(3*dim^2)),1);
    kappa=zeros(round(sqrt(3*dim^2)),1);
    bin_counter=zeros(round(sqrt(3*dim^2)),1);
%     E=zeros(uint64(maxdim),1);
%     kk=zeros(uint64(maxdim),1);
%     bin_counter=zeros(uint64(maxdim),1);
    for k=1:2*(dim-1)
        for j=1:2*(dim-1)
            for i=1:2*(dim-1)
                kx = i*pi/dx;
                ii = i;
                if (i > dim); 
                    kx=(2*(dim)-i)*pi/dx;
                    ii=(2*(dim)-i);
                end
                ky = j*pi/dy;
                jj = j;
                if (j > dim); 
                    ky=(2*(dim)-j)*pi/dy;
                    jj=(2*(dim)-j);
                end
                kz = k*pi/dz;
                kk = k;
                if (k > dim); 
                    kz=(2*(dim)-k)*pi/dz;
                    kk=(2*(dim)-k);
                end
                kappa_pos = round(sqrt(ii^2+jj^2+kk^2)); 
                kappa(kappa_pos) = sqrt(kx^2+ky^2+kz^2);
                E(kappa_pos) = E(kappa_pos) + phi(i,j,k);
                bin_counter(kappa_pos) = bin_counter(kappa_pos) + 1;
            end
        end
    end
    E1=E*4*pi./bin_counter.*kappa.^2;
%     E2=E;
% E=E.*kk.^2;
end
if (strcmp('2D',flag))
    dim = size(phi,1);
    maxdim = sqrt(dim^2*(2*pi/Lx)^2+dim^2*(2*pi/Ly)^2);
    E=zeros(uint64(sqrt(2*dim^2)),1);
    kk=zeros(uint64(sqrt(2*dim^2)),1);
    bin_counter=zeros(uint64(sqrt(2*dim^2)),1);
    for j=1:dim
        for i=1:dim
            kappa=sqrt(i*i*(2*pi/Lx)^2+j*j*(2*pi/Ly)^2);
            kappa_pos=uint64(sqrt(i*i+j*j));
            E(kappa_pos) = E(kappa_pos) + phi(i,j);
			bin_counter(kappa_pos) = bin_counter(kappa_pos) + 1;
            kk(kappa_pos) = kappa;
        end
    end
    E=E*2*pi.*kk./bin_counter;
%     EEE = E*2*pi.*kk;
end
%% Compute 1D spectrum
%
slope=1.5*664092^(2/3)*(kk.^(-5/3));
% test=importdata('INPUT/2D/CTRL_TURB_ENERGY');
%
% dissip=664092;
dissip=664092;
up=17;
L=Lx;
kkke=kk./(2*pi)*L;
kkkd=kk./(2*pi*100)*L;
VKP = 1.5*17^5/dissip.*(kkke).^4./(1+kkke.^2).^(17/6).*exp(-3/2*1.5.*(kkkd).^(4/3));
%
loglog(kk,slope,kk,VKP,kk(2:end),E(2:end))
ylim([1e-14 10]);
h=legend('Kolmogorov','VKP','Computed');
set(h,'Location','SouthWest')
