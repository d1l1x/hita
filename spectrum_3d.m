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
path('./functions',path)

%% Read data files
% Read in the data files and measure the time for reading. The output of
% the tic/toc block is in seconds. What you should get from the tic/toc
% block is that most of the time is spend during data I/O. The actual
% computation needs only ??? of the time of the I/O operations.
%
[uvel,vvel,wvel,time_read] = read_data(datadir,flag);
%% Set some neccessary parameters
% For further computations it is important to define some parmeters of the
% DNS simulation such as domain size, grid spacing, and the number of grid
% points.
%
    dim=256; % number of points in one dimension
    Lx=5e-3; % domain size
    Ly=Lx;
    Lz=Lx;
    dx=Lx/dim; % grid spacing
    dy=dx;
    dz=dx;
    nu=1.7e-5; % viscosity
    u=reshape(uvel,dim,dim,dim); % reshape arrays to have them in 3D
    v=reshape(vvel,dim,dim,dim);
    w=reshape(wvel,dim,dim,dim);
    clear uvel vvel wvel
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
[spectrum,k] = power_spec(u,v,w,Lx);
%% Compute dissipation and turbulent kinetic energy
[Dissipation,kinetic_E] = spec_prop(spectrum,k,nu);

%% Verify computation of kinetic Energy

up = sqrt(1/3*(u.^2+v.^2+w.^2));
kinetic_Up = sum(sum(sum(3/2*up.^2)))/size(up,1)^3;
%% Kolmogrov properties
eta = (nu^3/Dissipation)^(1/4);
u_eta = (nu*Dissipation)^(1/4);
tau = (nu/Dissipation)^(1/2);
            
%% Compute model spectra
%
% Von Karman-Pao Spektren
kd = k(end);
ke = pi/Lx;
A = 1.5;
up = mean2(up);

VKP1 = A*up^5/Dissipation.*(k./ke).^4./(1+(k./ke).^2).^(17/6).*exp(-3/2*A.*(k./kd).^(4/3));

kd = 1./eta;
VKP2 = 1.5*(k./kd).^(-5/3)./(Dissipation*nu^5)^(-1/4).*exp(-1.5*1.5.*(k./kd).^(4/3));

% Kolmogorov Spektrum
Kolmo=1.5*Dissipation^(2/3)*(k.^(-5/3));

% Plot spectra
loglog(k,Kolmo,k,VKP1,k,VKP2,k,spectrum)

h=legend('Kolmogorov','VKP1','VKP2','Computed');
set(h,'Location','SouthWest')
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
[R11,R22,r,R1,R2,R3]=correlation(u,v,w,Lx);