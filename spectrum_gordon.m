function [ENVP1,S1,k,L11,RMS]=spectrum_gordon(x_velocity,y_velocity,z_velocity,L,nx,ny,nz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Input: x_velocity, y_velocity: velocity components;
%        L : length of domain
%        nx,ny, number of points in x and y direction
%     
% Output:
%     ENVP1: Two-dimensional energy spectrum;
%     S1: turbulent kinetic energy
%     L11: longitudinal integral length scale
%     RMS: RMS of velocity fluctuations
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long g

L11=0;
RMS=0;

C=1/(nx*nx*ny*ny*nz*nz);

FVP1=ones(nx,ny,nz);
GVP1=ones(nx,ny,nz);

dx=2*pi/L;
dy=2*pi/L;
dz=2*pi/L;

% Mesh size in physical space
X1=[0:(L/(nx-1)):L];
Y1=[0:(L/(ny-1)):L];
Z1=[0:(L/(nz-1)):L];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% velocities fluctuations fields created
FVP1=x_velocity-sum((sum(sum(x_velocity)))/(nx*ny*nz));
GVP1=y_velocity-sum((sum(sum(y_velocity)))/(nx*ny*nz));
ZVP1=z_velocity-sum((sum(sum(z_velocity)))/(nx*ny*nz));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D FFT - Not normalized

FFT2_FVP1=fftn(FVP1);
FFT2_GVP1=fftn(GVP1);
FFT2_ZVP1=fftn(ZVP1);

% 2D FFT centering - Not normalized
% Normalization is done with energy

BVP1=fftshift((FFT2_FVP1.*conj(FFT2_FVP1)));
CVP1=fftshift((FFT2_GVP1.*conj(FFT2_GVP1)));
DVP1=fftshift((FFT2_ZVP1.*conj(FFT2_ZVP1))); 

XE=[1:nx];
YE=[1:ny];
ZE=[1:nz];

%figure 
%contourf(XE,YE,(BVP1+CVP1),25)
%
% if the total number of points is unpair : the center is one point of the
% mesh, else the center is the center of the square located in the center
% of the total mesh 

% in the following, we determine the energy contained in different rings
% (between Radius1 and Radius2) centered in the center of the spectral mesh.
% For that, we need to evaluate the distance (radius of circles) between all the
% points and the center of the mesh

% (I;J) : coordinates of points 
% mx my : coordinates of the center of the grid
% (X0(I); Y0(J)) : distances (with sign) from the center of the mesh

mx=(nx+1)/2;
my=(ny+1)/2;
mz=(nz+1)/2;

for I=1:nx
       X0(I)=(I-mx)*dx; 
end

for J=1:ny
       Y0(J)=(J-my)*dy;                                  
end

for K=1:nz
       Z0(J)=(K-mz)*dz;                                  
end

for I=1:nx
    for J=1:ny
        for K=1:nz
            R(I,J,K)=sqrt(X0(I)*X0(I)+Y0(J)*Y0(J)+Z0(K)*Z0(K));
        end
    end
end

clear T_UprimVP1 T_VprimVP1 T_EVP1 ENVP1

% Test to know the number of circles (Nmax) needed to cover all the grid
P=mod(nx,2);
if (P < 1)
    Nmax=mx-0.5;
else
    Nmax=mx;
end


for N=1:Nmax
    T_UprimVP1=0;
    T_VprimVP1=0;
    T_WprimVP1=0;
    T_EVP1=0;
    
    Radius1=sqrt(3)*(N-1)*dx;
    Radius2=sqrt(3)*N*dx;
    
    for I=1:nx
      for J=1:ny
          for K=1:nz
              if (Radius1 <= R(I,J,K)) && (R(I,J,K) < Radius2)
             
                T_UprimVP1=T_UprimVP1+BVP1(I,J,K);
                T_VprimVP1=T_VprimVP1+CVP1(I,J,K);
                T_WprimVP1=T_WprimVP1+DVP1(I,J,K);
                T_EVP1=T_UprimVP1+T_VprimVP1+T_WprimVP1;                                  
              else         
                T_UprimVP1=T_UprimVP1;
                T_VprimVP1=T_VprimVP1;
                T_WprimVP1=T_WprimVP1;
                T_EVP1=T_UprimVP1+T_VprimVP1+T_WprimVP1;                                                  
              end    
          end
      end
    end        
    ENVP1(N)=0.25.*T_EVP1;                
end

S1=sum(ENVP1);

% Normalisation with C
S=C*[S1];
% Factor 2 to consider total space
s=sqrt(2*S);

% T : time
T=[0.22572e-3];

k=[1:Nmax].*(2*pi/(L));

% % % figure
% % % hold on
% % % plot(log10(k),log10(C.*ENVP1),'b')
% % % set(gca,'fontname', 'new century schoolbook','fontweight','bold');
% % % title('Energy')
% % % xlabel('log10 (k)','fontweight','bold','fontname', 'new century schoolbook','fontsize',11)
% % % ylabel('log10 (E)','fontweight','bold','fontname', 'new century schoolbook','fontsize',11)
% % % hold off

% % % %%%%%%%%%%%%%%%%%% RMS values of velocity fluctuations %%%%%%%%%%%%
% % % 
% % % Mtime1=sum(sum(x_velocity))/(nx*ny);
% % % Uprimetime1=x_velocity-Mtime1;
% % % Uprimetime12=Uprimetime1.^2;
% % % Urms_time1=sum(sum(Uprimetime12))/(nx*ny);
% % % 
% % % Ntime1=sum(sum(y_velocity))/(nx*ny);
% % % 
% % % Vprimetime1=y_velocity-Ntime1;
% % % Vprimetime12=Vprimetime1.^2;
% % % 
% % % Vrms_time1=sum(sum(Vprimetime12))/(nx*ny);
% % % RMS_time1=sqrt(0.5*(Urms_time1+Vrms_time1));
% % % RMS=[RMS_time1]
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%% RMS values of velocity fluctuations - Spectral calculation %%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % RMS_spec_VP1=sqrt(sum((4/2).*C.*ENVP1))
% % % 
% % % %%%%%%%%%%%%%%%%%% L11 - Spectral and Pseudo spectral calculations 
% % % 
% % % DVP1=sum(2.*C.*ENVP1./k);
% % % L11_VPspec1 = (2/(RMS_spec_VP1^2)).*DVP1;
% % % L11VP1 = (2/(RMS_time1^2)).*DVP1;
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % L11_VPspec=[L11_VPspec1]
% % % L11VP=[L11VP1]
% % % L11=L11_VPspec;
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % RMS=[RMS_time1]





