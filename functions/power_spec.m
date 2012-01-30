function [spectrum,k] = power_spec(u,v,w,L)
nx = size(u,1);
ny = size(u,2);
nz = size(u,3);

NFFT = 2.^nextpow2(size(u)); % next power of 2 fitting the length of u
u_fft=fftn(u,NFFT);
v_fft=fftn(v,NFFT);
w_fft=fftn(w,NFFT);

% Calculate the numberof unique points
NumUniquePts = ceil((NFFT+1)/2);

% FFT is symmetric, throw away second half 
u_fft = u_fft(1:NumUniquePts,1:NumUniquePts,1:NumUniquePts);
v_fft = v_fft(1:NumUniquePts,1:NumUniquePts,1:NumUniquePts);
w_fft = w_fft(1:NumUniquePts,1:NumUniquePts,1:NumUniquePts);

mu = abs(u_fft)/length(u)^3;
mv = abs(v_fft)/length(v)^3;
mw = abs(w_fft)/length(w)^3;

% Take the square of the magnitude of fft of x. 
mu = mu.^2;
mv = mv.^2;
mw = mw.^2;

% Since we dropped half the FFT, we multiply mx by 2 to keep the same energy.
% The DC component and Nyquist component, if it exists, are unique and should not be multiplied by 2.

if rem(NFFT, 2) % odd nfft excludes Nyquist point
  mu(2:end) = mu(2:end)*2;
  mv(2:end) = mv(2:end)*2;
  mw(2:end) = mw(2:end)*2;
else
  mu(2:end -1) = mu(2:end -1)*2;
  mv(2:end -1) = mv(2:end -1)*2;
  mw(2:end -1) = mw(2:end -1)*2;
end
% Compute the radius vector along which the energies are sumed
mx=NumUniquePts;
my=NumUniquePts;
mz=NumUniquePts;

dx=2*pi/L;
dy=2*pi/L;
dz=2*pi/L;

for I=1:mx
       X0(I)=(I)*dx; 
end

for J=1:my
       Y0(J)=(J)*dy;                                  
end

for K=1:mz
       Z0(K)=(K)*dz;                                  
end

for I=1:mx
    for J=1:my
        for K=1:mz
            R(I,J,K)=sqrt(X0(I)*X0(I)+Y0(J)*Y0(J)+Z0(K)*Z0(K));
        end
    end
end

P=mod(nx,2);
if (P < 1)
    Nmax=mx-0.5;
else
    Nmax=mx;
end

% % % for N=1:Nmax
% % %     T_UprimVP1=0;
% % %     T_VprimVP1=0;
% % %     T_WprimVP1=0;
% % %     T_EVP1=0;
% % %     
% % %     Radius1=sqrt(3)*(N-1)*dx;
% % %     Radius2=sqrt(3)*N*dx;
% % %     
% % %     for I=1:mx
% % %       for J=1:my
% % %           for K=1:mz
% % %               if (Radius1 <= R(I,J,K)) && (R(I,J,K) < Radius2)
% % %              
% % %                 T_UprimVP1=T_UprimVP1+mu(I,J,K);
% % %                 T_VprimVP1=T_VprimVP1+mv(I,J,K);
% % %                 T_WprimVP1=T_WprimVP1+mw(I,J,K);
% % %                 T_EVP1=T_UprimVP1+T_VprimVP1+T_WprimVP1;                                  
% % %               else         
% % %                 T_UprimVP1=T_UprimVP1;
% % %                 T_VprimVP1=T_VprimVP1;
% % %                 T_WprimVP1=T_WprimVP1;
% % %                 T_EVP1=T_UprimVP1+T_VprimVP1+T_WprimVP1;                                                  
% % %               end    
% % %           end
% % %       end
% % %     end        
% % %     spectrum(N)=0.25.*T_EVP1;                
% % % end
spectrum=zeros(Nmax,1);
for N=1:Nmax

    Radius1=sqrt(3)*(N-1)*dx; %lower radius bound
    Radius2=sqrt(3)*N*dx; %upper radius bound
    % bild logical index for selecting values lying on the shell
    logical = (Radius1 <= R(:,:,:)) & (R(:,:,:) < Radius2);
    % build summation over shelle over all components
    T_EVP1=sum(mu(logical))+sum(mv(logical))+sum(mw(logical));
    % put them at position N in the spectrum
    spectrum(N)=0.25.*T_EVP1;                
end
k=[1:Nmax].*(2*pi/(L));

end