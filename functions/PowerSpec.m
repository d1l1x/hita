function [spectrum,k,mu,mv,mw,time] = power_spec(u,v,w,L)
% nx = size(u,1);
% ny = size(u,2);
% nz = size(u,3);
tic;
NFFT = 2.^nextpow2(size(u)); % next power of 2 fitting the length of u
u_fft=fftn(u,NFFT);
v_fft=fftn(v,NFFT);
w_fft=fftn(w,NFFT);

% Calculate the numberof unique points
NumUniquePts = ceil((NFFT(1)+1)/2);

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
% The Nyquist component, if it exists, is unique and should not be multiplied by 2.

if rem(NFFT, 2) % odd nfft excludes Nyquist point
  mu(2:end,2:end,2:end) = mu(2:end,2:end,2:end)*2;
  mv(2:end,2:end,2:end) = mv(2:end,2:end,2:end)*2;
  mw(2:end,2:end,2:end) = mw(2:end,2:end,2:end)*2;
else
  mu(2:end -1,2:end -1,2:end -1) = mu(2:end -1,2:end -1,2:end -1)*2;
  mv(2:end -1,2:end -1,2:end -1) = mv(2:end -1,2:end -1,2:end -1)*2;
  mw(2:end -1,2:end -1,2:end -1) = mw(2:end -1,2:end -1,2:end -1)*2;
end
% Compute the radius vector along which the energies are sumed
mx=NumUniquePts;
my=NumUniquePts;
mz=NumUniquePts;

dx=pi/L;
dy=pi/L;
dz=pi/L;
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

% P=mod(nx,2);
% if (P < 1)
%     Nmax=mx-0.5;
% else
    Nmax=mx;
% end
spectrum=zeros(Nmax,1);
for N=1:Nmax
    Radius1=sqrt(3)*(N-1)*dx; %lower radius bound
    Radius2=sqrt(3)*N*dx; %upper radius bound
    % bild logical index for selecting values lying on the shell
    logical = (Radius1 <= R(:,:,:)) & (R(:,:,:) < Radius2);
    % build summation over shell components
    T_EVP1=sum(mu(logical))+sum(mv(logical))+sum(mw(logical));
    % put them at position N in the spectrum
    spectrum(N)=T_EVP1./2;               
end
k=[1:Nmax].*dx;
spectrum = 1./(2*pi)^3.*spectrum;
time=toc;
end