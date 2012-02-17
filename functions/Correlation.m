function [R11,R22,r,R1,R2,R3]=Correlation(u,v,w,Lx,dim)
    scaling = 1;
    NFFT = 2.^nextpow2(size(u)); %power of 2 fitting length of u
    u_fft=fftn(u,NFFT)./scaling; 
    NFFT = 2.^nextpow2(size(v));
    v_fft=fftn(v,NFFT)./scaling;
    NFFT = 2.^nextpow2(size(w));
    w_fft=fftn(w,NFFT)./scaling;

    Rij_x=(u_fft.*conj(u_fft)); % compute velo. correlation tensor
    Rij_y=(v_fft.*conj(v_fft));
    Rij_z=(w_fft.*conj(w_fft));
    
    % x-component
    NFFT = 2.^nextpow2(size(u_fft));
    R1=ifftn(Rij_x,NFFT)/std2(u)^2/dim^3;

    % y-component
    NFFT = 2.^nextpow2(size(v_fft));
    R2=ifftn(Rij_y,NFFT)/std2(v)^2./dim^3;
    % z-component
    NFFT = 2.^nextpow2(size(w_fft));
    R3=ifftn(Rij_z,NFFT)/std2(w)^2./dim^3;
    
    R11 = (reshape(R3(1,1,:),NFFT(1),1)+R2(1,:,1)'+R1(:,1,1))/3;
    R11 = R11(1:size(u_fft)/2+1);
    %
    R1_22 = (R1(1,:,1)+R3(1,:,1))/2;
    R2_22 = (R2(:,1,1)+R3(:,1,1))/2;
    R3_22 = (reshape(R1(1,1,:),size(u_fft,1),1)+...
             reshape(R2(1,1,:),size(u_fft,1),1))/2;
    
    R22 = (R1_22'+R2_22+R3_22)/3;
    R22 = R22(1:size(u_fft)/2+1);

    r = linspace(0,Lx/2,size(u_fft,1)/2+1)/(Lx/2);
end