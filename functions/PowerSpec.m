function [spectrum,k,time] = PowerSpec(u,v,w,L,dim)
	tic;
	NFFT = 2.^nextpow2(size(u)); % next power of 2 fitting the length of u
	%u_fft=fftn(u,NFFT);
	%v_fft=fftn(v,NFFT);
	%w_fft=fftn(w,NFFT);
	% NFFT=33;
	uu_fft=fftn(u,NFFT);
	vv_fft=fftn(v,NFFT);
	ww_fft=fftn(w,NFFT);

	% Calculate the numberof unique points
	%NumUniquePts = ceil((NFFT(1)+1)/2);

	muu = abs(uu_fft)/length(u)^3;
	mvv = abs(vv_fft)/length(v)^3;
	mww = abs(ww_fft)/length(w)^3;

	% Take the square of the magnitude of fft of x. 
	%mu = mu.^2;
	%mv = mv.^2;
	%mw = mw.^2;
	muu = muu.^2;
	mvv = mvv.^2;
	mww = mww.^2;


% % % % % 	for i=1:dim-1
% % % % % 		xx(i) = i-(dim+1)/2;
% % % % % 		yy(i) = i-(dim+1)/2;
% % % % % 		zz(i) = i-(dim+1)/2;
% % % % %     end
    % equivalent see above
    dx=2*pi/L;
    k=[1:NFFT(1)/2].*dx;
	rx=[0:1:NFFT(1)-1] - NFFT(1)/2+1;
    ry=[0:1:NFFT(1)-1] - NFFT(1)/2+1;
    rz=[0:1:NFFT(1)-1] - NFFT(1)/2+1;
    
    
	test_x=circshift(rx',[NFFT(1)/2+1 1])*dx;
	test_y=circshift(ry',[NFFT(1)/2+1 1])*dx;
	test_z=circshift(rz',[NFFT(1)/2+1 1])*dx;
	
	[X,Y,Z]= meshgrid(test_x,test_y,test_z);
	r=(sqrt(X.^2+Y.^2+Z.^2));

    spectrum=zeros(NFFT(1)/2,1);
	for N=2:size(k,2)-1
        picker = (r(:,:,:) <= (k(N+1) + k(N))/2) & ...
                 (r(:,:,:) > (k(N) + k(N-1))/2);
		spectrum(N) = sum(muu(picker))+sum(mvv(picker))+sum(mww(picker));
    end
    % special handling for first and last energy value necessary
    picker = (r(:,:,:) <= (k(2) + k(1))/2);
    spectrum(1) = sum(muu(picker))+sum(mvv(picker))+sum(mww(picker));
    picker = (r(:,:,:) > (k(end) + k(end-1))/2);
    spectrum(end) = sum(muu(picker))+sum(mvv(picker))+sum(mww(picker));
	spectrum=0.5*spectrum./(2*pi/L);%(2*pi)^3;%


	time=toc;
end
