function [spectrum,k,time] = PowerSpec(u,v,w,L,dim)
	tic;
	NFFT = 2.^nextpow2(size(u)); % next power of 2 fitting the length of u

	uu_fft=fftn(u);
	vv_fft=fftn(v);
	ww_fft=fftn(w);
	
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
	rx=[0:1:dim-1] - (dim-1)/2;
    ry=[0:1:dim-1] - (dim-1)/2;
    rz=[0:1:dim-1] - (dim-1)/2;
    
    
	test_x=circshift(rx',[(dim+1)/2 1]);
	test_y=circshift(ry',[(dim+1)/2 1]);
	test_z=circshift(rz',[(dim+1)/2 1]);
	
	[X,Y,Z]= meshgrid(test_x,test_y,test_z);
	r=(sqrt(X.^2+Y.^2+Z.^2));

    dx=2*pi/L;
    k=[1:(dim-1)/2].*dx;
    spectrum=zeros(size(k,2),1);
	for N=2:(dim-1)/2-1

        picker = (r(:,:,:)*dx <= (k(N+1) + k(N))/2) & ...
                 (r(:,:,:)*dx > (k(N) + k(N-1))/2);
		spectrum(N) = sum(muu(picker))+sum(mvv(picker))+sum(mww(picker));

    end
    % special handling for first and last energy value necessary
    picker = (r(:,:,:)*dx <= (k(2) + k(1))/2);
    spectrum(1) = sum(muu(picker))+sum(mvv(picker))+sum(mww(picker));
    picker = (r(:,:,:)*dx > (k(end) + k(end-1))/2);
    spectrum(end) = sum(muu(picker))+sum(mvv(picker))+sum(mww(picker));
	spectrum=0.5*spectrum./(2*pi/L);%(2*pi)^3;%

	time=toc;
end
