function [y] = genFringes2(T, m, n,fringeN,intensity)

%% one shot no markers
% y = [zeros(1, Fs/2) ones(1, Fs) zeros(1, Fs/2)]*255;
% y = repmat(y,[m n/size(y,2)]);
% figure, imagesc(y), colormap(gray);
% y = uint8(y);
% y = repmat(y,[1 1 3]);

%% oneshot 

% tt = tt/n;
% fh = 28.44444444444444445*3; 
% psi1  = 2*pi*fh*tt + pi;
% A = 128;
% B1 = 128;1024






tt = 1:n;
y = (1+cos(2*pi*tt/T + fringeN*2*pi/3))*intensity;
y = repmat(y,[m 1]);
y = uint8(y);



% 
% tt = 1:m;
% y = rem(tt,255);
% y = repmat(y',[1 n]);
% y = uint8(y);



% b = Fs/2;
% 
% y = mod(tt-1, Fs);
% ndx1 = find(y < b);
% ndx2 = find(y >= b);
% y(ndx1) = 1;
% y(ndx2) = -1; 
% A = 127;
% B1 = 128;
% y = A + B1 * y;
% y2 = circshift(y',x)';
% y = uint8(round(repmat(y,   [m 1]))); 
% y2 = uint8(round(repmat(y2,   [m 1]))); 

% figure, imagesc(y), colormap(gray);

