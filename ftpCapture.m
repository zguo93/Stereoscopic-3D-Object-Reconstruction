% close all; 
% clc;
% clear;
%% capture

% [f fMarker ] = genFringes2(1, 768, 1024);


% dummy
prefix = 'captures/';
prefix2 = 'rec';


%% Capture image
num=1;
dispid=2;
acceptforMARK=0;
width=16;%fringe width
intensity_ref=180;
intensity_rec=110;


% while acceptforMARK==0,
% [fMarker ] = genFringes2(width, 768, 1024,-1,intensity_rec);% -120 degree
% fMarker = repmat(fMarker,[1 1 3]);
% figure,imagesc(fMarker);
% [x,y]=ginput(1);
% x=round(x);
% y=round(y);
% 
%     fMarker(y-1:y+1,x-1:x+1,1)=255;
%     fMarker(y-1:y+1,x-1:x+1,2)=255;
%     fMarker(y-1:y+1,x-1:x+1,3)=255;
% fullscreen(fMarker, dispid);
% prompt = 'Is this position acceptable, y=1 or n=0?';
% acceptforMARK = input(prompt);
% end
% pause(2);
% %take rec image with a mark point
% prefix2='AMB_rec';
% fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_f', 1);
% CanonControl(fnameTexture);




%take 4 pictures with obj
prompt = 'continue to take 4 rec images, y=1 or n=0?';
accept = input(prompt);
if accept==1
    prefix2 = 'rec';
    
    pause(3);
    bg = ones(768,1024)  * 150;
    bg = repmat(uint8(bg), [1 1 3]);
    bgleft=bg;
    bgleft(:,:,2)=uint8(zeros(768,1024));
    bgleft(:,:,3)=uint8(zeros(768,1024));
    
    bgright=bg;
    bgright(:,:,1)=uint8(zeros(768,1024));
    bgright(:,:,2)=uint8(zeros(768,1024));
    fullscreen(bgleft, 3);
    fullscreen2(bgright, 2);
    pause(1);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, 'bg', '_R', 1);
    CameraCapture(0,fnameTexture);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, 'bg', '_L', 1);
    CameraCapture(1,fnameTexture);
    
    answer=0;
    while answer==0
    [ fMarker ] = genFringes2(width, 768, 1024,-1,intensity_rec);%-120 degrees
    fMarker = repmat(fMarker,[1 1 3]);
    fMarkerleft=fMarker;
    fMarkerleft(:,:,2)=uint8(zeros(768,1024));
    fMarkerleft(:,:,3)=uint8(zeros(768,1024));
    

    
    imagesc(fMarkerleft)
    [Markx_L,Marky_L]=ginput(1);
    round(Markx_L)
    round(Marky_L)

    fMarkerright=fMarker;
    
    fMarkerright(:,:,1)=uint8(zeros(768,1024));
    fMarkerright(:,:,2)=uint8(zeros(768,1024));

    
    imagesc(fMarkerright)
    [Markx_R,Marky_R]=ginput(1);
    round(Markx_R)
    round(Marky_R)

    
    fullscreen(fMarkerleft, 3);
    fullscreen2(fMarkerright, 2);
    prompt = 'continue to take 4 rec images, y=1 or n=0?';
    answer = input(prompt);
    end
    
    pause(1);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_R', 1);
    CameraCapture(0,fnameTexture);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_L', 1);
    CameraCapture(1,fnameTexture);
    pause(1);

    
    
    [ fMarker ] = genFringes2(width, 768, 1024,0,intensity_rec);%0 degrees
    fMarker = repmat(fMarker,[1 1 3]);
    fMarkerleft=fMarker;
    fMarkerleft(:,:,2)=uint8(zeros(768,1024));
    fMarkerleft(:,:,3)=uint8(zeros(768,1024));
    

    
    fMarkerright=fMarker;
    fMarkerright(:,:,1)=uint8(zeros(768,1024));
    fMarkerright(:,:,2)=uint8(zeros(768,1024));
    

    
    fullscreen(fMarkerleft, 3);
    fullscreen2(fMarkerright, 2);
    pause(1);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_R', 2);
    CameraCapture(0,fnameTexture);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_L', 2);
    CameraCapture(1,fnameTexture);
    pause(1);

    
    
    
    [fMarker ] = genFringes2(width, 768, 1024,1,intensity_rec);% move toward left for 6 unit, correspoding to +120 degree(THIS IS FOR FS=18)
    fMarker = repmat(fMarker,[1 1 3]);
    fMarkerleft=fMarker;
    fMarkerleft(:,:,2)=uint8(zeros(768,1024));
    fMarkerleft(:,:,3)=uint8(zeros(768,1024));
    

    
    fMarkerright=fMarker;
    fMarkerright(:,:,1)=uint8(zeros(768,1024));
    fMarkerright(:,:,2)=uint8(zeros(768,1024));

    
    fullscreen(fMarkerleft, 3);
    fullscreen2(fMarkerright, 2);
    pause(1);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_R', 3);
    CameraCapture(0,fnameTexture);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_L', 3);
    CameraCapture(1,fnameTexture);
    pause(1);
    
end




prompt = 'continue to take 4 ref images, y=1 or n=0?';
accept = input(prompt);
if accept==1
    prefix2 = 'ref';
     pause(2);
    
    [ fMarker ] = genFringes2(width, 768, 1024,-1,intensity_ref);%-120 degrees
    fMarker = repmat(fMarker,[1 1 3]);
    fMarkerleft=fMarker;
    fMarkerleft(:,:,2)=uint8(zeros(768,1024));
    fMarkerleft(:,:,3)=uint8(zeros(768,1024));
    

    
    fMarkerleft(Marky_L-2:Marky_L+2,Markx_L-2:Markx_L+2,1)=255;
    
    fMarkerright=fMarker;
    fMarkerright(:,:,1)=uint8(zeros(768,1024));
    fMarkerright(:,:,2)=uint8(zeros(768,1024));

    
    fMarkerright(Marky_R-2:Marky_R+2,Markx_R-2:Markx_R+2,3)=255;
    fullscreen(fMarkerleft, 3);
    fullscreen2(fMarkerright, 2);
    pause(1);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_R', 1);
    CameraCapture(0,fnameTexture);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_L', 1);
    CameraCapture(1,fnameTexture);
    pause(1);

    
    
    
    
    [ fMarker ] = genFringes2(width, 768, 1024,0,intensity_ref);
    fMarker = repmat(fMarker,[1 1 3]);
    fMarkerleft=fMarker;
    fMarkerleft(:,:,2)=uint8(zeros(768,1024));
    fMarkerleft(:,:,3)=uint8(zeros(768,1024));
    

    
    fMarkerright=fMarker;
    fMarkerright(:,:,1)=uint8(zeros(768,1024));
    fMarkerright(:,:,2)=uint8(zeros(768,1024));
    

    
    fullscreen(fMarkerleft, 3);
    fullscreen2(fMarkerright, 2);
    pause(1);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_R', 2);
    CameraCapture(0,fnameTexture);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_L', 2);
    CameraCapture(1,fnameTexture);
    pause(1);
    
    
    
    

    [fMarker ] = genFringes2(width, 768, 1024,1,intensity_ref);% move toward left for 6 unit, correspoding to +120 degree(THIS IS FOR FS=18)
    fMarker = repmat(fMarker,[1 1 3]);
    fMarkerleft=fMarker;
    fMarkerleft(:,:,2)=uint8(zeros(768,1024));
    fMarkerleft(:,:,3)=uint8(zeros(768,1024));
    
    
    fMarkerright=fMarker;
    fMarkerright(:,:,1)=uint8(zeros(768,1024));
    fMarkerright(:,:,2)=uint8(zeros(768,1024));
    

    
    fullscreen(fMarkerleft, 3);
    fullscreen2(fMarkerright, 2);
    pause(1);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_R', 3);
    CameraCapture(0,fnameTexture);
    fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_L', 3);
    CameraCapture(1,fnameTexture);
    pause(1);
    
    
%     [fMarker ] = genFringes2(width, 768, 1024,-1,intensity_ref);% -120 degree
%     fMarker = repmat(fMarker,[1 1 3]);
%     fMarker(y-1:y+1,x-1:x+1,1)=255;
%     fMarker(y-1:y+1,x-1:x+1,2)=255;
%     fMarker(y-1:y+1,x-1:x+1,3)=255;
%     fullscreen(fMarker, dispid);
%     prefix2='AMB_ref';
%     fnameTexture = sprintf('%s%s%s%d.jpg', prefix, prefix2, '_f', 1);
%     CanonControl(fnameTexture);
end
% temp = imread('captures/ref_f.jpg');
% 
[ final ] = stereorecon( );

% temp2(:, :) = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);

%  temp = imread('captures\rec_L1.jpg');
% obj1 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% temp = imread('captures\rec_L2.jpg');
% obj2 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% temp = imread('captures\rec_L3.jpg');
% obj3 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% temp = imread('captures\ref_L1.jpg');
% ref1 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% temp = imread('captures\ref_L2.jpg');
% ref2 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% temp = imread('captures\ref_L3.jpg');
% ref3 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% yUp=600;
% yDown=1300;
% xRight=200;
% xLeft=800;
% obj1=obj1([yUp:yDown],[xRight:xLeft]);
% obj2=obj2([yUp:yDown],[xRight:xLeft]);
% obj3=obj3([yUp:yDown],[xRight:xLeft]);
% ref3=ref3([yUp:yDown],[xRight:xLeft]);
% ref2=ref2([yUp:yDown],[xRight:xLeft]);
% ref1=ref1([yUp:yDown],[xRight:xLeft]);
% img=double(image);
% temp = imread('captures\bg_L1.jpg');
% image = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% bgL(:,:,1)=temp([yUp:yDown],[xRight:xLeft],1);
% bgL(:,:,2)=temp([yUp:yDown],[xRight:xLeft],2);
% bgL(:,:,3)=temp([yUp:yDown],[xRight:xLeft],3);
% image=image([yUp:yDown],[xRight:xLeft]);
% img=double(image);
% img(img<80)=0;
% img(img~=0)=1;
% [ Q datamodulation backgroundmask phase angle angle1] = phaseunwrap( double(obj1),double(obj2),double(obj3),double(ref1),double(ref2),double(ref3));
% angle1=angle1-6*2*pi+16;
% [unwphx_qual_angle] = run_unwrapGold3(angle, img);
% [unwphx_qual_angle1] = run_unwrapGold3(angle1, img);
% Dimension=unwphx_qual_angle-unwphx_qual_angle1;
% L=Dimension;
% figure,imagesc(Dimension)
% 
% 
% 
% 
% temp = imread('captures\rec_R1.jpg');
% obj1 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% temp = imread('captures\rec_R2.jpg');
% obj2 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% temp = imread('captures\rec_R3.jpg');
% obj3 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% temp = imread('captures\ref_R1.jpg');
% ref1 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% temp = imread('captures\ref_R2.jpg');
% ref2 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% temp = imread('captures\ref_R3.jpg');
% ref3 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% yUp=750;
% yDown=1450;
% xRight=1600;
% xLeft=2200;
% obj1=obj1([yUp:yDown],[xRight:xLeft]);
% obj2=obj2([yUp:yDown],[xRight:xLeft]);
% obj3=obj3([yUp:yDown],[xRight:xLeft]);
% ref3=ref3([yUp:yDown],[xRight:xLeft]);
% ref2=ref2([yUp:yDown],[xRight:xLeft]);
% ref1=ref1([yUp:yDown],[xRight:xLeft]);
% img=double(image);
% temp = imread('captures\bg_R1.jpg');
% image = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
% bgR(:,:,1)=temp([yUp:yDown],[xRight:xLeft],1);
% bgR(:,:,2)=temp([yUp:yDown],[xRight:xLeft],2);
% bgR(:,:,3)=temp([yUp:yDown],[xRight:xLeft],3);
% image=image([yUp:yDown],[xRight:xLeft]);
% img=double(image);
% img(img<30)=0;
% img(img~=0)=1;
% [ Q datamodulation backgroundmask phase angle angle1] = phaseunwrap( double(obj1),double(obj2),double(obj3),double(ref1),double(ref2),double(ref3));
% angle=angle-6*2*pi;
% [unwphx_qual_angle] = run_unwrapGold3(angle, img);
% [unwphx_qual_angle1] = run_unwrapGold3(angle1, img);
% Dimension=unwphx_qual_angle1-unwphx_qual_angle;
% R=Dimension;
% figure,imagesc(Dimension)


% final=zeros(1104,1271);
% final(104:1904,1:630)=L(1:1801,1:630);
% final(1:1801,631:1271)=R(1:1801,361:1001);