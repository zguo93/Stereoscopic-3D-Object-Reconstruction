function [ fuse ] = stereorecon( )
tic;
%define the path here
addpath('D:\Final_Year\FYP-infrared\FYP\matlab\imagematch');

 temp = imread('captures\rec_L1.jpg');
obj1 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\rec_L2.jpg');
obj2 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\rec_L3.jpg');
obj3 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\ref_L1.jpg');
ref1 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\ref_L2.jpg');
ref2 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\ref_L3.jpg');
ref3 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
yUp=600;
yDown=1400;
xRight=1100;
xLeft=1900;
obj1=obj1([yUp:yDown],[xRight:xLeft]);
obj2=obj2([yUp:yDown],[xRight:xLeft]);
obj3=obj3([yUp:yDown],[xRight:xLeft]);
ref3=ref3([yUp:yDown],[xRight:xLeft]);
ref2=ref2([yUp:yDown],[xRight:xLeft]);
ref1=ref1([yUp:yDown],[xRight:xLeft]);

temp = imread('captures\bg_L1.jpg');
BimageL = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
bgL(:,:,1)=temp([yUp:yDown],[xRight:xLeft],1);
bgL(:,:,2)=temp([yUp:yDown],[xRight:xLeft],2);
bgL(:,:,3)=temp([yUp:yDown],[xRight:xLeft],3);
BimageL=BimageL([yUp:yDown],[xRight:xLeft]);
imgL=double(BimageL);
imgL(imgL<30)=0;
imgL(imgL~=0)=1;
[ Q datamodulation backgroundmask phase angle angle1] = phaseunwrap( double(obj1),double(obj2),double(obj3),double(ref1),double(ref2),double(ref3));
angle1=angle1-7*2*pi;
[unwphx_qual_angle] = run_unwrapGold3(angle, imgL);
[unwphx_qual_angle1] = run_unwrapGold3(angle1, imgL);
Dimension=unwphx_qual_angle-unwphx_qual_angle1;
L=Dimension;
% figure,imagesc(Dimension)


ref = imread('captures\bg_L1.jpg');
ref_trunc = ref(yUp:yDown,xRight:xLeft );



temp = imread('captures\rec_R1.jpg');
obj1 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\rec_R2.jpg');
obj2 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\rec_R3.jpg');
obj3 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\ref_R1.jpg');
ref1 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\ref_R2.jpg');
ref2 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\ref_R3.jpg');
ref3 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
yUp=700;
yDown=1500;
xRight=400;
xLeft=1200;
obj1=obj1([yUp:yDown],[xRight:xLeft]);
obj2=obj2([yUp:yDown],[xRight:xLeft]);
obj3=obj3([yUp:yDown],[xRight:xLeft]);
ref3=ref3([yUp:yDown],[xRight:xLeft]);
ref2=ref2([yUp:yDown],[xRight:xLeft]);
ref1=ref1([yUp:yDown],[xRight:xLeft]);

temp = imread('captures\bg_R1.jpg');
BimageR = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
bgR(:,:,1)=temp([yUp:yDown],[xRight:xLeft],1);
bgR(:,:,2)=temp([yUp:yDown],[xRight:xLeft],2);
bgR(:,:,3)=temp([yUp:yDown],[xRight:xLeft],3);
BimageR=BimageR([yUp:yDown],[xRight:xLeft]);
imgR=double(BimageR);
imgR(imgR<80)=0;
imgR(imgR~=0)=1;
[ Q datamodulation backgroundmask phase angle angle1] = phaseunwrap( double(obj1),double(obj2),double(obj3),double(ref1),double(ref2),double(ref3));
angle=angle-6*2*pi;
[unwphx_qual_angle] = run_unwrapGold3(angle, imgR);
[unwphx_qual_angle1] = run_unwrapGold3(angle1, imgR);
Dimension=unwphx_qual_angle1-unwphx_qual_angle;
R=Dimension;
% figure,imagesc(Dimension)



% ref = imread('captures\bg_L1.jpg');
% ref_trunc = ref(600:1400,1100:1900 );

reference = imagematch;
reference.set_image(ref_trunc);
reference.get_featurepoints();
reference.get_descriptors();
reference.plot_featurepoints();
cur = imread('captures\bg_R1.jpg');
cur_trunc = cur(yUp:yDown,xRight:xLeft );

current = imagematch;
current.set_image(cur_trunc);
current.get_featurepoints();
current.get_descriptors();
current.plot_featurepoints();
matchpoints = imagematch.get_matchpoints(reference,current);
imagematch.plot_matchpoints(reference,current,matchpoints);
[bestmatrix BestDistance bestInIdx bestInNum avX avY angle amp] = ransac_demo(matchpoints,5,2500,500,0.30,reference,current);
[fuse output1 output2] = fuseiamge( round(avX),round(avY),L,R,0,BimageL,BimageR,imgL,imgR);

toc;


end

