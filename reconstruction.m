clear
temp = imread('captures\rec_f1.jpg');
obj1 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\rec_f2.jpg');
obj2 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\rec_f3.jpg');
obj3 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\ref_f1.jpg');
ref1 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\ref_f2.jpg');
ref2 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
temp = imread('captures\ref_f3.jpg');
ref3 = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
yUp=1100;
yDown=1600;
xRight=1600;
xLeft=2100;
obj1=obj1([yUp:yDown],[xRight:xLeft]);
obj2=obj2([yUp:yDown],[xRight:xLeft]);
obj3=obj3([yUp:yDown],[xRight:xLeft]);
ref3=ref3([yUp:yDown],[xRight:xLeft]);
ref2=ref2([yUp:yDown],[xRight:xLeft]);
ref1=ref1([yUp:yDown],[xRight:xLeft]);

temp = imread('captures\bg_f1.jpg');
image = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
yUp=1100;
yDown=1600;
xRight=1600;
xLeft=2100;
image=image([yUp:yDown],[xRight:xLeft]);
bg(:,:,1)=temp([yUp:yDown],[xRight:xLeft],1);
bg(:,:,2)=temp([yUp:yDown],[xRight:xLeft],2);
bg(:,:,3)=temp([yUp:yDown],[xRight:xLeft],3);

img=double(image);
img(img<100)=0;
img(img~=0)=1;



temp=imread('captures\AMB_rec_f1.jpg');
AMB_rec = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
yUp=1100;
yDown=1600;
xRight=1600;
xLeft=2100;
AMB_rec=AMB_rec([yUp:yDown],[xRight:xLeft]);
figure, imagesc(uint8(double(AMB_rec).*img));
[Markx,Marky]=ginput(1);
round(Markx)

round(Marky)
prompt = 'count the fringes number of object ';
k1 = input(prompt);

temp=imread('captures\AMB_ref_f1.jpg');
AMB_ref= .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);
yUp=1100;
yDown=1600;
xRight=1500;
xLeft=2000;
AMB_ref=AMB_ref([yUp:yDown],[xRight:xLeft]);
figure, imagesc(AMB_ref);
[Markx,Marky]=ginput(1);
AMB_ref(round(Marky)-2:round(Marky)+2,round(Markx)-2:round(Markx)+2)=255;
figure, imagesc(uint8(double(AMB_ref).*img));

round(Markx)

round(Marky)
prompt = 'count the fringes number of reference ';
k2 = input(prompt);

% figure,imagesc(ref1);
% [Markx,Marky]=ginput(1);
% Marky=round(Marky);
% Markx=round(Markx);
% ref1(Marky,Markx)=sum(sum(ref1(Marky-3:Marky+3,Markx-3:Markx+3)))/49;
% ref2(Marky,Markx)=sum(sum(ref1(Marky-3:Marky+3,Markx-3:Markx+3)))/49;
% ref3(Marky,Markx)=sum(sum(ref1(Marky-3:Marky+3,Markx-3:Markx+3)))/49;
% 
% 
% 
% figure,imagesc(obj1);
% [Markx_obj,Marky_obj]=ginput(1);
% Marky_obj=round(Marky_obj);
% Markx_obj=round(Markx_obj);
% obj1(Marky,Markx)=sum(sum(obj1(Marky-3:Marky+3,Markx-3:Markx+3)))/49;
% obj2(Marky,Markx)=sum(sum(obj2(Marky-3:Marky+3,Markx-3:Markx+3)))/49;
% obj3(Marky,Markx)=sum(sum(obj3(Marky-3:Marky+3,Markx-3:Markx+3)))/49;


[ Q datamodulation backgroundmask phase angle angle1] = phaseunwrap( double(obj1),double(obj2),double(obj3),double(ref1),double(ref2),double(ref3));

% angle=angle-angle(Marky_obj,Markx_obj);
% angle1=angle1-angle1(Marky,Markx);
if k1>=k2
    angle1=angle1-(k1-k2)*2*pi;
else
    angle=angle-(k1-k2)*2*pi;
end
[unwphx_qual_angle] = run_unwrapGold3(angle, img);
[unwphx_qual_angle1] = run_unwrapGold3(angle1, img);
Dimension=unwphx_qual_angle-unwphx_qual_angle1;
figure,imagesc(Dimension)
figure,warp(Dimension,bg)
