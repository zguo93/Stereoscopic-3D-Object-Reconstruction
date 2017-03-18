function [ Q datamodulation backgroundmask phase angle angle1] = phaseunwrap( ref1,ref2,ref3,org1,org2,org3 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[r,c]=size(org1);
angle=zeros(r,c);%reference ringe pattern angle
angle1=zeros(r,c);%real object pattern
partialY=zeros(r,c);
partialX=zeros(r,c);
Q=zeros(r,c);%quality of the pixels
mask=zeros(r,c,3);
backgroundmask=zeros(r,c);
intensityofref=zeros(r,c);%I of the reference
intensityoforg=zeros(r,c);%I of the orginal object pattern

sobelx=[-8 0 8; -8 0 8; -8 0 8];
sobely=[8 8 8; 0 0 0;-8 -8 -8];
edgex=zeros(r,c);
edgey=zeros(r,c);

% ref1=zeros(r,c);
% ref2=zeros(r,c);
% ref3=zeros(r,c);
% org1=zeros(r,c);
% org2=zeros(r,c);
% org3=zeros(r,c);
% filter=[ones(1,0.04*c) zeros(1,c*(1-0.04)+1)];
% mfilter=ifft(filter,c);
% Lfilter=mfilter(1:c);
% for i=1:1:r
%     ref1(i,:)=ifft(fft(O1(i,:)).*filter);
% end
% for i=1:1:r
%     ref2(i,:)=ifft(fft(O2(i,:)).*filter);
% end
% for i=1:1:r
%     ref3(i,:)=ifft(fft(O3(i,:)).*filter);
% end
% for i=1:1:r
%     org1(i,:)=ifft(fft(R1(i,:)).*filter);
% end
% for i=1:1:r
%     org2(i,:)=ifft(fft(R2(i,:)).*filter);
% end
% for i=1:1:r
%     org3(i,:)=ifft(fft(R3(i,:)).*filter);
% end

datamodulation=(3*(double(ref1)-double(ref3)).^2+(2*double(ref2)-double(ref1)-double(ref3)).^2).^0.5./(double(ref1)+double(ref2)+double(ref3));


for i=1:1:r%the similar process but for org1 to org3
    for k=1:1:c%calculate the wrapped angle of original pattern
        if ((org2(i,k)>=org1(i,k))&&(org1(i,k)>=org3(i,k))&&(org2(i,k)>=org3(i,k)))
            N=1;
            intensityoforg(i,k)=(org1(i,k)-org3(i,k))/(org2(i,k)-org3(i,k));
            
            angle1(i,k)=pi/3*(2*round((N-1)/2)+intensityoforg(i,k)*(-1)^(N-1));
            
        else if ((org1(i,k)>=org2(i,k))&&(org1(i,k)>=org3(i,k))&&(org2(i,k)>=org3(i,k)))
                N=2;
                intensityoforg(i,k)=(org2(i,k)-org3(i,k))/(org1(i,k)-org3(i,k));
                
                angle1(i,k)=pi/3*(2*round((N-1)/2)+intensityoforg(i,k)*(-1)^(N-1));
                
            else if ((org1(i,k)>=org2(i,k))&&(org1(i,k)>=org3(i,k))&&(org3(i,k)>=org2(i,k)))
                    N=3;
                    intensityoforg(i,k)=(org3(i,k)-org2(i,k))/(org1(i,k)-org2(i,k));
                    
                    angle1(i,k)=pi/3*(2*round((N-1)/2)+intensityoforg(i,k)*(-1)^(N-1));
                    
                else if ((org3(i,k)>=org2(i,k))&&(org1(i,k)>=org2(i,k))&&(org3(i,k)>=org1(i,k)))%3>1>2
                        N=4;
                        intensityoforg(i,k)=(org1(i,k)-org2(i,k))/(org3(i,k)-org2(i,k));
                        
                        angle1(i,k)=pi/3*(2*round((N-1)/2)+intensityoforg(i,k)*(-1)^(N-1));
                        
                    else if ((org3(i,k)>=org2(i,k))&&(org2(i,k)>=org1(i,k))&&(org3(i,k)>=org1(i,k)))
                            N=5;
                            intensityoforg(i,k)=(org2(i,k)-org1(i,k))/(org3(i,k)-org1(i,k));
                            
                            angle1(i,k)=pi/3*(2*round((N-1)/2)+intensityoforg(i,k)*(-1)^(N-1));
                            
                        else if ((org2(i,k)>=org3(i,k))&&(org2(i,k)>=org1(i,k))&&(org3(i,k)>=org1(i,k)))
                                N=6;
                                intensityoforg(i,k)=(org3(i,k)-org1(i,k))/(org2(i,k)-org1(i,k));
                                
                                angle1(i,k)=pi/3*(2*round((N-1)/2)+intensityoforg(i,k)*(-1)^(N-1));
                                
                            end
                        end
                    end
                end
            end
        end
    end
end
for i=1:1:r% the process for ref1 to ref3
    for k=1:1:c%discuss the 6 situations to calculate the phase wrapped.
        if ((ref2(i,k)>=ref1(i,k))&&(ref1(i,k)>=ref3(i,k))&&(ref2(i,k)>=ref3(i,k)))
            N=1;
            intensityofref(i,k)=(ref1(i,k)-ref3(i,k))/(ref2(i,k)-ref3(i,k));
            
            angle(i,k)=pi/3*(2*round((N-1)/2)+intensityofref(i,k)*(-1)^(N-1));
            
        else if ((ref1(i,k)>=ref2(i,k))&&(ref1(i,k)>=ref3(i,k))&&(ref2(i,k)>=ref3(i,k)))
                N=2;
                intensityofref(i,k)=(ref2(i,k)-ref3(i,k))/(ref1(i,k)-ref3(i,k));
                
                angle(i,k)=pi/3*(2*round((N-1)/2)+intensityofref(i,k)*(-1)^(N-1));
                
            else if ((ref1(i,k)>=ref2(i,k))&&(ref1(i,k)>=ref3(i,k))&&(ref3(i,k)>=ref2(i,k)))
                    N=3;
                    intensityofref(i,k)=(ref3(i,k)-ref2(i,k))/(ref1(i,k)-ref2(i,k));
                    
                    angle(i,k)=pi/3*(2*round((N-1)/2)+intensityofref(i,k)*(-1)^(N-1));
                    
                else if ((ref3(i,k)>=ref2(i,k))&&(ref1(i,k)>=ref2(i,k))&&(ref3(i,k)>=ref1(i,k)))
                        N=4;
                        intensityofref(i,k)=(ref1(i,k)-ref2(i,k))/(ref3(i,k)-ref2(i,k));
                        
                        angle(i,k)=pi/3*(2*round((N-1)/2)+intensityofref(i,k)*(-1)^(N-1));
                        
                    else if ((ref3(i,k)>=ref2(i,k))&&(ref2(i,k)>=ref1(i,k))&&(ref3(i,k)>=ref1(i,k)))
                            N=5;
                            intensityofref(i,k)=(ref2(i,k)-ref1(i,k))/(ref3(i,k)-ref1(i,k));
                            
                            angle(i,k)=pi/3*(2*round((N-1)/2)+intensityofref(i,k)*(-1)^(N-1));
                            
                        else if ((ref2(i,k)>=ref3(i,k))&&(ref2(i,k)>=ref1(i,k))&&(ref3(i,k)>=ref1(i,k)))
                                N=6;
                                intensityofref(i,k)=(ref3(i,k)-ref1(i,k))/(ref2(i,k)-ref1(i,k));
                                
                                angle(i,k)=pi/3*(2*round((N-1)/2)+intensityofref(i,k)*(-1)^(N-1));
                                
                            end
                        end
                    end
                end
            end
        end    
    end
end
% for i=1:1:r
%     for k=1:1:c
%         angle(i,k)=atan(3^0.5*(ref1(i,k)-ref3(i,k))./(2*ref2(i,k)-ref1(i,k)-ref3(i,k)))+pi/2;
%     end
% end
% for i=1:1:r
%     for k=1:1:c
%         angle1(i,k)=atan(3^0.5*(org1(i,k)-org3(i,k))./(2*org2(i,k)-org1(i,k)-org3(i,k)))+pi/2;
%     end
% end

normalizedP=angle./2*pi;
for i=2:1:r-1
    for k=1:1:c
        partialY(i,k)=max(abs(normalizedP(i,k)-normalizedP(i-1,k)-round(normalizedP(i,k)-normalizedP(i-1,k))),abs(normalizedP(i+1,k)-normalizedP(i,k)-round(normalizedP(i+1,k)-normalizedP(i,k))));
    end
end
partialY(1,:)=normalizedP(2,:)-normalizedP(1,:)-round(normalizedP(2,:)-normalizedP(1,:));
partialY(r,:)=normalizedP(r,:)-normalizedP(r-1,:)-round(normalizedP(r,:)-normalizedP(r-1,:));
for i=1:1:r
    for k=2:1:c-1
        partialX(i,k)=max(abs(normalizedP(i,k)-normalizedP(i,k-1)-round(normalizedP(i,k)-normalizedP(i,k-1))),abs(normalizedP(i,k+1)-normalizedP(i,k)-round(normalizedP(i,k+1)-normalizedP(i,k))));
    end
end
partialX(:,1)=normalizedP(:,2)-normalizedP(:,1)-round(normalizedP(:,2)-normalizedP(:,1));
partialX(:,c)=normalizedP(:,c)-normalizedP(:,c-1)-round(normalizedP(:,c)-normalizedP(:,c-1));
for i=1:1:r
    for k=1:1:c
        Q(i,k)=max(partialX(i,k),partialY(i,k));
    end
end
Qmean=sum(sum(Q))/r/c;
QstandardD=sum(sum(((Q-Qmean).*(Q-Qmean)).^0.5))/r/c;
threshold=[Qmean Qmean+2*QstandardD Qmean+4*QstandardD];
for i=1:1:r
    for k=1:1:c
        if (datamodulation(i,k)<0.705&&datamodulation(i,k)>=0.17&&ref1(i,k)>5)
            backgroundmask(i,k)=1;
        end
    end
end
for i=1:1:r
    for k=1:1:c
        if (Q(i,k)<threshold(1)&&backgroundmask(i,k)==1)
            mask(i,k,1)=1;
        end
    end
end
for i=1:1:r
    for k=1:1:c
        if (Q(i,k)>=threshold(1)&&Q(i,k)<threshold(2)&&backgroundmask(i,k)==1)
            mask(i,k,2)=1;
        end
    end
end
for i=1:1:r
    for k=1:1:c
        if (Q(i,k)>=threshold(1)&&Q(i,k)>=threshold(2)&&backgroundmask(i,k)==1)
            mask(i,k,3)=1;
        end
    end
end

xCenter=round(c/2);
yCenter=round(r/2);
for i=yCenter-4:1:yCenter+4
    for k=xCenter-4:1:xCenter+4
        if(mask(i,k,1)==1)
            Centerx=k;
            Centery=i;
        end
    end
end

for i=2:1:r-1
    for k=2:1:c-1
        edgex(i,k)=sum(sum(angle(i-1:i+1,k-1:k+1).*sobelx));%calculate the edge by sobel operator in horizontal direction
    end
end
for i=2:1:r-1
    for k=2:1:c-1
        edgey(i,k)=sum(sum(angle(i-1:i+1,k-1:k+1).*sobely));%calculate the edge by sobel operator in vertical direction
    end
end





phase=angle-angle1;









end

