function [fuse record_output1 record_output2] = fuseiamge( X,Y,L,R,windowsize,bgL,bgR,imgL,imgR)
%UNTITLED Summary of this function goes here
%   [fuse output1 output2] = mindifference( round(avX),round(avY),L,R );
[r,c]=size(L);        output1=zeros(2*r,2*c);
output2=zeros(2*r,2*c);
bg1=uint8(zeros(2*r,2*c));
bg2=uint8(zeros(2*r,2*c));
bgL=bgL+44;%compensate for color intensity
if Y>=0 && X<=c

output1(1+Y:r+Y,1+X:c+X)=L;
output2(1:r,(c+1):(2*c))=R;
bg1(1+Y:r+Y,1+X:c+X)=bgL;
bg2(1:r,(c+1):(2*c))=bgR;



end
        if Y<0 && X<=c
            Y=-Y;
            output1((r+1-Y):(2*r-Y),(1+X):(c+X))=L;
            output2((1+r):(2*r),(c+1):(2*c))=R;
            bg1((r+1-Y):(2*r-Y),(1+X):(c+X))=bgL;
            bg2((1+r):(2*r),(c+1):(2*c))=bgR;

        end
        if Y>=0 && X>c

            output1(1+Y:r+Y,X-c+1:X)=L;
            output2(1:r,1:c)=R;
            bg1(1+Y:r+Y,X-c+1:X)=bgL;
            bg2(1:r,1:c)=bgR;

        end
        if Y<0 && X>c
            Y=-Y;
            output1((r+1-Y):(2*r-Y),X-c+1:X)=L;
            output2((r+1):(2*r),1:c)=R;
            bg1((r+1-Y):(2*r-Y),X-c+1:X)=bgL;
            bg2((r+1):(2*r),1:c)=bgR;

        end

        fuse=zeros(2*r,2*c);
        fuse_bg=uint8(zeros(2*r,2*c));
        difference=abs(output1-output2);

        difference_no_0=difference(find(difference~=0));%remove the 0 lvl
        totaldifference_ref=sum(difference_no_0(:));
        record_difference=difference_no_0;
        record_output1=output1;
        record_output2=output2;
        
for searchX=max([X-windowsize 1]):1:min([X+windowsize 2*r])
    for searchY=max([Y-windowsize 1]):1:min([Y+windowsize 2*c])
        

        output1=zeros(2*r,2*c);
        output2=zeros(2*r,2*c);
        
        if Y>=0 && X<=c

            output1(1+Y:r+Y,1+X:c+X)=L;
            output2(1:r,(c+1):(2*c))=R;

        end
        if Y<0 && X<=c
            Y=-Y;
            output1((r+1-Y):(2*r-Y),(1+X):(c+X))=L;
            output2((1+r):(2*r),(c+1):(2*c))=R;

        end
        if Y>=0 && X>c

            output1(1+Y:r+Y,X-c+1:X)=L;
            output2(1:r,1:c)=R;

        end
        if Y<0 && X>c
            Y=-Y;
            output1((r+1-Y):(2*r-Y),X-c+1:X)=L;
            output2((r+1):(2*r),1:c)=R;

        end

        fuse=zeros(2*r,2*c);

        difference=abs(output1-output2);

        difference_no_0=difference(find(difference~=0));%remove the 0 lvl
        totaldifference_cur=sum(difference_no_0(:));
        if totaldifference_cur<totaldifference_ref
            totaldifference_ref=totaldifference_cur;
            record_difference=difference_no_0;
            record_x=X;
            record_y=y;
            record_output1=output1;
            record_output2=output2;
        end

    end
end


[hist_diff,center]=hist(record_difference(:),30);%30 bins
% figure,hist(record_difference(:),30);
peaks=findpeaks(hist_diff);% the maximum peaks refers the heigh level difference in overlapping areas( assume the overlapping area occupy large percentage)
threshold=center(find(hist_diff==max(peaks))+2);%assume the center value of right side 4th bin is the threshold to distinguish overlapping part and sided part
compensation=center(hist_diff==max(peaks));%heigh level difference



indices=find(difference>(compensation-1)&difference<(compensation+1));
sum_1=mean(output1(indices(100:200)));
sum_2=mean(output2(indices(100:200)));



for i=1:2*r
    for k=1:2*c
        if(difference(i,k)<threshold)% overlapping part uses avg to estimate the height
            fuse(i,k)=(record_output1(i,k)+record_output2(i,k))/2;
            fuse_bg(i,k)=(bg1(i,k)+bg2(i,k))/2;
        else if output1(i,k)>output2(i,k)%sided part can be estimated by either model
                if sum_1>sum_2
                    fuse(i,k)=record_output1(i,k)-compensation/2;%compensate for the level difference of 2 models
                    fuse_bg(i,k)=bg1(i,k);
                else
                    fuse(i,k)=record_output1(i,k)+compensation/2;
                    fuse_bg(i,k)=bg1(i,k);
                end
            elseif sum_1>sum_2             
                fuse(i,k)=record_output2(i,k)+compensation/2;
                fuse_bg(i,k)=bg2(i,k);
            else
                fuse(i,k)=record_output2(i,k)-compensation/2;
                fuse_bg(i,k)=bg2(i,k);
            end
           
        end
    end
end
 fuse_BG(:,:,1)=fuse_bg;
fuse_BG(:,:,2)=fuse_bg;
fuse_BG(:,:,3)=fuse_bg;
% figure,imagesc(record_output1);
% figure,imagesc(record_output2);
% figure,mesh(fuse);
% figure,imagesc(fuse);
figure,imagesc(difference);
figure,warp(fuse,fuse_BG);



end

