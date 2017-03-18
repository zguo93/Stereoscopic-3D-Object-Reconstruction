function [ bestidx bestInNum bestangle ] = Line_ransac( angle,amp,num,iter,threshDist,threshDist_amp,inlierRatio)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
number = size(angle,2); % Total number of points
bestidx=[];
bestInNum=0;
bestangle=0;
 for i=1:iter
 %% Randomly select 2 points
     idx = randperm(number);   
 %% Compute the distances between all points with the fitting line 
     avg=sum(angle(idx(1:num))/num);
     distance=angle-avg;
     avg_amp=sum(amp(idx(1:num))/num);
     distance_amp=amp-avg_amp;
 %% Compute the inliers with distances smaller than the threshold
%      inlierIdx = find(abs(distance)<=threshDist && abs(distance_amp)<=threshDist_amp);
     inlierIdx=[];
     [r,c]=size(angle);
     for i=1:c
         
         if abs(distance(i))<=threshDist && abs(distance_amp(i))<=threshDist_amp %%check the k and amplitude of  vector
             inlierIdx=[inlierIdx i];
         end
     end
     
     
     inlierNum = length(inlierIdx);
 %% Update the number of inliers and fitting model if better model is found     
     if inlierNum>=round(inlierRatio*number) && inlierNum>bestInNum
         bestInNum = inlierNum;
         bestidx=inlierIdx;
         bestangle=avg;
     end
 end

end

