function [bestmatrix BestDistance bestInIdx bestInNum avX avY angle amp] = ransac_demo(matchpoints,num,iter,threshDist,inlierRatio,reference,current)
 % data: a 2xn dataset with #n data points
 % num: the minimum number of points. For line fitting problem, num=2
 % iter: the number of iterations
 % threshDist: the threshold of the distances between points and the fitting line
 % inlierRatio: the threshold of the numer of inliers 
% [bestmatrix BestDistance bestInIdx bestInNum avX avY angle amp] = ransac_demo(matchpoints,5,1000,80,0.35,reference,current);
 %% Plot the data points
bestInIdx=[];
 number = size(matchpoints,1); % Total number of points
 bestInNum = 0; % Best fitting line with largest number of inliers
bestmatrix=[];% transform ref to cur. ref will also DISPLAY on the left plane in ploting.
 offset_x_current = size(reference.gs,2);                
          combinedimage = zeros(max(size(reference.gs,1),size(current.gs,1)),size(reference.gs,2)+size(current.gs,2));
          combinedimage(1:size(reference.gs,1),1:size(reference.gs,2)) = reference.gs;
          combinedimage(1:size(current.gs,1),(1:size(current.gs,2))+offset_x_current) = current.gs;
 for i=1:iter
 %% Randomly select num points
    
     idx = randperm(number);   
     
     
     
    
 %% Compute the distances between all points with the fitting line 
     src_points=[];
    target_points=[];
    for i=1:num
        src_points=[src_points' [matchpoints(idx(i)).x_ref matchpoints(idx(i)).y_ref]']';
        target_points=[target_points' [matchpoints(idx(i)).x_cur matchpoints(idx(i)).y_cur]']';
    end
    t = cp2tform(src_points, target_points, 'affine');
    distance=[];
    for k=1:number
        distance(k) = norm([matchpoints(k).x_ref matchpoints(k).y_ref 1] * t.tdata.T-[matchpoints(k).x_cur matchpoints(k).y_cur 1] * t.tdata.T); % q = [u ,v, 1]
        
    end
    

 %% Compute the inliers with distances smaller than the threshold

     inlierIdx = find(abs(distance)<=threshDist);
     inlierNum = length(inlierIdx);
 %% Update the number of inliers and fitting model if better model is found     
     if inlierNum>=round(inlierRatio*number) && inlierNum>bestInNum
         bestInNum = inlierNum;
         bestInIdx=inlierIdx;
         BestDistance=distance;
         
          offset_x_current = size(reference.gs,2);                
          combinedimage = zeros(max(size(reference.gs,1),size(current.gs,1)),size(reference.gs,2)+size(current.gs,2));
          combinedimage(1:size(reference.gs,1),1:size(reference.gs,2)) = reference.gs;
          combinedimage(1:size(current.gs,1),(1:size(current.gs,2))+offset_x_current) = current.gs;
                
          % Plot
%           figure, imshow(combinedimage,[]);
%           for i = inlierIdx(1:inlierNum)
%               hline=line([matchpoints(i).x_ref matchpoints(i).x_cur+offset_x_current],[matchpoints(i).y_ref matchpoints(i).y_cur]);                    
%           end      
        
         bestmatrix=t.tdata.T;

     end
 end
 figure, imshow(combinedimage,[]);
 x_axis=0;
 y_axis=0;
 angle=[];
 amp=[];
          for i = bestInIdx(1:bestInNum)
              hline=line([matchpoints(i).x_ref matchpoints(i).x_cur+offset_x_current],[matchpoints(i).y_ref matchpoints(i).y_cur]);
              x_axis=x_axis+(matchpoints(i).x_cur+offset_x_current-matchpoints(i).x_ref);
              y_axis=y_axis+(matchpoints(i).y_cur-matchpoints(i).y_ref);
              amp=[amp sqrt((matchpoints(i).x_cur+offset_x_current-matchpoints(i).x_ref)^2+(matchpoints(i).y_cur-matchpoints(i).y_ref)^2)];
              angle=[angle (matchpoints(i).y_cur-matchpoints(i).y_ref)/(matchpoints(i).x_cur+offset_x_current-matchpoints(i).x_ref)];
          end
          avX=x_axis/bestInNum;
          avY=y_axis/bestInNum;
 %% Plot the best fitting line
 
 [ Line_bestidx Line_bestInNum Line_bestangle ] = Line_ransac( angle,amp,3,1500,0.02,50,0.1);
 
allowed_index= bestInIdx(Line_bestidx);
allowed_index_num=size(allowed_index,2);
figure, imshow(combinedimage,[]);
x_axis=0;
 y_axis=0;
 angle=[];
 amp=[];
for i = allowed_index(1:allowed_index_num)
              hline=line([matchpoints(i).x_ref matchpoints(i).x_cur+offset_x_current],[matchpoints(i).y_ref matchpoints(i).y_cur]);
              x_axis=[x_axis (matchpoints(i).x_cur+offset_x_current-matchpoints(i).x_ref)];
              y_axis=[y_axis (matchpoints(i).y_cur-matchpoints(i).y_ref)];
              angle=[angle (matchpoints(i).y_cur-matchpoints(i).y_ref)/(matchpoints(i).x_cur+offset_x_current-matchpoints(i).x_ref)];
              amp=[amp sqrt((matchpoints(i).x_cur+offset_x_current-matchpoints(i).x_ref)^2+(matchpoints(i).y_cur-matchpoints(i).y_ref)^2)];
end
avX=sum(x_axis)/allowed_index_num;
avY=sum(y_axis)/allowed_index_num;










end





