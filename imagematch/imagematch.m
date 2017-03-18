classdef imagematch < handle
% This is a matlab class which performs the SIFT (scale invariant feature
% transform) on two images (called the reference and current images),
% matches similar feature points, and then peforms a least square pose
% estimation to determine where the reference image lies within the current
% image. This is a side project I worked on and is far from a
% proper/efficient/robust implementation. I intend it to be a nice
% introductory learning resource for people who have heard of SIFT and 
% would like to see the basics of how it works and tinker with it within 
% the matlab environment. 
%
%
% Anyway, the most important methods are:
% 
% 1) set_image(data) - This sets the image for the imagematch object to use
% 2) get_featurepoints() - This uses the DOG approach to find suitable
%    featurepoints. This is described in David Lowe's paper on SIFT. It
%    also uses "affine invariance" which I found on slides 23-24 here:
%    http://www.cs.illinois.edu/~slazebni/spring11/lec08_blob.pdf
%    when calculating the "sift patches." Lastly, orientations (possibly 
%    more than one) are recorded using a histogram approach.
% 3) get_descriptors() - This method calculates the "sift vectors" from the
%    sift patches for each possible orientation
% 4) get_matchpoints(reference,current) - This static method matches
%    featurepoints that have similar sift vectors.
% 5) get_pose(reference,current,matchpoints) - This static method performs
%    a least squares pose estimation based on the matched points
%
% 
% Lastly, there are many plotting methods (which being with "plot" - i.e. 
% plot_featurepoints()) to demonstrate what is going on in each step.
%
    % Properties ---------------------------------------------------------%
    properties(SetAccess = private)
        gs;                 % double    
        gaussianpyramid;    % cell
        logpyramid;         % cell
        gradpyramid;        % cell
        
        maximaplots;        % cell
        featurepoints;      % struct
        descriptors;        % struct
    end
   
    % Methods ------------------------------------------------------------%
    methods(Access = public)
        % Constructor
        function obj = imagematch   
            obj.gs = [];               
            obj.gaussianpyramid = cell(0);
            obj.logpyramid = cell(0); 
            obj.gradpyramid = cell(0); 
            obj.maximaplots = cell(0);
            obj.featurepoints = struct('x',{},'y',{},'scale',{},'ellipse_affine',{},'orientation',struct('patch_sift',{},'mat_trans',{}));
            obj.descriptors = struct('orientation',struct('patch_sift',{})); 
        end
        
        function set_image(obj,image_i)        
            if(isreal(image_i))
                if((strcmp(class(image_i),'uint8') || strcmp(class(image_i),'uint16') || strcmp(class(image_i),'double')) && ...
                  (length(size(image_i)) == 2 || length(size(image_i)) == 3) && (size(image_i,3) == 1 || size(image_i,3) == 3) && ...
                  size(image_i,1) >= 16 && size(image_i,2) >= 16)
                    % Format grayscale values
                    if(size(image_i,3) == 3)
                        % Image is RGB-color
                        % Set gs:
                        obj.gs = im2double(rgb2gray(image_i));
                    else
                        % Image is monochrome
                        % Set gs:
                        obj.gs = im2double(image_i);
                    end                   
                    
                    % Get gaussian/LOG/grad pyramid
                    obj.set_pyramids();       
                else
                    error('Image input must be mxn or mxnx3 and must be uint8, uint16 or double. It also must be greater than or equal to 16x16 in resolution.');
                end
            else
                error('Image must be real.');
            end
        end
        
        function get_featurepoints(obj)                
            if(~isempty(obj.gs))  
                % Calculates maximaplots and featurepoints
                obj.maximaplots = cell(0);
                obj.featurepoints(:) = [];
                for i = 0:3 % Cycle over octaves - 4 total                    
                    for j = 1:2 % Cycle middle two LOG scales
                        % Cycle inner portion of of LOG plane - take absolute
                        % value to limit search to maxima
                        scale1 = abs(obj.logpyramid{i+1,j});
                        scale2 = abs(obj.logpyramid{i+1,j+1});                  
                        scale3 = abs(obj.logpyramid{i+1,j+2});

                        % Get boolean plot that shows which points are maxima
                        combscale = zeros(size(scale1,1),size(scale1,2),3);
                        combscale(:,:,1) = scale1;
                        combscale(:,:,2) = scale2;
                        combscale(:,:,3) = scale3;
                        
                        msk = true(3,3,3);
                        msk(2,2,2) = false;
                        % Assign, to every voxel, the maximum of its
                        % neighbors - From Jonas on stackoverflow
                        combscale_dil = imdilate(combscale,msk);
                        localmax = combscale > combscale_dil;
                        
                        % Get Plot
                        localmax = localmax(:,:,2);
                        localmax(1,:) = false;
                        localmax(end,:) = false;
                        localmax(:,1) = false;
                        localmax(:,end) = false;
                        
                        % Store plot
                        obj.maximaplots{i+1,j} = localmax;                        

                        if(any(obj.maximaplots{i+1,j}(:)))
                            % Get locations of maxima
                            [y_maxima x_maxima] = find(obj.maximaplots{i+1,j});
                            % Convert to 0-based indexing
                            y_maxima = y_maxima-1;
                            x_maxima = x_maxima-1;

                            % Go through points
                            for k = 0:length(y_maxima)-1                                                                                                
                                % Get subpixel location - only do one iteration
                                % using gradient information
                                sigma_0 = (2^i)*(sqrt(2)^(j-2));
                                sigma_1 = (2^i)*(sqrt(2)^(j-1));
                                sigma_2 = (2^i)*(sqrt(2)^(j));
                                dsigma_0 = sigma_1-sigma_0; 
                                dsigma_1 = sigma_2-sigma_1;

                                % dx and dy are both 1

                                % Calculate Gradient - [dD/dx dD/dy dD/dsigma]
                                dg_dp(1) = (obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)+2)-obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)))/2;
                                dg_dp(2) = (obj.logpyramid{i+1,j+1}(y_maxima(k+1)+2,x_maxima(k+1)+1)-obj.logpyramid{i+1,j+1}(y_maxima(k+1),x_maxima(k+1)+1))/2;                            
                                dg_dp(3) = ((obj.logpyramid{i+1,j+2}(y_maxima(k+1)+1,x_maxima(k+1)+1)-obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)+1))/dsigma_1 + ...
                                           (obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)+1)-obj.logpyramid{i+1,j}(y_maxima(k+1)+1,x_maxima(k+1)+1))/dsigma_0)/2;

                                % Calculate Hessian - only upper half
                                ddg_dp2(1,1) = obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)+2) - ...
                                               2*obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)+1) + ...
                                               obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)); %dxdx
                                ddg_dp2(1,2) = ((obj.logpyramid{i+1,j+1}(y_maxima(k+1)+2,x_maxima(k+1)+2) - obj.logpyramid{i+1,j+1}(y_maxima(k+1),x_maxima(k+1)+2))/2 - ...
                                               (obj.logpyramid{i+1,j+1}(y_maxima(k+1)+2,x_maxima(k+1)) - obj.logpyramid{i+1,j+1}(y_maxima(k+1),x_maxima(k+1)))/2)/2; %dydx
                                ddg_dp2(2,2) = obj.logpyramid{i+1,j+1}(y_maxima(k+1)+2,x_maxima(k+1)+1) - ...
                                               2*obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)+1) + ...
                                               obj.logpyramid{i+1,j+1}(y_maxima(k+1),x_maxima(k+1)+1); %dydy
                                ddg_dp2(1,3) = (((obj.logpyramid{i+1,j+2}(y_maxima(k+1)+1,x_maxima(k+1)+2) - obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)+2))/dsigma_1 + ...
                                               (obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)+2) - obj.logpyramid{i+1,j}(y_maxima(k+1)+1,x_maxima(k+1)+2))/dsigma_0)/2 - ...
                                               ((obj.logpyramid{i+1,j+2}(y_maxima(k+1)+1,x_maxima(k+1)) - obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)))/dsigma_1 + ...
                                               (obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)) - obj.logpyramid{i+1,j}(y_maxima(k+1)+1,x_maxima(k+1)))/dsigma_0)/2)/2; %dsigmadx
                                ddg_dp2(2,3) = (((obj.logpyramid{i+1,j+2}(y_maxima(k+1)+2,x_maxima(k+1)+1) - obj.logpyramid{i+1,j+1}(y_maxima(k+1)+2,x_maxima(k+1)+1))/dsigma_1 + ...
                                               (obj.logpyramid{i+1,j+1}(y_maxima(k+1)+2,x_maxima(k+1)+1) - obj.logpyramid{i+1,j}(y_maxima(k+1)+2,x_maxima(k+1)+1))/dsigma_0)/2 - ...
                                               ((obj.logpyramid{i+1,j+2}(y_maxima(k+1),x_maxima(k+1)+1) - obj.logpyramid{i+1,j+1}(y_maxima(k+1),x_maxima(k+1)+1))/dsigma_1 + ...
                                               (obj.logpyramid{i+1,j+1}(y_maxima(k+1),x_maxima(k+1)+1) - obj.logpyramid{i+1,j}(y_maxima(k+1),x_maxima(k+1)+1))/dsigma_0)/2)/2; %dsigmady
                                ddg_dp2(3,3) = ((obj.logpyramid{i+1,j+2}(y_maxima(k+1)+1,x_maxima(k+1)+1) - obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)+1))/dsigma_1 - ...
                                               (obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)+1) - obj.logpyramid{i+1,j}(y_maxima(k+1)+1,x_maxima(k+1)+1))/dsigma_0)/((dsigma_1+dsigma_0)/2); %dsigmadsigma

                                % Fill lower half of hessian
                                ddg_dp2(2,1) = ddg_dp2(1,2);
                                ddg_dp2(3,1) = ddg_dp2(1,3);
                                ddg_dp2(3,2) = ddg_dp2(2,3);      
                                
                                % Find incremental parameters
                                deltap = -ddg_dp2^-1*dg_dp';   

                                x_maxima_subpixel = x_maxima(k+1)+deltap(1);
                                y_maxima_subpixel = y_maxima(k+1)+deltap(2);
                                scale_maxima_subpixel = (2^i)*(sqrt(2)^(j-1)) + deltap(3);  

                                % Check to make sure contrast is high,
                                % coordinates have not been moved out of bounds,
                                % deltap remains at least within the 26
                                % neighbor box - this assumes large changes
                                % imply point is unstable, and use the
                                % hessian to calculate the edge response to
                                % remove points which lie on edges
                                if(abs(deltap(1)) < 1 && abs(deltap(2)) < 1 && abs(deltap(3)) < dsigma_1 &&...
                                   abs(obj.logpyramid{i+1,j+1}(y_maxima(k+1)+1,x_maxima(k+1)+1) + ...
                                       (1/2)*dg_dp*deltap) > 0.03 && ...
                                   x_maxima_subpixel >= 0 && ...
                                   x_maxima_subpixel < size(obj.logpyramid{i+1,j+1},2)-1 && ...
                                   y_maxima_subpixel >= 0 && ...
                                   y_maxima_subpixel < size(obj.logpyramid{i+1,j+1},1)-1 && ...
                                   scale_maxima_subpixel > 0 && ...
                                   trace(ddg_dp2(1:2,1:2))^2/det(ddg_dp2(1:2,1:2)) < 12.1) % Less than (r+1)^2/r where r is set to 10
                                    % Now perform affine invariance;
                                    % use bilinear interpolation to
                                    % interpolate gradients at center; do 
                                    % not use subpixel scale; only subpixel
                                    % x and y      
                                    radius_maxima_octavescaled = ceil(sqrt(2)*scale_maxima_subpixel/(2^i)); % Get radius from sigma conversion, and then account for octave
                                    side_kernel = 6*radius_maxima_octavescaled+1; % Must be add to be centered around maxima
                                    halfside_kernel = (side_kernel-1)/2;
                                                                          
                                    % Obtain local image gradients for
                                    % affine invariance
                                    kernel_gauss = fspecial('gaussian',[side_kernel side_kernel],1.5*scale_maxima_subpixel/(2^i));
                                    M = zeros(2);

                                    % Debug ------------------------------%
                                    % debug_gradx = zeros(size(kernel_gauss));
                                    % debug_grady = zeros(size(kernel_gauss));
                                    % ------------------------------------%
                                    
                                    for l = x_maxima(k+1)-(halfside_kernel):x_maxima(k+1)+(halfside_kernel)
                                        x_ref = l;
                                        x_transformed = x_ref+deltap(1);
                                        for m = y_maxima(k+1)-(halfside_kernel):y_maxima(k+1)+(halfside_kernel)
                                            y_ref = m;
                                            y_transformed = y_ref+deltap(2);
                                            
                                            if(x_transformed < 1 || y_transformed < 1 || x_transformed > size(obj.gaussianpyramid{i+1,j+1},2)-2 || y_transformed > size(obj.gaussianpyramid{i+1,j+1},1)-2)
                                                % Coordinate are out of bounds, skip
                                                % (this will set the value to zero)
                                                continue;
                                            end   

                                            y_transformed_floor = floor(y_transformed);
                                            x_transformed_floor = floor(x_transformed);

                                            y_transformed_tilda = y_transformed-y_transformed_floor;
                                            x_transformed_tilda = x_transformed-x_transformed_floor;

                                            % Interpolate gradients -
                                            % use scale at {i+1,j+1}
                                            val_kernel = kernel_gauss(y_ref-(y_maxima(k+1)-(halfside_kernel))+1,x_ref-(x_maxima(k+1)-(halfside_kernel))+1);
                                            grad_x = val_kernel*[1-y_transformed_tilda y_transformed_tilda]*obj.gradpyramid{i+1,j+1}.x(y_transformed_floor+1:y_transformed_floor+2,x_transformed_floor+1:x_transformed_floor+2)*[1-x_transformed_tilda x_transformed_tilda]';
                                            grad_y = val_kernel*[1-y_transformed_tilda y_transformed_tilda]*obj.gradpyramid{i+1,j+1}.y(y_transformed_floor+1:y_transformed_floor+2,x_transformed_floor+1:x_transformed_floor+2)*[1-x_transformed_tilda x_transformed_tilda]';

                                            % Debug ----------------------%
                                            % debug_gradx(m-(y_maxima(k+1)-(halfside_kernel))+1,l-(x_maxima(k+1)-(halfside_kernel))+1) = grad_x;
                                            % debug_grady(m-(y_maxima(k+1)-(halfside_kernel))+1,l-(x_maxima(k+1)-(halfside_kernel))+1) = grad_y;
                                            % ----------------------------%
                                            
                                            M(1,1) = M(1,1) + grad_x^2;
                                            M(1,2) = M(1,2) + grad_x*grad_y;
                                            M(2,2) = M(2,2) + grad_y^2; 
                                        end
                                    end
                                    % Finish M
                                    M(2,1) = M(1,2);   

                                    % Affine Invariance
                                    [V D] = eig(M);
                                    % Format V to turn it into rotation
                                    % matrix
                                    V = flipud(V);
                                    D = rot90(D,2);                          

                                    eig_1 = D(1,1);
                                    eig_2 = D(2,2);

                                    % Sample and normalize image patch
                                    % Get new transformation matrix
                                    mat_dilate = [2*sqrt(eig_2)/(sqrt(eig_1)+sqrt(eig_2)) 0; 
                                                  0                                       2*sqrt(eig_1)/(sqrt(eig_1)+sqrt(eig_2))];
                                    mat_trans = V^-1*mat_dilate; % rotation * dilation

                                    % Initialize normalized gradients for
                                    % orientation calculation - note
                                    % that these cant be used for sift
                                    % descriptors because orientation
                                    % has not been set yet
                                    normalizedpatch = zeros(side_kernel);
                                    for l = x_maxima(k+1)-(halfside_kernel):x_maxima(k+1)+(halfside_kernel)
                                        x_ref = l;
                                        x_transformed = x_ref+deltap(1);
                                        for m = y_maxima(k+1)-(halfside_kernel):y_maxima(k+1)+(halfside_kernel)
                                            y_ref = m;
                                            y_transformed = y_ref+deltap(2);                                                                                                

                                            % Transform coordinates
                                            y_transformed_o = y_maxima_subpixel + mat_trans(2,1)*(x_transformed-x_maxima_subpixel) + mat_trans(2,2)*(y_transformed-y_maxima_subpixel); 
                                            x_transformed_o = x_maxima_subpixel + mat_trans(1,1)*(x_transformed-x_maxima_subpixel) + mat_trans(1,2)*(y_transformed-y_maxima_subpixel);

                                            if(x_transformed_o < 1 || y_transformed_o < 1 || x_transformed_o > size(obj.gaussianpyramid{i+1,j+1},2)-2 || y_transformed_o > size(obj.gaussianpyramid{i+1,j+1},1)-2)
                                                % Coordinate are out of bounds, skip
                                                % (this will set the value to zero)
                                                continue;
                                            end   

                                            y_transformed_floor = floor(y_transformed_o);
                                            x_transformed_floor = floor(x_transformed_o);

                                            y_transformed_tilda = y_transformed_o-y_transformed_floor;
                                            x_transformed_tilda = x_transformed_o-x_transformed_floor;

                                            % Interpolate                                                
                                            normalizedpatch(m-(y_maxima(k+1)-(halfside_kernel))+1,l-(x_maxima(k+1)-(halfside_kernel))+1) = ... 
                                                [1-y_transformed_tilda y_transformed_tilda]*obj.gaussianpyramid{i+1,j+1}(y_transformed_floor+1:y_transformed_floor+2,x_transformed_floor+1:x_transformed_floor+2)*[1-x_transformed_tilda x_transformed_tilda]';
                                        end
                                    end

                                    % Debug ------------------------------%
                                    %{
                                    figure;
                                    subplot(1,2,1)                                        
                                    imshow(normalizedpatch,[]); 
                                    subplot(1,2,2) 
                                    halfside_kernel = 5;
                                    imshow(obj.gaussianpyramid{i+1,j+1}((y_maxima(k+1)-(halfside_kernel):y_maxima(k+1)+(halfside_kernel))+1,(x_maxima(k+1)-(halfside_kernel):x_maxima(k+1)+(halfside_kernel))+1),[])
                                    %}
                                    % ------------------------------------%

                                    % Now perform orientation assignment on
                                    % normalized patch. Take gradient by
                                    % using finite difference (this leaves
                                    % 1 pixel boundary of empty values
                                    
                                    % Get histogram - 36 bins
                                    gradspacing = linspace(-pi,pi,37);
                                    gradhist = zeros(1,length(gradspacing)-1);

                                    % Use finite difference to calculate
                                    % gradient of normalized patch. These
                                    % are weighted by gaussian kernel. Also
                                    % leave a 1 pixel border.
                                    for l = 1:side_kernel-2
                                        for m = 1:side_kernel-2
                                            normalizedpatch_dx = (normalizedpatch(m+1,l+2)-normalizedpatch(m+1,l))/2; 
                                            normalizedpatch_dy = (normalizedpatch(m+2,l+1)-normalizedpatch(m,l+1))/2;
                                            
                                            % Get angle
                                            gradangle = atan2(normalizedpatch_dy,normalizedpatch_dx);
                                            
                                            % Get gaussian weighted
                                            % magnitude
                                            val_kernel = kernel_gauss(m+1,l+1);
                                            gradmag = val_kernel*sqrt(normalizedpatch_dx^2 + normalizedpatch_dy^2);

                                            % Store in histogram
                                            if(gradangle == pi)
                                                % This is a special case
                                                gradhist(end) = gradhist(end) + gradmag;                                     
                                            else
                                                % This only works for
                                                % evenly spaced bins
                                                bin = floor((gradangle+pi)/(2*pi/length(gradhist)));
                                                gradhist(bin+1) = gradhist(bin+1) + gradmag;                                        
                                            end                                            
                                        end
                                    end
                                                                        
                                    % Find subpixel maximum, use 3
                                    % neighbors and circular. 
                                    maximaloc = [];
                                    % Edge case
                                    if(gradhist(1) > gradhist(end) && gradhist(1) > gradhist(2))    
                                        maximaloc(end+1) = 0;
                                    end
                                    % Middle
                                    for l = 1:length(gradhist)-2    
                                        if(gradhist(l+1) > gradhist(l) && gradhist(l+1) > gradhist(l+2))
                                            maximaloc(end+1) = l;                       
                                        end                                  
                                    end
                                    % Edge case
                                    if(gradhist(end) > gradhist(end-1) && gradhist(end) > gradhist(1))
                                        maximaloc(end+1) = length(gradhist)-1;
                                    end
                                    if(isempty(maximaloc))
                                        % This can happen is the gradhist
                                        % is flat (i.e. there are no
                                        % maxima)
                                        maximaloc = 0;
                                    end

                                    % Only keep maximalocations greater than 80% of
                                    % the max
                                    maximaloc_filtered = maximaloc(gradhist(maximaloc+1) > 0.8*max(gradhist));                                    
                                    
                                    % Cycle over maximaloc_filtered - may
                                    % create additional feature points
                                    orientation = struct('patch_sift',{},'mat_trans',{}); 
                                    for l = 0:length(maximaloc_filtered)-1
                                        % Interpolate value
                                        if(maximaloc_filtered(l+1) == 0)
                                            maximaloc_filtered(l+1) = maximaloc_filtered(l+1) - ((gradhist(2)-gradhist(end))/2)/(gradhist(2)-2*gradhist(1)+gradhist(end));
                                        elseif(maximaloc_filtered(l+1) == length(gradhist)-1)
                                            maximaloc_filtered(l+1) = maximaloc_filtered(l+1) - ((gradhist(1)-gradhist(end-1))/2)/(gradhist(1)-2*gradhist(end)+gradhist(end-1));
                                        else
                                            maximaloc_filtered(l+1) = maximaloc_filtered(l+1) - ((gradhist(maximaloc_filtered(l+1)+2)-gradhist(maximaloc_filtered(l+1)))/2)/(gradhist(maximaloc_filtered(l+1)+2)-2*gradhist(maximaloc_filtered(l+1)+1)+gradhist(maximaloc_filtered(l+1)));
                                        end
                                        % Make sure maximaloc_filtered did not
                                        % overshoot bounds
                                        maximaloc_filtered(l+1) = min(max(maximaloc_filtered(l+1),-1),length(gradhist));

                                        % maximaloc_filtered is subpixel
                                        % now; obtain the angle from the
                                        % location
                                        gradangles = linspace(-((length(gradhist)-1)/length(gradhist))*pi,((length(gradhist)-1)/length(gradhist))*pi,length(gradhist));                                    
                                        maximaloc_floor = floor(maximaloc_filtered(l+1));
                                        maximaloc_delta = maximaloc_filtered(l+1) - floor(maximaloc_filtered(l+1));
                                        if((maximaloc_filtered(l+1) >= -1 && maximaloc_filtered(l+1) < 0) || (maximaloc_filtered(l+1) >= length(gradhist)-1 && maximaloc_filtered(l+1) <= length(gradhist)))
                                            normalizationangle = gradangles(end)+(gradangles(1)-gradangles(end-1))/2*maximaloc_delta+1/2*(gradangles(1)-2*gradangles(end)+gradangles(end-1))*maximaloc_delta^2;
                                        elseif(maximaloc_filtered(l+1) >= 0 && maximaloc_filtered(l+1) < 1)
                                            normalizationangle = gradangles(1)+(gradangles(2)-gradangles(end))/2*maximaloc_delta+1/2*(gradangles(2)-2*gradangles(1)+gradangles(end))*maximaloc_delta^2;
                                        else
                                            normalizationangle = gradangles(maximaloc_floor+1)+(gradangles(maximaloc_floor+2)-gradangles(maximaloc_floor))/2*maximaloc_delta+1/2*(gradangles(maximaloc_floor+2)-2*gradangles(maximaloc_floor+1)+gradangles(maximaloc_floor))*maximaloc_delta^2;
                                        end                                    

                                        % Convert normalization angle to matrix tranformation  
                                        mat_normalizationangle = [cos(normalizationangle) -sin(normalizationangle);
                                                                  sin(normalizationangle) cos(normalizationangle)];

                                        % Reform transformation matrix
                                        mat_trans_rotinvar = V^-1*mat_dilate*mat_normalizationangle; % rotation * dilation * rotation                                        
                                        
                                        % Now sample normalized patch; this
                                        % will be used to calculate the gradients
                                        % to form the sift descriptors.
                                        % Also interpolate grayscale values. Use
                                        % a 18x18 window (for 1 pixel
                                        % border)
                                        width_siftwindow = 2*(2*sqrt(2)*scale_maxima_subpixel/(2^i)); % Twice 3*octavenormalized radius
                                        width_18_half = (18-1)/2;
                                        patch_sift = zeros(18);
                                        
                                        % Debug --------------------------%
                                        % debug_gs = obj.gaussianpyramid{i+1,j+1};
                                        % --------------------------------%
                                        
                                        for m = 0:17
                                            x_delta = ((m-width_18_half)/17)*width_siftwindow;
                                            for n = 0:17
                                                y_delta = ((n-width_18_half)/17)*width_siftwindow;
                                                % Apply transformation to
                                                % x_delta and y_delta then
                                                % shift by the x and y
                                                % subpixel maxima locations
                                                y_transformed = y_maxima_subpixel + mat_trans_rotinvar(2,1)*x_delta + mat_trans_rotinvar(2,2)*y_delta; 
                                                x_transformed = x_maxima_subpixel + mat_trans_rotinvar(1,1)*x_delta + mat_trans_rotinvar(1,2)*y_delta;
                                                
                                                if(x_transformed < 1 || y_transformed < 1 || x_transformed > size(obj.gaussianpyramid{i+1,j+1},2)-2 || y_transformed > size(obj.gaussianpyramid{i+1,j+1},1)-2)
                                                    % Coordinate are out of bounds, skip
                                                    % (this will set the value to zero)
                                                    continue;
                                                end   
                                                
                                                % Debug ------------------%
                                                % debug_gs(floor(y_transformed)+1,floor(x_transformed)+1) = 1;
                                                % ------------------------%                                                
                                                
                                                y_transformed_floor = floor(y_transformed);
                                                x_transformed_floor = floor(x_transformed);

                                                y_transformed_tilda = y_transformed-y_transformed_floor;
                                                x_transformed_tilda = x_transformed-x_transformed_floor;

                                                % Interpolate patch -
                                                % use scale at {i+1,j+1}
                                                patch_sift(n+1,m+1) = [1-y_transformed_tilda y_transformed_tilda]*obj.gaussianpyramid{i+1,j+1}(y_transformed_floor+1:y_transformed_floor+2,x_transformed_floor+1:x_transformed_floor+2)*[1-x_transformed_tilda x_transformed_tilda]';
                                            end
                                        end
                                        
                                        % Debug --------------------------%
                                        % figure, imshow(debug_gs,[]);
                                        % --------------------------------%
                                                                                
                                        % This is a good featurepoint;
                                        % append it
                                        x_maxima_ref = min(max(round((x_maxima_subpixel)*2^(i)),0),size(obj.gs,2)-1);
                                        y_maxima_ref = min(max(round((y_maxima_subpixel)*2^(i)),0),size(obj.gs,1)-1);          
                                        radius_maxima_ref = ceil(sqrt(2)*scale_maxima_subpixel);
                                        scaledaxes = 2*radius_maxima_ref/(sqrt(eig_1)+sqrt(eig_2));

                                        orientation_template.patch_sift = patch_sift;
                                        orientation_template.mat_trans = mat_trans_rotinvar;
                                        
                                        % Store
                                        orientation = vertcat(orientation,orientation_template);
                                    end
                                    
                                    featurepoints_template.x = x_maxima_ref;
                                    featurepoints_template.y = y_maxima_ref;
                                    featurepoints_template.scale = scale_maxima_subpixel;
                                    featurepoints_template.ellipse_affine = [sqrt(eig_1)*scaledaxes sqrt(eig_2)*scaledaxes atan2(V(1,1),V(2,1))]; % [axis1 axis2 angle]
                                    featurepoints_template.orientation = orientation;

                                    % Store
                                    obj.featurepoints = vertcat(obj.featurepoints,featurepoints_template);  
                                end
                            end
                        end
                    end
                end 
            else
                error('Image has not been set yet.');
            end
        end
        
        function get_descriptors(obj)
            if(~isempty(obj.featurepoints))
                obj.descriptors(:) = [];
                % Cycle over sift patches and calculate the sift
                % descriptors. 
                kernel_gauss_16 = fspecial('gaussian',[16 16],8);                
                for i = 0:length(obj.featurepoints)-1
                    orientation = struct('sift_vec',{});
                    for j = 0:length(obj.featurepoints(i+1).orientation)-1
                        grad_x = zeros(16);
                        grad_y = zeros(16);
                        for k = 1:16
                            for l = 1:16
                                val_kernel = kernel_gauss_16(l,k);
                                grad_x(l,k) = val_kernel*(obj.featurepoints(i+1).orientation(j+1).patch_sift(l+1,k+2)-obj.featurepoints(i+1).orientation(j+1).patch_sift(l+1,k))/2;
                                grad_y(l,k) = val_kernel*(obj.featurepoints(i+1).orientation(j+1).patch_sift(l+2,k+1)-obj.featurepoints(i+1).orientation(j+1).patch_sift(l,k+1))/2;
                            end
                        end

                        % Now calculate sift descriptor
                        sift_vec = zeros(1,128);
                        for k = 0:3
                            for l = 0:3
                                sift_vec(((l+4*k)*8:(l+4*k+1)*8-1)+1) = sifthist(grad_x((l*4:(l+1)*4-1)+1,(k*4:(k+1)*4-1)+1),...
                                                                                 grad_y((l*4:(l+1)*4-1)+1,(k*4:(k+1)*4-1)+1));
                            end
                        end        

                        % Normalize
                        sift_vec = sift_vec/norm(sift_vec);

                        % Threshold any values above 0.2, then renormalize
                        sift_vec(sift_vec > 0.2) = 0.2;
                        sift_vec = sift_vec/norm(sift_vec);

                        orientation_template.sift_vec = sift_vec;
                        
                        % Store
                        orientation = vertcat(orientation,orientation_template);
                    end

                    descriptors_template.orientation = orientation;

                    % Store;
                    obj.descriptors = vertcat(obj.descriptors,descriptors_template);
                end            
            else
                error('Feature points have not been determined yet or none were found.');
            end    
            
            function hist_vec = sifthist(grad_x_sub,grad_y_sub)                
                % Get histogram - 8 bins
                % Gradients already gauss weighted
                gradspacing = linspace(-pi,pi,9);
                hist_vec = zeros(1,length(gradspacing)-1);
                
                for m = 0:3
                    for n = 0:3
                        % Get angle
                        gradangle = atan2(grad_y_sub(n+1,m+1),grad_x_sub(n+1,m+1));

                        % Get magnitude
                        gradmag = sqrt(grad_x_sub(n+1,m+1)^2 + grad_y_sub(n+1,m+1)^2);

                        % Store in histogram
                        if(gradangle == pi)
                            % This is a special case
                            hist_vec(end) = hist_vec(end) + gradmag;                                     
                        else
                            % This only works for
                            % evenly spaced bins
                            bin = floor((gradangle+pi)/(2*pi/length(hist_vec)));
                            hist_vec(bin+1) = hist_vec(bin+1) + gradmag;                                        
                        end                                            
                    end
                end                
            end           
        end
        
        function plot_featurepoints(obj)               
            if(~isempty(obj.featurepoints))  
                figure, imshow(obj.gs,[]);
                for i = 0:size(obj.featurepoints,1)-1
                    ellipse(obj.featurepoints(i+1).ellipse_affine(1), ...   % axis2
                            obj.featurepoints(i+1).ellipse_affine(2), ...   % axis1
                            obj.featurepoints(i+1).ellipse_affine(3), ...   % angle
                            obj.featurepoints(i+1).x+1, ...                 % X coordinate, converted to 1 based indexing
                            obj.featurepoints(i+1).y+1, ...                 % Y coordinate, converted to 1 based indexing
                            'r');
                    text(obj.featurepoints(i+1).x+1,obj.featurepoints(i+1).y+1,num2str(i+1));
                end   
            else
                error('Featurepoints have not been calculated yet');
            end
        end
                
        function plot_pyramid(obj)          
            if(~isempty(obj.featurepoints))  
                % Get total width and height    
                width_total = 0;                       
                height_total = 0;             
                for i = 0:3 % Cycle over scales - 4 total
                    width_total = width_total+size(obj.gaussianpyramid{i+1,1},2);
                    height_total = height_total+size(obj.gaussianpyramid{i+1,1},1);
                end

                compositegaussian = zeros(width_total,height_total);   
                ind_x = 0;
                for i = 0:3 % Cycle over octaves - 4 total                   
                    for j = 0:4 % Cycle over scales - 5 total
                        ind_y = j*size(obj.gaussianpyramid{i+1,j+1},1);
                        compositegaussian(ind_y+1:ind_y+size(obj.gaussianpyramid{i+1,j+1},1),ind_x+1:ind_x+size(obj.gaussianpyramid{i+1,j+1},2)) = obj.gaussianpyramid{i+1,j+1};
                    end
                    ind_x = ind_x + size(obj.gaussianpyramid{i+1,1},2); 
                end 

                figure, imshow(compositegaussian,[]);                
            else
                error('Featurepoints have not been calculated yet.');
            end
        end
    end  
    
    methods(Access = private)        
        function set_pyramids(obj)
            if(~isempty(obj.gs))                
                obj.gaussianpyramid = cell(0);
                obj.logpyramid = cell(0);   
                kernel_sobel_x = -fspecial('sobel')'/8;
                kernel_sobel_y = -fspecial('sobel')/8;
                for i = 0:3 % Cycle over octaves - 4 total                    
                    for j = 0:4 % Cycle over scales - 5 total
                        if(i == 0 && j == 0)
                            % This is the first octave and first octave
                            obj.gaussianpyramid{i+1,j+1} = imagematch.filt_gaussian(obj.gs,1/sqrt(2));
                        elseif(i > 0 && j == 0)
                            % This is the first scale after the first
                            % octave, take 3rd image from previous row
                            obj.gaussianpyramid{i+1,j+1} = obj.gaussianpyramid{i,3}(1:2:end,1:2:end);
                        else
                            % These are second + scales. Take image from
                            % previous column and apply the same scale
                            obj.gaussianpyramid{i+1,j+1} = imagematch.filt_gaussian(obj.gaussianpyramid{i+1,j},sqrt(2)^(j-2));   

                            % Store logpyramid
                            obj.logpyramid{i+1,j} = obj.gaussianpyramid{i+1,j+1} - obj.gaussianpyramid{i+1,j};
                        end         
                        
                        obj.gradpyramid{i+1,j+1}.x = imfilter(obj.gaussianpyramid{i+1,j+1},kernel_sobel_x,'same','replicate');
                        obj.gradpyramid{i+1,j+1}.y = imfilter(obj.gaussianpyramid{i+1,j+1},kernel_sobel_y,'same','replicate');
                    end
                end
            else
                error('Image has not been set yet');
            end
        end  
    end
    
    methods(Static)        
        function data_gs_padded_filt = filt_gaussian(data_gs_padded,sigma)              
            if(isnumeric(data_gs_padded) && isnumeric(sigma) && isreal(sigma) && sigma > 0)                 
                % window MUST BE ODD. ceil insures integer, and
                % multiplication of an integer by an even number insures
                % another even number. Adding by 1 insures an odd number.
                window = 6*ceil(sigma)+1;
                data_gs_padded_filt = imfilter(data_gs_padded,fspecial('gaussian',[window window],sigma),'same','replicate');
            else
                error('Inputs are incorrect.');
            end
        end
    
        function matchpoints = get_matchpoints(reference,current)
            if(isa(reference,'imagematch') && isa(current,'imagematch') && length(reference.descriptors) >= 1 && length(current.descriptors) >= 2)
                % Initialize matchpoints
                matchpoints = struct('x_ref',{},'y_ref',{},'descriptorind_ref',{},'orientationind_ref',{}, ...
                                     'x_cur',{},'y_cur',{},'descriptorind_cur',{},'orientationind_cur',{}, ...
                                     'dist',{});      
                
                % Do exhaustive search to find best match. Cycle over 
                % reference descriptors and find the best orientation with
                % the best current descriptor that best matches it. Find the minimum
                % of the squared euclidean distance. Pick a the orientation
                % with the best match, and also check to make sure
                % coordinates for the current point aren't already in the
                % matchpoints list before appending it; if there's a
                % duplicate point then pick the best one, these precautions
                % ensure a one-to-one correspondence
                for i = 0:length(reference.descriptors)-1
                    % Cycle through reference orientations and pick the one
                    % with the best match
                    best_match_overall.descriptorind_ref = i;
                    best_match_overall.orientationind_ref = -1;
                    best_match_overall.descriptorind_cur = -1;
                    best_match_overall.orientationind_cur = -1;
                    best_match_overall.dist = inf;
                    second_match_overall.dist = inf;
                    for j = 0:length(reference.descriptors(i+1).orientation)-1    
                        % Get reference descriptor
                        descriptor_ref = reference.descriptors(i+1).orientation(j+1).sift_vec;

                        % Cycle over all current descriptors and get the
                        % best and second best matches for this reference
                        % orientation                   
                        best_match_orientation.descriptorind_cur = -1;
                        best_match_orientation.orientationind_cur = -1;
                        best_match_orientation.dist = inf;
                        second_match_orientation.dist = inf;
                        for k = 0:length(current.descriptors)-1
                            for l = 0:length(current.descriptors(k+1).orientation)-1                                
                                descriptor_cur = current.descriptors(k+1).orientation(l+1).sift_vec;                       

                                dist = sum((descriptor_ref-descriptor_cur).^2);
                                if(dist < best_match_orientation.dist)
                                    % Set second_match to previous best_match
                                    second_match_orientation.dist = best_match_orientation.dist;

                                    % set best_match to current best match
                                    best_match_orientation.descriptorind_cur = k;
                                    best_match_orientation.orientationind_cur = l;
                                    best_match_orientation.dist = dist;
                                end      
                            end
                        end
                        
                        % See if the best current descriptor for this
                        % reference orientation is less than the last best
                        if(best_match_orientation.dist < best_match_overall.dist)
                            % This orientation has a better match then the
                            % previous one; replace the previous one
                            best_match_overall.orientationind_ref = j;  
                            best_match_overall.descriptorind_cur = best_match_orientation.descriptorind_cur;    
                            best_match_overall.orientationind_cur = best_match_orientation.orientationind_cur;  
                            best_match_overall.dist = best_match_orientation.dist;  
                            second_match_overall.dist = second_match_orientation.dist;
                        end  
                    end

                    % Compare best and second match distances. If they are
                    % very close then toss the match. If distances are
                    % close then their ratios will be close to 1.                    
                    if(best_match_overall.dist/second_match_overall.dist < 0.8)
                        % This is a good point, but check list and make
                        % sure current location is unique, and if it's not
                        % then pick the best one;
                        x_ref = reference.featurepoints(i+1).x;
                        y_ref = reference.featurepoints(i+1).y;
                        x_cur = current.featurepoints(best_match_overall.descriptorind_cur+1).x;
                        y_cur = current.featurepoints(best_match_overall.descriptorind_cur+1).y;
                        
                        addpoint = true;
                        for j = 0:length(matchpoints)-1
                            if(matchpoints(j+1).x_cur == x_cur && matchpoints(j+1).y_cur == y_cur)
                                % Match found, pick the current point with
                                % the lowest distance
                                if(best_match_overall.dist < matchpoints(j+1).dist)
                                    % This new point is better, delete the
                                    % old point
                                    matchpoints(j+1) = [];   
                                    break;
                                else
                                    % Old point is better, do not add new
                                    % point
                                    addpoint = false;
                                    break;
                                end
                            end
                        end
                        
                        if(addpoint)
                            matchpoints_template.dist = best_match_overall.dist;
                            matchpoints_template.x_ref = x_ref;
                            matchpoints_template.y_ref = y_ref;
                            matchpoints_template.descriptorind_ref = best_match_overall.descriptorind_ref;
                            matchpoints_template.orientationind_ref = best_match_overall.orientationind_ref;                            
                            matchpoints_template.x_cur = x_cur;
                            matchpoints_template.y_cur = y_cur;
                            matchpoints_template.descriptorind_cur = best_match_overall.descriptorind_cur;
                            matchpoints_template.orientationind_cur = best_match_overall.orientationind_cur;
                            
                            % Store
                            matchpoints = vertcat(matchpoints,matchpoints_template);
                        end                        
                    end       
                end
            else
                error('Inputs must be of class "imagematch." The reference image must have at least one descriptor while the current image must have at least two.');
            end
        end         
        
        function mat_pose = get_pose(reference,current,matchpoints)
            if(isa(reference,'imagematch') && isa(current,'imagematch') && isa(matchpoints,'struct') && ~isempty(matchpoints))
                % Cycle over match points and form hough transform
                num_houghbins = 60; % MUST BE EVEN
                hough_val = zeros(num_houghbins,num_houghbins);
                hough_ind = cell(num_houghbins,num_houghbins);
                range_hough_x = size(current.gs,2)-1;
                range_hough_y = size(current.gs,1)-1;
                for i = 0:length(matchpoints)-1
                    % Get composition transformation matrix, apply
                    % transform to coordinates, then see the offset in
                    % coordinate values
                    mat_trans_ref = reference.featurepoints(matchpoints(i+1).descriptorind_ref+1).orientation(matchpoints(i+1).orientationind_ref+1).mat_trans;
                    scale_ref = reference.featurepoints(matchpoints(i+1).descriptorind_ref+1).scale;
                    x_ref = reference.featurepoints(matchpoints(i+1).descriptorind_ref+1).x;
                    y_ref = reference.featurepoints(matchpoints(i+1).descriptorind_ref+1).y;
                    x_ref_centroid = (size(reference.gs,2)-1)/2;
                    y_ref_centroid = (size(reference.gs,1)-1)/2;
                    
                    mat_trans_cur = current.featurepoints(matchpoints(i+1).descriptorind_cur+1).orientation(matchpoints(i+1).orientationind_cur+1).mat_trans;
                    scale_cur = current.featurepoints(matchpoints(i+1).descriptorind_cur+1).scale;
                    x_cur = current.featurepoints(matchpoints(i+1).descriptorind_cur+1).x;
                    y_cur = current.featurepoints(matchpoints(i+1).descriptorind_cur+1).y;
                    x_cur_centroid = (size(current.gs,2)-1)/2;
                    y_cur_centroid = (size(current.gs,1)-1)/2;
                    
                    % Get composite transformation
                    mat_trans_comp = mat_trans_cur*scale_cur/scale_ref*mat_trans_ref^-1;                    
                    coords_trans = mat_trans_comp*[x_ref-x_ref_centroid;y_ref-y_ref_centroid];
                                        
                    x_offset = coords_trans(1) - (x_cur-x_cur_centroid);
                    y_offset = coords_trans(2) - (y_cur-y_cur_centroid);    
                    
                    % Do hough for (x,y) coordinates, also store the
                    % indice of the matchpoint for the hough
                    ind_hough_x = floor((x_offset + range_hough_x/2)*(num_houghbins/range_hough_x));
                    ind_hough_y = floor((y_offset + range_hough_y/2)*(num_houghbins/range_hough_y));
                    
                    if(ind_hough_x > 1 && ind_hough_x < num_houghbins-2 && ind_hough_y > 1 && ind_hough_y < num_houghbins-2)
                        % If the index is on the boundary of the hough,
                        % then just discard the value as its probably wrong
                        
                        % Store with a +- 2 bin area
                        for j = -2:2
                            for k = -2:2
                                hough_val(ind_hough_y+k+1,ind_hough_x+j+1) = hough_val(ind_hough_y+k+1,ind_hough_x+j+1) + 1;
                                hough_ind{ind_hough_y+k+1,ind_hough_x+j+1} = vertcat(hough_ind{ind_hough_y+k+1,ind_hough_x+j+1},i);
                            end
                        end  
                    end                      
                end
                
                % Get indices of max hough value
                [ind_row ind_col] = find(hough_val == max(hough_val(:)));

                % Only use the first one
                ind_houghmatch = hough_ind{ind_row(1),ind_col(1)};
                
                if(length(ind_houghmatch)  >= 3)                
                    % Use these indices for least squares fit for pose
                    % estimation                
                    mat_ls = zeros(2*length(ind_houghmatch),6);
                    vec_ls = zeros(2*length(ind_houghmatch),1);
                    for i = 0:length(ind_houghmatch)-1
                        % Get composition transformation matrix, apply
                        % transform to coordinates, then see the offset in
                        % coordinate values
                        x_ref = reference.featurepoints(matchpoints(ind_houghmatch(i+1)+1).descriptorind_ref+1).x;
                        y_ref = reference.featurepoints(matchpoints(ind_houghmatch(i+1)+1).descriptorind_ref+1).y;
                        x_ref_centroid = (size(reference.gs,2)-1)/2;
                        y_ref_centroid = (size(reference.gs,1)-1)/2;

                        x_cur = current.featurepoints(matchpoints(ind_houghmatch(i+1)+1).descriptorind_cur+1).x;
                        y_cur = current.featurepoints(matchpoints(ind_houghmatch(i+1)+1).descriptorind_cur+1).y;
                        x_cur_centroid = (size(current.gs,2)-1)/2;
                        y_cur_centroid = (size(current.gs,1)-1)/2;

                        x_ref_normalized = x_ref-x_ref_centroid;
                        y_ref_normalized = y_ref-y_ref_centroid;
                        x_cur_normalized = x_cur-x_cur_centroid;
                        y_cur_normalized = y_cur-y_cur_centroid;

                        % Fill mat_ls and vec_ls
                        mat_ls(2*i+1,:) = [x_ref_normalized y_ref_normalized 0 0 1 0];
                        mat_ls(2*i+2,:) = [0 0 x_ref_normalized y_ref_normalized 0 1];

                        vec_ls(2*i+1) = x_cur_normalized;
                        vec_ls(2*i+2) = y_cur_normalized;                    
                    end

                    % Get LS solution
                    mat_trans_ls = ((mat_ls'*mat_ls)^-1)*mat_ls'*vec_ls;

                    % Reform into matrix                
                    mat_pose_prelim = zeros(3);
                    mat_pose_prelim(1,1) = mat_trans_ls(1);
                    mat_pose_prelim(2,1) = mat_trans_ls(3);
                    mat_pose_prelim(1,2) = mat_trans_ls(2);
                    mat_pose_prelim(2,2) = mat_trans_ls(4);
                    mat_pose_prelim(1,3) = mat_trans_ls(5);
                    mat_pose_prelim(2,3) = mat_trans_ls(6);

                    % Possibly perform RANSAC or some variant at this point
                    % to exclude possible bad datapoints   
                    dist_cutoff = sqrt(sum((size(current.gs)/num_houghbins).^2));                    
                    % Find points greater than cutoff
                    badpointexist = true;
                    ind_good = ind_houghmatch;
                    while(badpointexist)
                        ind_bad = [];
                        for i = 0:length(ind_good)-1
                            % Get composition transformation matrix, apply
                            % transform to coordinates, then see the offset in
                            % coordinate values
                            x_ref = reference.featurepoints(matchpoints(ind_good(i+1)+1).descriptorind_ref+1).x;
                            y_ref = reference.featurepoints(matchpoints(ind_good(i+1)+1).descriptorind_ref+1).y;
                            x_ref_centroid = (size(reference.gs,2)-1)/2;
                            y_ref_centroid = (size(reference.gs,1)-1)/2;

                            x_cur = current.featurepoints(matchpoints(ind_good(i+1)+1).descriptorind_cur+1).x;
                            y_cur = current.featurepoints(matchpoints(ind_good(i+1)+1).descriptorind_cur+1).y;
                            x_cur_centroid = (size(current.gs,2)-1)/2;
                            y_cur_centroid = (size(current.gs,1)-1)/2;

                            x_ref_normalized = x_ref-x_ref_centroid;
                            y_ref_normalized = y_ref-y_ref_centroid;
                            x_cur_normalized = x_cur-x_cur_centroid;
                            y_cur_normalized = y_cur-y_cur_centroid;

                            % Get transformed coordinates                
                            coords_trans = mat_pose_prelim*[x_ref-x_ref_centroid;y_ref-y_ref_centroid;1];

                            x_offset = coords_trans(1) - (x_cur-x_cur_centroid);                        
                            y_offset = coords_trans(2) - (y_cur-y_cur_centroid);

                            dist_offset = sqrt(x_offset^2+y_offset^2);

                            if(dist_offset > dist_cutoff)
                               ind_bad = vertcat(ind_bad,i);
                            end                        
                        end    
                                                
                        % Delete bad points from ind_good
                        if(length(ind_bad) > 0)
                            ind_good(ind_bad+1) = [];
                            
                            % See if at least 3 points exist
                            if(length(ind_good) >= 3)
                                % Recompute transformation matrix
                                mat_ls = zeros(2*length(ind_good),6);
                                vec_ls = zeros(2*length(ind_good),1);
                                for i = 0:length(ind_good)-1
                                    % Get composition transformation matrix, apply
                                    % transform to coordinates, then see the offset in
                                    % coordinate values
                                    x_ref = reference.featurepoints(matchpoints(ind_good(i+1)+1).descriptorind_ref+1).x;
                                    y_ref = reference.featurepoints(matchpoints(ind_good(i+1)+1).descriptorind_ref+1).y;
                                    x_ref_centroid = (size(reference.gs,2)-1)/2;
                                    y_ref_centroid = (size(reference.gs,1)-1)/2;

                                    x_cur = current.featurepoints(matchpoints(ind_good(i+1)+1).descriptorind_cur+1).x;
                                    y_cur = current.featurepoints(matchpoints(ind_good(i+1)+1).descriptorind_cur+1).y;
                                    x_cur_centroid = (size(current.gs,2)-1)/2;
                                    y_cur_centroid = (size(current.gs,1)-1)/2;

                                    x_ref_normalized = x_ref-x_ref_centroid;
                                    y_ref_normalized = y_ref-y_ref_centroid;
                                    x_cur_normalized = x_cur-x_cur_centroid;
                                    y_cur_normalized = y_cur-y_cur_centroid;

                                    % Fill mat_ls and vec_ls
                                    mat_ls(2*i+1,:) = [x_ref_normalized y_ref_normalized 0 0 1 0];
                                    mat_ls(2*i+2,:) = [0 0 x_ref_normalized y_ref_normalized 0 1];

                                    vec_ls(2*i+1) = x_cur_normalized;
                                    vec_ls(2*i+2) = y_cur_normalized;                    
                                end

                                % Get LS solution
                                mat_trans_ls = ((mat_ls'*mat_ls)^-1)*mat_ls'*vec_ls;

                                % Reform into matrix                
                                mat_pose_prelim = zeros(3);
                                mat_pose_prelim(1,1) = mat_trans_ls(1);
                                mat_pose_prelim(2,1) = mat_trans_ls(3);
                                mat_pose_prelim(1,2) = mat_trans_ls(2);
                                mat_pose_prelim(2,2) = mat_trans_ls(4);
                                mat_pose_prelim(1,3) = mat_trans_ls(5);
                                mat_pose_prelim(2,3) = mat_trans_ls(6);
                            else
                                % Exit loop
                                break;
                            end      
                        else
                            % No bad points exist - exit loop
                            badpointexist = false;
                        end
                    end
                    
                    if(length(ind_good) >= 3)
                        % Assign output
                        mat_pose = mat_pose_prelim;
                    else
                        error('At least three matchpoints with agreeable pose during post hough filtering were not identified.');
                    end 
                else
                    error('At least three matchpoints with consistent pose from the hough were not identified.');
                end
            else
                error('Inputs must be of class "imagematch" and a structure for the matchpoints. Matchpoints cant be empty either.');
            end
        end
        
        function plot_matchpoints(reference,current,matchpoints)
            if(isa(reference,'imagematch') && isa(current,'imagematch') && isa(matchpoints,'struct') && ~isempty(reference.gs) && ~isempty(current.gs))
                offset_x_current = size(reference.gs,2);                
                
                combinedimage = zeros(max(size(reference.gs,1),size(current.gs,1)),size(reference.gs,2)+size(current.gs,2));
                combinedimage(1:size(reference.gs,1),1:size(reference.gs,2)) = reference.gs;
                combinedimage(1:size(current.gs,1),(1:size(current.gs,2))+offset_x_current) = current.gs;
                
                % Plot
                figure, imshow(combinedimage,[]);
                for i = 0:length(matchpoints)-1
                    line([matchpoints(i+1).x_ref+1 matchpoints(i+1).x_cur+offset_x_current+1],[matchpoints(i+1).y_ref+1 matchpoints(i+1).y_cur+1]);                    
                end
            else
                error('Inputs must be of class "imagematch" and a structure for the matchpoints. The reference and current images must be loaded as well.');
            end
        end
        
        function plot_pose(reference,current,mat_pose)
            if(isa(reference,'imagematch') && isa(current,'imagematch') && isnumeric(mat_pose) && ~isempty(reference.gs) && ~isempty(current.gs) && isequal(size(rand(3)),[3 3]))
                offset_x_current = size(reference.gs,2);                
                
                combinedimage = zeros(max(size(reference.gs,1),size(current.gs,1)),size(reference.gs,2)+size(current.gs,2));
                combinedimage(1:size(reference.gs,1),1:size(reference.gs,2)) = reference.gs;
                combinedimage(1:size(current.gs,1),(1:size(current.gs,2))+offset_x_current) = current.gs;
                
                x_ref_centroid = (size(reference.gs,2)-1)/2;
                y_ref_centroid = (size(reference.gs,1)-1)/2;
                x_cur_centroid = (size(current.gs,2)-1)/2;
                y_cur_centroid = (size(current.gs,1)-1)/2;
                
                % Perform transformation on image boundary
                line_x = [-x_ref_centroid x_ref_centroid, ...
                          x_ref_centroid x_ref_centroid, ...
                          x_ref_centroid -x_ref_centroid, ...
                          -x_ref_centroid -x_ref_centroid];   
                
                line_y = [-y_ref_centroid -y_ref_centroid, ...
                          -y_ref_centroid y_ref_centroid, ...
                          y_ref_centroid y_ref_centroid, ...
                          y_ref_centroid -y_ref_centroid];
                
                line_trans = mat_pose*vertcat(line_x,line_y,ones(1,length(line_x)));
                
                % Now offset according to centroid
                offset_x = offset_x_current+x_cur_centroid;
                offset_y = y_cur_centroid; 
                
                line_cur_x = line_trans(1,:) + offset_x;
                line_cur_y = line_trans(2,:) + offset_y;
                
                % Plot
                line_ref_x = [0 2*x_ref_centroid, ...
                              2*x_ref_centroid 2*x_ref_centroid, ...
                              2*x_ref_centroid 0, ...
                              0 0];   
                
                line_ref_y = [0 0, ...
                              0 2*y_ref_centroid, ...
                              2*y_ref_centroid 2*y_ref_centroid, ...
                              2*y_ref_centroid 0];
                figure, imshow(combinedimage,[]);
                for i = 0:3
                    % plot original box
                    line([line_ref_x(2*i+1:2*i+2)+1],[line_ref_y(2*i+1:2*i+2)+1],'LineWidth',2)
                    
                    % plot connector line
                    line([line_ref_x(2*i+1)+1 line_cur_x(2*i+1)+1],[line_ref_y(2*i+1)+1 line_cur_y(2*i+1)+1])
                    line([line_ref_x(2*i+2)+1 line_cur_x(2*i+2)+1],[line_ref_y(2*i+2)+1 line_cur_y(2*i+2)+1])
                    
                    % plot transformed box
                    line([line_cur_x(2*i+1:2*i+2)+1],[line_cur_y(2*i+1:2*i+2)+1],'LineWidth',2);                    
                end
            else
                error('Inputs must be of class "imagematch" and a 3x3 matrix for the mat_pose. The reference and current images must be loaded as well.');
            end
        end
    end
end