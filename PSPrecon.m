function [ output] = PSPrecon(ref1,ref2,ref3, org1,org2,org3 )
% UNTITLED2 Summary of this function goes here

%example as follow:
% [ref1 ref2 ref3 org1 org2 org3 shape] = PSPGenerateRefAndConeFringe(512, 1, 16, 'cone', 0);
% h=PSPrecon(ref1,ref2,ref3,org1,org2,org3);
% mesh(h);
% ref1 refers to object pattern with phase shift -2*pi/3,while ref2 no phase
% shift, ref3 +2*pi/3 phase shift
% org means the reference pattern with the similar phase shift with series
% of ref1-ref3.
%   Detailed explanation goes here
[r,c]=size(org1);
intensityofref=zeros(r,c);%I of the reference
intensityoforg=zeros(r,c);%I of the orginal object pattern
angle=zeros(r,c);%reference ringe pattern angle
angle1=zeros(r,c);%real object pattern
offset=4.0;%phase unwrapping when 2 straight numbers difference larger than offset, then add 2*pi to the letter one



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
    angle(i,:)=unwrap(angle(i,:));
end
% for k=1:1:c
%     angle(:,k)=unwrap(angle(:,k));
% end
        

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
angle1(i,:)=unwrap(angle1(i,:));
end

output=angle-angle1;% in order to eliminate the factor of 02.*pi.*fo/fs.*X inside the COS() function, after that get the pure value of the height information


end
                            
                                    
                    
        
        
        
        
        
        
        
        
        
        
        
        


