function [ transformed_cur,transformed_ref ] = transform( bestmatrix,reference,current,ref,cur )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[r,c]=size(ref);

offset_x_current = size(reference.gs,2);   
combinedimage = zeros(max(size(reference.gs,1),size(current.gs,1)),size(reference.gs,2)+size(current.gs,2));
combinedimage(1:size(reference.gs,1),1:size(reference.gs,2)) = reference.gs;
combinedimage(1:size(current.gs,1),(1:size(current.gs,2))+offset_x_current) = current.gs;
for i=1:r
    for k=1:c
        Changed=[i k 1]*bestmatrix;
        transformed_ref(round(Changed(1))+800,round(Changed(2)))=ref(i,k);
        
    end
end
for i=1:r
    for k=1:c

        transformed_cur(i+500,k+500)=cur(i,k);
        
    end
end

end

