  function [ edgeR,edgeL ] = match( R,L,smoothfactor )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% R_value=double(R);
% L_value=double(L);

% for i=1+smoothfactor:1:r-smoothfactor
%     for k=1+smoothfactor:1:c-smoothfactor
%         R_value(i,k)=0.5*(mean(double(R(i,[k-smoothfactor:k+smoothfactor])))+mean(double(R(i-smoothfactor:i+smoothfactor,k))));
%     end
% end
% 
% [r,c]=size(L);
% for i=1+smoothfactor:1:r-smoothfactor
%     for k=1+smoothfactor:1:c-smoothfactor
%         L_value(i,k)=0.5*(mean(double(L(i,[k-smoothfactor:k+smoothfactor])))+mean(double(L(i-smoothfactor:i+smoothfactor,k))));
%     end
% end
[r,c]=size(R);
highfilter=[1 0 -1 0 1; 0 1 -1 1 0; -1 -1 0 -1 -1;0 1 -1 1 0;1 0 -1 0 1 ];
Diagfilter=[4 -1 -1 -1 -1; -1 4 -1 -1 -1; -1 -1 4 -1 -1; -1 -1 -1 4 -1; -1 -1 -1 -1 4];
DiagLfilter=[-1 -1 -1 -1 4; -1 -1 -1 4 -1; -1 -1 4 -1 -1; -1 4 -1 -1 -1;4 -1 -1 -1 -1];
R(r:r+4,c:c+4)=0;
for i=1:r
    for k=1:c
        R(i,k)=sum(sum(R(i:i+4,k:k+4).*Diagfilter));
        
    end
end



R_value=double(R);
L_value=double(L);





[r,c]=size(R);

magnitudeR=zeros(r,c);
orientationR=zeros(r,c);
edgeR=zeros(r,c);
for i=2:1:r-1
    for k=2:1:c-1
        if (R_value(i,k)~=0)&&(R_value(i,k+1)~=0)&&(R_value(i+1,k)~=0)&&(R_value(i,k+1)~=0)&&(R_value(i,k-1)~=0)
            edgeR(i,k)=abs(R_value(i-1,k)-R_value(i,k))+abs(R_value(i,k)-R_value(i+1,k))+abs(R_value(i,k-1)-R_value(i,k))+abs(R_value(i,k)-R_value(i,k+1));
            magnitudeR(i,k)=100*((R_value(i,k)-R_value(i,k+1))^2+(R_value(i,k)-R_value(i+1,k))^2)^0.5;
            orientationR(i,k)=atan2(R_value(i,k)-R_value(i,k+1),R_value(i,k)-R_value(i+1,k));
        end 
    end
end

[r,c]=size(L);

magnitudeL=zeros(r,c);
orientationL=zeros(r,c);
for i=2:1:r-1
    for k=2:1:c-1
        if (L_value(i,k)~=0)&&(L_value(i,k+1)~=0 )&&( L_value(i+1,k)~=0)&&(L_value(i,k+1)~=0)&&(L_value(i,k-1)~=0)
            edgeL(i,k)=abs(L_value(i-1,k)-L_value(i,k))+abs(L_value(i,k)-L_value(i+1,k))+abs(L_value(i,k-1)-L_value(i,k))+abs(L_value(i,k)-L_value(i,k+1));
            magnitudeL(i,k)=100*((L_value(i,k)-L_value(i,k+1))^2+(L_value(i,k)-L_value(i+1,k))^2)^0.5;
            orientationL(i,k)=atan2(L_value(i,k)-L_value(i,k+1),L_value(i,k)-L_value(i+1,k));
        end 
    end
end



end

