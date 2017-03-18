function [ output ] = Togray( temp )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
output = .299*temp(:,:,1) + .587*temp(:,:,2) + .114*temp(:,:,3);

end

