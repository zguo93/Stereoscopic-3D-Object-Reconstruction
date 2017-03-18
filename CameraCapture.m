%
%
% function  capture(camera, fileName)
% 
% Budianto
% August 27, 2014.
% 
% 
function  CameraCapture(camera, fileName)



% % % % % % % % % % % % % % % % % %%%%%%%%%%%
% Compose a command line
% % % % % % % % % % % % % % % % % %%%%%%%%%%%
% unwrap_gold
commandline = [ 'C:\Users\fyp2\Desktop\Stereoscopic\matlab\CameraControl ' num2str(camera) ' ' fileName ]
tic;            
system(commandline);
disp(['Capture time: ', num2str(toc)]);



return