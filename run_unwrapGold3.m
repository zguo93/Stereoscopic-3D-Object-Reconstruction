%
%
% [unwphx_qual, brcx, resx] = run_unwrapGold3(thephase, bmask)
% 
% Enable the dipole removal procedure to yes.
% 
% thephase: input phase image
% unwphx_qual: unwrapped phase image
% brcx: branch cut maps.
% resx: residues.
% 
% Script for running the Goldstein's method for phase unwrapping.
% The native program is compiled from the code supplied from the book:
% 
% D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
% Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.
%
% * Use tempname to decrease probabilty of error when running parallel.
% 
% Hsung Tai-Chiu, Richard.
% April 29, 2009.
% Nov 12, 2009.
% 
function [unwphx_qual] = run_unwrapGold3(thephase, bmask)

temp_prefix = tempname;

% % % % % % % % % % % % % % % % % %%%%%%%%%%%
% Write file for input
% % % % % % % % % % % % % % % % % %%%%%%%%%%%
szphase = size(thephase);
fname = [temp_prefix 'im_noisy.phase'];
outfname = [temp_prefix 'im_noisy.uw'];
mname = [temp_prefix 'im_noisy.bmask'];
fid = fopen(fname, 'w');
fwrite(fid,thephase+pi,'float');
fclose(fid);
fid = fopen(mname, 'wb');
fwrite(fid,bmask,'uint8');
fclose(fid);

% % % % % % % % % % % % % % % % % %%%%%%%%%%%
% Compose a command line
% % % % % % % % % % % % % % % % % %%%%%%%%%%%
% unwrap_gold
commandline = [ 'GoldStein -input ', fname, ' -format float -output ', outfname,...
                ' -xsize ', num2str(szphase(1)), ...
                ' -ysize ', num2str(szphase(2)), ...
                ' -mask ', mname, ...
                ' -dipole yes', ...
                ' -debug yes' ]
tic;            
system(commandline);
disp(['Goldstein Time Command Line: ', num2str(toc)]);

% % % % % % % % % % % % % % % % % %%%%%%%%%%%
% Load the results
% % % % % % % % % % % % % % % % % %%%%%%%%%%%
fid = fopen(outfname,'r');
unwphx_qual  = fread(fid,szphase,'float32');
fclose(fid);
% brcname = [outfname, '.brc'];
% fid = fopen(brcname,'r');
% brcx  = fread(fid,szphase,'uint8');
% fclose(fid);
% resname = [outfname, '.res'];
% fid = fopen(resname,'r');
% resx  = fread(fid,szphase,'uint8');
% fclose(fid);
% 
% 
% % % % % % % % % % % % % % % % % % %%%%%%%%%%%
% % clean!
% % % % % % % % % % % % % % % % % % %%%%%%%%%%%
% commandline = ['del ', outfname,' ', fname,' ', mname,' ', brcname, ' ', resname]
commandline = ['del ', outfname,' ', fname,' ', mname]
 system(commandline);


return