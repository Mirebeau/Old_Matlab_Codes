% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% This file illustrates the Edge Enhancing Diffusion (cEED) and 
% Coherence Enhancing Diffusion (cCED) filters of J. Weickert,
% discretized using the monotony preserving explicit scheme published in :

% J. Fehrenbach, J.-M. Mirebeau, Sparse non-negative stencils for anisotropic diffusion,
% J. Math. Imag. Vis., vol. 49(1) (2014), pp. 123-147

% Main options fields :  
% - Weickert_lambda (edge detection threshold)
% - final_time (PDE evolution time)

% Secondary options fields : 
% - Weickert_choice ('cEED','cCED','EED','CED'. Choice of PDE) 
% - Weickert_alpha (diffusion tensors condition number is <=1/alpha)
% - Weickert_m (exponent in Weickert's tensors construction)

% - noise_filter, feature_filter (for structure tensor construction)
% - rescale for unit maximum trace (rescale structure tensors, true by default)

% - max_diff_iter (max number of time steps, and diffusion tensor updates)
% - max_inner_iter (number of inner time steps, between diffusion tensor updates)

% - verbose (true or false)

% Remark on performance: On 'large' cases, such as the MRI below, computation time 
% is dominated by the sparse matrix assembly "spmat(col,row,coef,n,n)". 
% In case of need, consider the following optimized C++ implementation designed for
% the Insight Toolkit (ITK) 

% J. Fehrenbach, J.-M. Mirebeau, L. Risser, S. Tobji,
% Anisotropic Diffusion in ITK, Insight Journal, 2015
% http://www.insight-journal.org/browse/publication/953

addpath('../GetOptions');
addpath('../AD-LBR');
addpath('TensorConstruction');
addpath('../Eig3Folder');

disp('----------------- Demo : MRI -----------------')
clear options;
img=double(hdf5read('ImageData/mrbrain_noisy_01.hdf5','/ITKImage/0/VoxelData'))/255;
%options.Weickert_choice = 'cEED'; %Edge enhancing diffusion (default)
options.Weickert_lambda = 0.003; %Edge detection threshold.
options.final_time=8; %PDE evolution time.
options.max_inner_iter=3;

smoothed=NonLinearDiffusion_3D(img,options);
imshow([img(:,:,50),smoothed(:,:,50)]);
pause();
imshow([squeeze(img(:,120,:)),squeeze(smoothed(:,120,:))]);
pause();
imshow([squeeze(img(100,:,:)),squeeze(smoothed(100,:,:))]);
pause();

disp('---------------- Demo : Cos3D ---------------')
clear options;
img=double(hdf5read('ImageData/Cos3D_Noisy.hdf5','/ITKImage/0/VoxelData'))/255;
options.Weickert_choice = 'cCED'; 
options.Weickert_lambda = 0.02; %Edge detection threshold.
options.final_time=10; %PDE evolution time.

options.noise_filter = fspecial('gaussian',[10,1],4);
options.feature_filter = fspecial('gaussian',[16,1],5);

smoothed=NonLinearDiffusion_3D(img,options);
imshow([img(:,:,90),smoothed(:,:,90)]);
pause();