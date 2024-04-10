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
% - rescale for unit maximum trace (rescale structure tensors)

% - max_diff_iter (max number of diffusion tensor updates)
% - max_inner_iter (number of inner time steps, between diffusion tensor updates)

% - verbose (true or false)

% Remark on performance: On 'large' cases, such as the MRI below, computation time 
% is dominated by the sparse matrix assembly "spmat(col,row,coef,n,n)". 
% In case of need, consider the following optimized C++ implementation designed for
% the Insight Toolkit (ITK) 
% http://www.insight-journal.org/browse/publication/953

addpath('../GetOptions');
addpath('../AD-LBR');
addpath('TensorConstruction');

disp('---------------- Demo : Color Lena image ---------------');

clear options;
img=double(imread('ImageData/lena.png'))/255;
options.Weickert_lambda=0.003;
options.Weickert_exponent=4;
options.final_time=2;

smoothed=NonLinearDiffusion_2D(img,options);
imshow([img,smoothed]);
imwrite(smoothed,'ImageResults/smoothed_lena.png');
pause();

disp('------------ Demo : FingerPrint ---------------')

clear options;
img=double(imread('ImageData/FingerPrint.png'))/255;
options.Weickert_choice = 'cCED'; %Coherence enhancing diffusion (only touches image edges)
options.Weickert_lambda=0.01;
options.final_time=7;

smoothed=NonLinearDiffusion_2D(img,options);
imshow([img,smoothed]);
imwrite(smoothed,'ImageResults/smoothed_FingerPrint.png');
pause();


disp('----------- Demo : Grayscale Pac-Man image ------------')
clear options;
img=double(imread('ImageData/noisy_pac_man.png'))/255;
options.Weickert_lambda = 0.1;
options.final_time = 12;

% Due to the high noise level, we use a wider smoothing kernel for the structure tensor.
options.noise_filter = fspecial('gaussian',9,3);

smoothed=NonLinearDiffusion_2D(img,options);
imshow([img,smoothed]);
imwrite(smoothed,'ImageResults/smoothed_pac_man.png')
pause();

disp('----- Demo : Coherence enhancing diffusion only touches image edges ----');

options.Weickert_choice = 'cCED'; 
smoothed=NonLinearDiffusion_2D(img,options);
imshow([img,smoothed]);
imwrite(smoothed,'ImageResults/smoothed_pac_man.png')
