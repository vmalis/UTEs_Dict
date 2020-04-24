function [Mask] = UTEsMask(Images,intensity_threshold)

%==========================================================================
% Subroutine to creat full axial slice mask
%==========================================================================
%
% INput:  Images                  image volume [size = xres, yres, nslices]
%         intensity_threshold     intensity threshold
%
% OUTput: Mask [size = xres, yres, nslices]
%
%--------------------------------------------------------------------------
% written by Vadim Malis
% 04/20 at UCSD RIL
%==========================================================================
thrsh = intensity_threshold;
I = Images;

N=size(I,3); %number of slices
Mask=zeros(size(I,1),size(I,2),N);

    for n=1:N
        bw=zeros(size(I,1),size(I,2));
        bw(I(:,:,n,1)>thrsh)=1;
        bw=logical(bw);
        
        % keep only one big object and fill holes
        temp = imfill((bwareafilt(bw,1)),'holes');
        
        % blur by convolution to smooth edges
            windowSize = 20;
            kernel = ones(windowSize) / windowSize ^ 2;
            blurryImage = conv2(single(temp), kernel, 'same');
        
        % store    
        Mask(:,:,n) = blurryImage > 0.5; % Rethreshold
       
    end

end