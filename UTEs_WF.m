%------------------------------------------
% UTEs water/fat fraction
%
%   Input:  5x1 column vector - utes signals S1...S5
%           precalculated dict          
%    
%   Output: arrays, 
%           w    - water
%           f    - fat
%           Freq - frequency
%           T2l  - T2 Long
%           
%------------------------------------------

%% read images and TEs
disp("Step 1 of 4:   reading images...")
[UTEs,UTEs_Struct] = dicom2struct(cd,'UTEs');
info_temp=[UTEs_Struct.header];
TEs=unique([info_temp.EchoTime])*1e-3;
N = size(UTEs,3); %number of slices
clear info_temp
disp("reading complete!")

%% generating mask for bckg pixels
disp("Step 2 of 4:   creating masks...")
UTE_E1=UTEs(:,:,:,1); %only 1st echo (highest signal intensity)
Mask = UTEsMask(UTE_E1,100);
disp("masks created!")


%% optional to keep only echose 2 to last
UTEs(:,:,:,1)=[];
TEs(1)=[];


%% generate  dictionary
disp("Step 3 of 4:   generating dictionary")    
    
    %-------dictionary range----------
    %frequencies range
    max_Freq    = 60;   %Hz
    min_Freq    = -60;  %Hz
    delta_Freq  = 1;    %Hz
    %T2s range
    max_T2  = 60e-3;    %s
    min_T2  = 10e-3;    %s
    delta_T2= 0.25e-3;  %s
    
    FREQ=min_Freq:delta_Freq:max_Freq;
    T2s=min_T2:delta_T2:max_T2;
    %----------------------------------
    
dict=genUTEsDict(TEs,FREQ,T2s);
disp("dictionary generated!")

%% finding optimum values for freq, T2, W and F
%---------------------------------------------
% Result -> Structure with  W, F, freq, T2 all same size 
disp("Step 4 of 4:   minimization")
WaitBar = waitbar(0,'Minimization...');

% number of voxels to process:
total=nnz(Mask);
processed=0;

% allocate arrays
W=zeros(size(UTE_E1));
F = W;
freq = W;
T2 = W;

% per voxel array allocation
num_freq = size(min_Freq:delta_Freq:max_Freq,2);
num_t2s = size(min_T2:delta_T2:max_T2,2);

w=zeros([num_freq, num_t2s]);
f=w;
% minimization matrix
Z = w;

%%loop through each pixel
for i=1:N
%per slice    
    for x=1:size(UTEs,1)
    %per column
 
        for y=1:size(UTEs,2)
        %per row
            
            %process only non-zero pixels
            if Mask(y,x,i)>0

                parfor frq=1:num_freq
                % parallel  per frequency offset  
                    for T=1:num_t2s
                    %per T2
                        
                        % Step I    get values for w and f
                        %           each is (frq x T) size matrix
                        [w(frq,T),f(frq,T)] = wf_decomp(UTEs(y,x,i,:), dict(frq,T,:,:));
               
                        % Step II   get sum of norms Equation (11) 
                        Z(frq,T) = UTEsDictSumNorm(double(squeeze(UTEs(y,x,i,:))), dict(frq,T,:,:), w(frq,T),f(frq,T));
                        
                    end %per frequency offset  
                end % per T2
                
                
                % Step III Get indices for min; extract freq, T2, W and F
                [minFrq,minT]=getMatMinIndex(Z);
                W(y,x,i)    =   w(minFrq,minT);
                F(y,x,i)    =   f(minFrq,minT);
                freq(y,x,i) =   FREQ(minFrq);
                T2(y,x,i)   =   T2s(minT);
                
                processed=processed+1;
                waitbar(processed/total,WaitBar,sprintf('Minimization: %.5f %% complete...',round(processed/total,5)));
                
            else
                W(y,x,i)    =   NaN;
                F(y,x,i)    =   NaN;
                freq(y,x,i) =   NaN;
                T2(y,x,i)   =   NaN;
                
            end
        
        end %per row
    end %per column
end %per slice    

clearvars -except W F freq T2 UTEs


%% 
function [w,f] = wf_decomp(UTEs, dict)
    % must be double and remove singleton dimensions
    dict=squeeze(dict);
    UTEs=double(squeeze(UTEs));
    x = dict\UTEs;
    w=x(1);
    f=x(2);
end

function x = UTEsDictSumNorm(UTEs, dict, w, f)
    % per echo difference
    %warning('off', 'MATLAB:rankDeficientMatrix')
    dict=squeeze(dict);
    d = (UTEs.^2 - abs(dict*[w;f]).^2).^(0.5);  %L2 norm
    x = abs(sum(d,1));
end
          