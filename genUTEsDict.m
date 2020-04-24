function dictionary = genUTEsDict(TEs,Freq,T2s)

%==========================================================================
% Subroutine to creat dictionary for UTE fitting
%==========================================================================
%
% INput:  TEs           column vector of UTEs TEs in ms
%         Freq          column vector for range of frequencies
%         T2s           column vector for range of T2s 
%
% OUTput: dictionary    matrix size = [sizeTEs x sizeFreq, sizeT2s, 2 ]
%
%--------------------------------------------------------------------------
% Original paper:
%
%   C. A. Araujo, E., Azzabou, N., Vignaud, A., Guillot, G. and Carlier, P. 
%   Quantitative ultrashort TE imaging of the short‐T2 components in 
%   skeletal muscle using an extended echo‐subtraction method. 
%   Magn. Reson. Med., 78: 997-1008. doi:10.1002/mrm.26489
%
%--------------------------------------------------------------------------
% written by Vadim Malis
% 04/20 at UCSD RIL
%==========================================================================

% function handles
f1 = @UTE_expTerm;
f2 = @fat_spectrum;


% lipids 
lipids_freq = [-404; -307; 91; -452; -221];      %Hz
lipids_peak = [0.66; 0.16; 0.09; 0.07; 0.01];    %unitless 

f_lipids = f2(lipids_peak,lipids_freq,TEs); % this only has n=number of TEs
% replicate by number of freq and T2s to allow pairwise multiplication
F_lipids = squeeze(repmat(f_lipids,[1,size(Freq),size(T2s)]));
F_lipids = permute(F_lipids,[2,1,3]);

% generate all possible combinations
[TE, FQ, T2] = meshgrid(TEs,Freq,T2s);

Sw=f1(TE,FQ,T2);
Sf=F_lipids.*Sw;

dictionary = cat(4,Sw,Sf);
dictionary = permute(dictionary,[1,3,2,4]);

end


%% exp factor
function x = UTE_expTerm(TE,freq,T2)
    x = exp((1i*2*pi*freq-1/T2).*TE);
end

%% fat
function f = fat_spectrum(rho,freq,t)
   temp = rho.*exp(1i*2*pi*(freq.*t));
   f=sum(temp,1);
   f=f';
end
