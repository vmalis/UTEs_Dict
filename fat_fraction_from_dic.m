
%UTEs structure, including first echo:

[UTEs,UTEs_Struct] = dicom2struct(pwd,'UTEs');

UTEs = double(UTEs); 

TEs = unique( [UTEs_Struct.TE] )*1e-3;

% Calculating fat dephasing factors for all TE's:

lipids_freq = [-404; -307; 91; -452; -221];      
lipids_peak = [0.66; 0.16; 0.09; 0.07; 0.01];

temp = lipids_peak.*exp(1i*2*pi*(lipids_freq.*TEs));
f=sum(temp,1);
f=f'; 
f = abs(f);

% Calculating exponential terms for first two echoes:

E_0_l = exp(-TEs(1)./T2); 

E_1_l = exp(-TEs(2)./T2);

% Extracting first two echoes:

S0 = UTEs(:,:,:,1);
S1 = UTEs(:,:,:,2); 

% Calculating Fat Fraction (FF) and macromolecular fraction  (CF):

FF = abs(F) .* ( E_0_l ./ S0 ) ; 
C_sharp = (S0./E_0_l) - (S1./E_1_l) - ( abs(F) .* (f(1) - f(2)) ) ;
CF = C_sharp .* ( E_0_l ./ S0 ) ; 























