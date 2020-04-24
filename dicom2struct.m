function [ARRAY,STRUCT] = dicom2struct(path,filename)

%==========================================================================
%  Subroutine to read dicom images to structure
%==========================================================================
%
% Read original images and create structure with Image, Location and header
% 
% INput:        path
% OUTput:       saves *.mat in path
%               outputs struct to workspace
%--------------------------------------------------------------------------
% written by Vadim Malis
% 04/16 at UCSD RIL
%==========================================================================
warning('off', 'images:dicominfo:fileVRDoesNotMatchDictionary')

cd(path)
dicomlist=dir('*.dcm');
    
% get header info for current set from 1st image header
dicom_header = dicominfo(dicomlist(1).name);
numim = dicom_header.ImagesInAcquisition;
r = dicom_header.Rows;
c = dicom_header.Columns;

for i=1:numim
	STRUCT(i).Image=dicomread(dicomlist(i).name);
	STRUCT(i).header=dicominfo(dicomlist(i).name);
    STRUCT(i).echoN = STRUCT(i).header.EchoNumber;
    STRUCT(i).TE = STRUCT(i).header.EchoTime;
	STRUCT(i).location=STRUCT(i).header.SliceLocation;
end

echos=max([STRUCT.echoN]);

ARRAY=reshape([STRUCT.Image],[r,c,numim/echos,echos]);
save(filename,'STRUCT')