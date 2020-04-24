function [a,b] = getMatMinIndex(c)
%==========================================================================
% Subroutine to get matrix min element indiceis
%==========================================================================
%
% INput:  c
%
% OUTput: a, b indicies
%
%--------------------------------------------------------------------------
% written by Vadim Malis
% 04/20 at UCSD RIL
%==========================================================================

as=size(c);
[~,I]=min(c(:));
r=rem(I,as(1));
a=r;
b=((I-a)/as(1))+1;
    if a==0
        a=as(1);
        b=b-1;
    else
        a=r;
        b=b;
    end
    
end