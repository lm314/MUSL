function fe = ElectronScatFactor( q2, element )
%ELECTRONSCATFACTOR gives the electron scatter factor in angstroms for the 
% elements with atomic numbers from 1 to 103 using the approximations 
% mentioned in Appendix C of Kirkland (2010) Advanced
% Computing in Electron Microscopy. The parameterization coefficients are
% stored in electronScatteringParameterization.mat and are cached in
% memory between function calls.

persistent paramMat
if(isempty(paramMat))
    %keep in memory locally between calls
    temp = struct2cell(load('electronScatteringParameterization.mat'));
    paramMat = temp{1};
end

if(element <= size(paramMat,1))
    a1 = paramMat(element,1);
    b1 = paramMat(element,2);
    a2 = paramMat(element,3);
    b2 = paramMat(element,4);
    a3 = paramMat(element,5);
    b3 = paramMat(element,6);
    c1 = paramMat(element,7);
    d1 = paramMat(element,8);
    c2 = paramMat(element,9);
    d2 = paramMat(element,10);
    c3 = paramMat(element,11);
    d3 = paramMat(element,12); 
else
    a1 = 0;
    b1 = 0;
    a2 = 0;
    b2 = 0;
    a3 = 0;
    b3 = 0;
    c1 = 0;
    d1 = 0;
    c2 = 0;
    d2 = 0;
    c3 = 0;
    d3 = 0; 
    warning('Parameterization does not exist for element number %d\n',element)
end

    fe = a1./(q2+b1) + a2./(q2+b2) + a3./(q2+b3) ...
        + c1*exp(-d1*q2) + c2*exp(-d2*q2) + c3*exp(-d3*q2);
end