function [ prop ] = PropMesh( lattice, numUnitCells, numPixels, zSpacing, E, angleX, angleY  )
%PROPMESH Outputs the propagator in momentum space. 
%The units for lattice and zSpacing are as before - angstroms, the accelerating 
%voltage E is in Volts, and the angles along x and y that the sample is displaced 
%are in mrad, though angles should be less than 1 degree for the approximation to work.

    lambda = 12.264306/(sqrt(E)*sqrt(1+0.97846707*10^(-6)*E));
    pixelSize = 1./lattice(1:2)./numUnitCells;
    
    %Assuming a square unit cell, the sapce in kSpace is the same in both
    %directions
    kSpace(1:(numPixels/2)) = ((1:(numPixels/2))-1);
    kSpace((numPixels/2):numPixels) = (((numPixels/2):numPixels) - numPixels-1);
    [Kx,Ky] = meshgrid(kSpace);
    Kx = Kx*pixelSize(1);
    Ky = Ky*pixelSize(2);
    K2 = (Kx.^2+Ky.^2);
    
    prop = exp(-1j*pi*lambda*K2*zSpacing + 2*pi*1j*zSpacing...
        *(Kx*tan(angleX*0.001)+Ky*tan(angleY*0.001)));
    
end

