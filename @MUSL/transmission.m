function [ trans ] = transmission( V, E, zSpacing )
%TRANSMISSION Given the projected potential V across some grid in real
%space and the accelerating voltage of the beam E in volts, the output 
%is the transmission function in real space

    %relativistic electron interactivion constant
    sigma = (0.25615739*(1+1.9569341*10^-6*E))/(sqrt(E)*sqrt(1+0.97846707*10^-6*E));
    
    trans = exp(1j*sigma*V*zSpacing);
    
    %scaling factor from equation 5.21 of Kirkland(2010) lambda*lorentz
    %factor
    %lambda = 12.264306/(sqrt(E)*sqrt(1+0.97846707*10^(-6)*E));
    %scale  = lambda*(1 + E/(0.5109989461*10^6));
    %trans = exp(1j*scale*V*zSpacing);
end

