function Vmesh = CrystalPotentialMax( unitCellSites, lattice, numPixels, numUnitCells,element, E, kMax, absorp, temperature)
%CRYSTALPOTENTIAL Calculates the projected potential (in Volts) for   
%the layer indicated by the array unitCellSites, which is an array 
%with the position of atoms in units of the unit cell length. 
%It outputs the potential in real space Vmesh by calculating the potential
%in k space for those points that fall withing kMax (angstroms^-1).
%The output is properly scaled to account for the inverse fft that gives
%the final output using the factor (numPixels)^2. The potential also
%includes random thermal motion of the atoms via the Debye Weller 
%temperature factor DW (in angstroms^2) as well as absorption.

    %The zero of the momentum space is in the upper left corner (to suit
    %FFT2 formatting. Keep in mind the rearranged coordinate system that
    %comes with the use of FFT2
    pixelSize = 1./lattice./numUnitCells;
    Vmesh = zeros(numPixels);
    phase = zeros(numPixels);
    [atomPositions,element] = MUSL.populateUnitCells(unitCellSites,element,numUnitCells,lattice);
    
    % Make sure that maximum reciprocal lattice vector kMax is less than the radius
    % of the maximum circle that can fit in simulation grid. We wish to have
    % rotational symmetry in the final result.
    if(kMax^2 > min((pixelSize(1:2)*numPixels/2).^2))
        kMax = min(abs((pixelSize(1:2)*numPixels/2)));
    end
    k2max = kMax^2;
    
    %fill coordinates in k-space, keeping in the origin is at the
    %upper left hand corner and the quadrants in the fourier transform are
    %shifted so that they are going from left to right top to bottom, 4th,
    %3rd, 1st, 2nd. We need only make one kSpace vector as it is the same
    %along both (assuming a square unit cell). Then make the 2d grid with 
    %kx,ky and finally calculate the k vector squared.
    kSpace(1:(numPixels/2)) = ((1:(numPixels/2))-1);
    kSpace((numPixels/2):numPixels) = (((numPixels/2):numPixels) - numPixels-1);
    [Kx,Ky] = meshgrid(kSpace);
    Kx = Kx*pixelSize(1);
    Ky = Ky*pixelSize(2);
    k2 = (Kx.^2+Ky.^2);
    
        
    %add extra negative sign to phase to compensate to for sign difference in fft2
    for atom = 1:size(atomPositions,1)  %sums over all atoms in particular layer
        DW = MUSL.DebyeWaller(element(atom),temperature);
        phase(k2 <= k2max) = exp(-2*pi*1j*(Kx(k2 <= k2max)*atomPositions(atom,1) ...
            +Ky(k2 <= k2max)*atomPositions(atom,2)));
        %Note that the factor of 1/4 is included since the Debye-Waller 
        %factor uses s as is variable and not q and s = q/2
        Vmesh(k2 <= k2max) = Vmesh((k2 <= k2max)) + ...
            + MUSL.ElectronScatFactor(k2(k2 <= k2max),element(atom)).*...
            exp(-DW/4*k2(k2 <= k2max)).*phase(k2 <= k2max);
        if(absorp)
            %The factor of 1/4 is included here as well for the same reason as
            %above
            Vmesh(k2 <= k2max) = Vmesh(k2 <= k2max) +...
                1j*MUSL.AbsorptivePotential( k2(k2 <= k2max)/4,element, E, temperature)...
                .*phase(k2 <= k2max);
        end
    end
    
    %Self et al (1983)
    Vmesh(k2 <= k2max) = 47.878009/prod(lattice)*Vmesh(k2 <= k2max);
    
    %Determine number of beams (ie pixels included in k2max area)
    %numBeams = sum(sum(k2 <= k2max));
    %fprintf('Maximum excitation error is %0.2f angstroms^-1\n', kMax)
    %fprintf('Number of beams: %d\n', numBeams)
        
    %renormalization to account for ifft2
    Vmesh = numPixels^2*Vmesh/numUnitCells^2;
    Vmesh = ifft2(Vmesh);
end

