%%----------------Example 3----------------
%   Run the multislice program to determine the approximate 2 beam 
%   extinction length of Si(001) and compare it to the relativistically 
%   corrected dynamical calculated value found in "Electron Microscopy of 
%   Thin Crystals" (1965) by P.B. Hirsch et al. To try to force the 
%   two-beam condition in the multislice program, we must make some 
%   adjustments to the propagator function. This is handled internally.
%   From Hirsch et al, we must also take into account that the original
%   calculations did not include the Debye-Waller temperature factor.
%   For a beam at 100keV, tilting to the Bragg conditions of the following
%   reflections should yield the following extinction lengths
%       (220)   75.7nm
%       (400)   126.8nm
%       (440)   209.3nm
%       
%       Silicon(100) with a 100keV beam


%% User Inputs
objLatticeInfo.Lattice = [5.431; 5.431; 5.431]; %lattice constants along each dimension [angstroms]
objLatticeInfo.Positions = [0,0,0;              %locations of the atoms in units of the lattice constants.
                            0.5,0.5,0;          %   we assume we can treat the unit cell as a single layer.
                            0.25,0.25,0;
                            0.75,0.75,0;
                            0.5,0,0;
                            0,0.5,0;
                            0.25,0.75,0;
                            0.75,0.25,0 ];

%atomic number of element at each point of the Positions matrix
objLatticeInfo.AtomicNum = repmat(14,size(objLatticeInfo.Positions,1),1);
objLatticeInfo.BeamEnergy = 100*10^3;   %accelerating potential [Volts]
objLatticeInfo.crystalThickness = 4000; %crystal thickness [angstroms]

%miller indices of the reflection we will be examining
h = 4;
k = 4;

%% Initialize the program
%   To make the multislice algorithm behave as though it were propogating 
%   two beams, we use the keyword pair 'TwoBeamCondition', [h,k], where 
%   h and k are the miller indices of the spot we are interested in. 
%   As was mentioned above, we must remove the Debye-Waller
%   Temperature factor, which can be done by using the keyword pair
%   'Temperature',-1.

objMUSL = MUSL(objLatticeInfo,'BravaisLattice','diamond','Temperature',-1,'TwoBeamCondition',[h,k]);

%% Run the multislice intensity calculation
%

%x and y tilts calculated from the miller indices [mrad]
AngleX = asin(h./2.*objMUSL.getLambda./objLatticeInfo.Lattice(1))*10^3;     %Tilt of the sample along X [mrad]
AngleY = asin(k./2.*objMUSL.getLambda./objLatticeInfo.Lattice(1))*10^3;     %Tilt of the sample along Y [mrad]

[intVal,reflList] = objMUSL.intensity(AngleX,AngleY,'depth');

%% Table with the reflections and intensities
% Form a table with the data

intensityData = table(reflList,permute(intVal,[1 3 2]));
intensityData.Properties.VariableNames = {'relList' 'intensity'};

%% Plot the data: depth vs intensity

reflNum1 = objMUSL.findSpot([0,0]);   %(000) reflection
reflNum2 = objMUSL.findSpot([h,k]);   %(040) reflection
figure
depth = linspace(0,objLatticeInfo.crystalThickness,length(intensityData.intensity(reflNum1,:)));
plot(depth,intensityData.intensity(reflNum1,:))
hold on
plot(depth,intensityData.intensity(reflNum2,:))
hold off
xlabel('depth (nm)')
ylabel('intensity (a.u.)')
legendcell = {'(000)',sprintf('(%d%d0)',h,k)};
legend(legendcell)