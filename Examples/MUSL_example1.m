%%----------------Example 1----------------
%   Run the multislice program without changing any of the defaults
%       Silicon(001) with a 200keV beam


%% User Inputs

% lattice constants along each dimension [angstroms]
objLatticeInfo.Lattice = [5.431; 5.431; 5.431]; 
% locations of the atoms in units of the lattice constants.
% we assume we can treat the unit cell as a single layer.
objLatticeInfo.Positions = [0,0,0;              
                            0.5,0.5,0;         
                            0.25,0.25,0;
                            0.75,0.75,0;
                            0.5,0,0;
                            0,0.5,0;
                            0.25,0.75,0;
                            0.75,0.25,0 ];

% atomic number of element at each point of the Positions matrix
objLatticeInfo.AtomicNum = repmat(14,size(objLatticeInfo.Positions,1),1);
objLatticeInfo.BeamEnergy = 200*10^3;   %accelerating potential [Volts]
objLatticeInfo.crystalThickness = 2000; %crystal thickness [angstroms]

AngleX = 0;     %Tilt of the sample along X [mrad]
AngleY = 0;     %Tilt of the sample along Y [mrad]

%% Run the program
%   running MUSL initializes the multislice program, while intensity
%   performs the multislice calculation

objMUSL = MUSL(objLatticeInfo);
[intVal,reflList] = objMUSL.intensity(AngleX,AngleY);

%% Table with the reflections and intensities
%

intensityData = table(reflList,permute(intVal,[1 3 2]));
intensityData.Properties.VariableNames = {'relList' 'intensity'};

%% Real space image

im = surf(abs(fftshift(ifft2(objMUSL.Wave))));
surf(abs(fftshift(ifft2(objMUSL.Wave))),'edgecolor','none')
