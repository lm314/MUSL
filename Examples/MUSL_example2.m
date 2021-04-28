%%----------------Example 2----------------
%   Run the multislice program to see intensity vs depth. Since silicon is
%   diamond cubic, keeping all the reflections as is in the defualt case of
%   a simple cubic lattice, is not necesary, so we use another feature of
%   MUSL at initialization to return only the allowed reflections.
%       Silicon(001) with a 200keV beam


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
objLatticeInfo.BeamEnergy = 200*10^3;   %accelerating potential [Volts]
objLatticeInfo.crystalThickness = 2000; %crystal thickness [angstroms]

AngleX = 0;     %Tilt of the sample along X [mrad]
AngleY = 0;     %Tilt of the sample along Y [mrad]

%% Run the program
%   running MUSL initializes the multislice program, while intensity
%   performs the multislice calculation. Note that we have added the
%   keyword 'depth' to the input of the intensity function. At 
%   initialization, we also included the keyword pair 'BraviasLattice' and 
%   'diamond' to reduce the number of reflections displayed.

objMUSL = MUSL(objLatticeInfo,'BravaisLattice','diamond');
[intVal,reflList] = objMUSL.intensity(AngleX,AngleY,'depth');

%% Table with the reflections and intensities
%

intensityData = table(reflList,permute(intVal,[1 3 2]));
intensityData.Properties.VariableNames = {'relList' 'intensity'};

%% Plot the data

reflNum = 1;

figure
depth = linspace(0,objLatticeInfo.crystalThickness,length(intensityData.intensity(reflNum,:)));
plot(depth,intensityData.intensity(reflNum,:))
xlabel('depth (nm)')
ylabel('intensity (a.u.)')