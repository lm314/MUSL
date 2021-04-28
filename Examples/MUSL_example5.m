%%----------------Example 5----------------
%   Run the multislice program to create an intensity Map
%       Silicon(001) with a 2 MeV beam


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
objLatticeInfo.BeamEnergy = 2*10^6;   %accelerating potential [Volts]
objLatticeInfo.crystalThickness = 3400; %crystal thickness [angstroms]

AngleX = linspace(0,deg2rad(1)/8,161)*10^3;     %Tilt of the sample along X [mrad]
AngleY = linspace(0,deg2rad(1)/2,161)*10^3;     %Tilt of the sample along Y [mrad]

[X,Y] = meshgrid(AngleX,AngleY);

%% Run the program
%   running MUSL initializes the multislice program, while intensity
%   performs the multislice calculation.

objMUSL = MUSL(objLatticeInfo,'BravaisLattice','diamond','RotationCystal',45);
intVal = cell(length(AngleX),length(AngleY));

for ii=1:numel(X)
    [intVal{ii},~] = objMUSL.intensity(X(ii),Y(ii));
end

%%
fontSize = 12;
figure('Position', [800 100 500 300])
objMUSL.plotMap(X,Y,intVal,[0,0],fontSize)
colorbar
print('plot000.png','-dpng','-r200')

figure('Position', [800 100 500 300])
objMUSL.plotMap(X,Y,intVal,[2,2],fontSize)
colorbar
print('plot220.png','-dpng','-r200')