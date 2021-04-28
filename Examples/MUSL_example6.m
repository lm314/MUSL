%%----------------Example 6----------------
%   Run the multislice program to create an intensity Map
%       Gold(001) with a 2.5 MeV beam


%% User Inputs
objLatticeInfo.Lattice = [4.065; 4.065; 4.065]; %lattice constants along each dimension [angstroms]
objLatticeInfo.Positions = [0,0,0;              %locations of the atoms in units of the lattice constants.
                            0.5,0.5,0;          %   we assume we can treat the unit cell as a single layer.
                            0.5,0,0;
                            0,0.5,0];

%atomic number of element at each point of the Positions matrix
objLatticeInfo.AtomicNum = repmat(79,size(objLatticeInfo.Positions,1),1);
objLatticeInfo.BeamEnergy = 2.5*10^6;   %accelerating potential [Volts]
objLatticeInfo.crystalThickness = 2500; %crystal thickness [angstroms]

extent = 0.25;   %degrees
AngleX = linspace(-deg2rad(extent),deg2rad(extent),81)*10^3;     %Tilt of the sample along X [mrad]
AngleY = linspace(-deg2rad(extent),deg2rad(extent),81)*10^3;     %Tilt of the sample along Y [mrad]

[X,Y] = meshgrid(AngleX,AngleY);

%% Run the program
%   running MUSL initializes the multislice program, while intensity
%   performs the multislice calculation.

objMUSL = MUSL(objLatticeInfo,'BravaisLattice','fcc','RotationCystal',-16.4);
intVal = cell(length(AngleX),length(AngleY));

parfor ii=1:numel(X)
    [intVal{ii},~] = objMUSL.intensity(X(ii),Y(ii));
end


%%
fontSize = 12;

spots = [0,0;
         -2,2;
         0,2;
         2,2;
         -2,0;
         2,0;
         -2,-2;
         0,-2;
         2,-2];
     
for ii = 1:size(spots,1)
    figure('Position', [800 100 500 300])
    objMUSL.plotMap(X,Y,intVal,spots(ii,:),fontSize)
    colorbar
    print(sprintf('plot%d%d0.png',spots(ii,1),spots(ii,2)),'-dpng','-r200')
end