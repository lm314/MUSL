%%----------------Example 7----------------
%   Run the multislice program for Gold(100) with a 2.5meV beam in a normal
%   geometry for various temperatures.


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

%tempVals = [65,70,75,77,80,90,100,110,120,150,200,250,290,293,295,320,400,500]; %temperatures [K]
tempVals = [250,290,293,295,320]; %temperatures [K]

%% Run the program
%   running MUSL initializes the multislice program, while intensity
%   performs the multislice calculation.

intVal = cell(numel(tempVals),1);

for ii=1:numel(tempVals)
    objMUSL = MUSL(objLatticeInfo,'BravaisLattice','fcc','RotationCystal',-16.4,'Temperature',tempVals(ii));
    [intVal{ii},~] = objMUSL.intensity(0,0);
end


%%

figure('Position', [800 100 700 700])
spots = [-2,2;
         0,2;
         2,2;
         -2,0;
         2,0;
         -2,-2;
         0,-2;
         2,-2];
ind = objMUSL.findSpot(spots);
IntMap = nan(numel(tempVals),size(spots,1));
for ii=1:numel(tempVals)
    IntMap(ii,:) = intVal{ii}(ind)';
end
IntMap = reshape(IntMap,length(tempVals),size(spots,1));
IntMap = IntMap./sum(IntMap,2);

for ii = 1:size(spots,1)
    figure('Position', [800 100 700 700])
    plot(tempVals,IntMap(:,ii))
    title(sprintf('(%d%d0)',spots(ii,1),spots(ii,2)))
    xlabel('T (K)')
    ylabel('Y tilt (mrad)')
    %caxis([0,maxVal])
    view(0,90)
    %print(sprintf('plot(%d%d0).png',spots(ii,1),spots(ii,2)),'-dpng','-r200')
end