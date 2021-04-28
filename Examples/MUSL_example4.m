%%----------------Ultramicroscopy Runs----------------
%   Run the multislice program to create an intensity Map
%       200nm thickness Silicon(100) with a 2.26 MeV beam 


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
objLatticeInfo.BeamEnergy = 2.26*10^6;   %accelerating potential [Volts]
objLatticeInfo.crystalThickness = 2000; %crystal thickness [angstroms]

AngleX = linspace(-2.75,3.24,70);     %Tilt of the sample along X [mrad]
AngleY = linspace(10.75,13.72,35);    %Tilt of the sample along Y [mrad]

[X,Y] = meshgrid(AngleX,AngleY);

%% Run the program
%   running MUSL initializes the multislice program, while intensity
%   performs the multislice calculation.

%objMUSL = MUSL(objLatticeInfo,'BravaisLattice','diamond','RotationCrystal',-45);
objMUSL = MUSL(objLatticeInfo,...
                'BravaisLattice','diamond',...
                'RotationCrystal',-45,...
                'partKMax',0.0133,...
                'partKExtent',3,...
                'NumUnitCells',8,...
                'NumPixels',512,...
                'UseGPU',true);
intVal = cell(length(AngleX),length(AngleY));

%timer to check progress of mutlslice step
tStart = tic;

for ii=1:numel(X)
    if(mod(ii,100)==0) %progress displayed
        tElapsed = toc(tStart);
        fprintf('On linear index %d of %d. Elapsed time: %0.1f minutes \n',...
            ii,numel(X),tElapsed/60)
    end
    [intVal{ii},~] = objMUSL.intensity(X(ii),Y(ii));
end

%% Plots of Intensity Maps
fontSize = 12;
figure('Position', [800 100 500 300])
objMUSL.plotMap(X,Y,intVal,[0,0],fontSize)
colorbar

print('intensityMapExample000','-dpng','-r150')

figure('Position', [800 100 500 300])
objMUSL.plotMap(X,Y,intVal,[2,-2],fontSize)
colorbar

print('intensityMapExample220','-dpng','-r150')

%save(sprintf('intensity_%s.mat',datestr(datetime('now'),'yy-mm-dd_HH_MM_SS')))