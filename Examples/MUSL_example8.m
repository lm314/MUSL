%%----------------Intensity Maps----------------
%   Run the MUSL multislice program to generate an 
%   intensity Map for 200 nm thick Si(001) with a 
%   2.26 MeV electron beam

%% User Inputs
% lattice constants along each dimension [angstroms]
objLatticeInfo.Lattice = [5.431; 5.431; 5.431]; 
% Locations of the atoms in units of the lattice constants
% We treat the unit cell as a single layer.
objLatticeInfo.Positions = [0,0,0;              
                            0.5,0.5,0;          
                            0.25,0.25,0;
                            0.75,0.75,0;
                            0.5,0,0;
                            0,0.5,0;
                            0.25,0.75,0;
                            0.75,0.25,0 ];

% Atomic number of element at each point
objLatticeInfo.AtomicNum = repmat(14,...
            size(objLatticeInfo.Positions,1),1);
% Accelerating potential [Volts]
objLatticeInfo.BeamEnergy = 2.26*10^6;   
% Crystal thickness [angstroms]
objLatticeInfo.crystalThickness = 2000; 

% Tilt of the sample along X [mrad]
AngleX = linspace(-2.75,3.24,70);  
% Tilt of the sample along Y [mrad]
AngleY = linspace(10.75,13.72,35);    

[X,Y] = meshgrid(AngleX,AngleY);

%% Calculating Intensity as a Function of Orientation
%   Running MUSL initializes the multislice program, while 
%   intensity performs the multislice calculation. At 
%   initialization, we also included the keyword pair 
%   'BraviasLattice' and 'diamond' to reduce the number of 
%   reflections displayed. The keyword 'partKMax' 
%   indicates the RMS spread of the reciprocal lattice 
%   vector in inverse angstroms, assuming the beam has 
%   a Gaussian distribution. The keyword 'partKExtent' 
%   indicates maxmimum extent of the simulated spread 
%   in units of standard deviations.

objMUSL = MUSL(objLatticeInfo,...
                'BravaisLattice','diamond',...
                'RotationCrystal',-45,...
                'partKMax',0.0133,...
                'partKExtent',3);
intVal = cell(length(AngleX),length(AngleY));

for ii=1:numel(X)
    [intVal{ii},~] = objMUSL.intensity(X(ii),Y(ii));
end

%% Plotting Intensity Maps

% (000) reflection
fontSize = 12;
figure('Position', [800 100 500 300])
objMUSL.plotMap(X,Y,intVal,[0,0],fontSize)

% (2-20) reflection
figure('Position', [800 100 500 300])
objMUSL.plotMap(X,Y,intVal,[2,-2],fontSize)