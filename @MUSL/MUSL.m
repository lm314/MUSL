classdef MUSL  < handle
%MUSL calcates the wave function as it traverses a cubic crystal and can
%   be used to calculate the final intensity of the diffraction spots at
%   exit using the multislice method. The MUSL object has the following 
%   properties
%
%   Lattice: lattice constant [angstroms] [1x3]
%   Positions: locations of the atoms in units of the Lattice constants
%       [nx3]. Due to the nature of the fft2 used, atoms that appear on the
%       edge also appear on the opposite edges of the unit cell, i.e. an
%       atom on the left edge appears on the right edge in the same
%       position, an atom on one corner appears on all corners. Thus, car
%       must be taken to not include unneccesary particles. This can be
%       checked by plotting the real space projected potential to see if 
%       there are any spikes larger than the rest, which would indicate
%       doubled or quadrupled particles.
%   BravaisLattice: type of Bravais Lattice of sample. This determines with 
%       hkl indices are output, that are nonzero. Can be:
%           simple-cubic (sc) - any integer miller indices allowed
%           body-centered cubic (bcc) - h+k+l even
%           face-centered cubic (fcc) - h,k,l all odd or all even
%           diamond-fcc (diamond) - all odd, or all even with h+k+l=4n
%   AtomicNum: atomic number of element at each point of the Positions 
%       matrix[nx1]
%   BeamEnergy: accelerating potential [Volts] [1x1].
%   CrystalThickness: crystal thickness [angstroms] [1x1]
%   Temperature: temperature of the sample defaults to 295K [Kelvin] [1x1]
%   Trans: transmission function of the sample [no units] [mxm]
%   Wave: wavefunctions of the electrons as it passes through sample [no
%       units] [mxm]
%   Lambda: electron wavelength [angstroms]
%   Absorption: indicates wheteher absorption is to be included [boolean]
%       [1x1]
%   RotationCrystal: rotation of the crystal plane relative to the x/y axes
%       of tilt (ccw) [degrees] [1x1]
%   KMax: maximum reciprocal lattice vector; determines the number of beams 
%       included in the simulation via the potential/transmission function 
%       [angstroms^-1] [1x1]
%   partKMax: for use when partial coherence is included, i.e. beam with 
%       angular spread; 0 indicates there is no angular spread [angstroms^-1]
%       ParCoIndFind does not function when this value is greater than the
%       separation between diffraction spots. Assumes beam has Gaussian
%       spatial distribution and that the angular spread is axially 
%       symmetric. partKMax is the standard deviation.
%   partKExtent: Assumes beam has Gaussian spatial distribution. Is the
%       number of standard deviations out the beam extends. 
%   KGrid: the Kx and Ky values of the grid which fall within partKMax
%   xPixelSize: Pixel size in real space (x,y) [angstroms]
%   kPixelSize: Pixel size in inverse momentum space space [angstrom^-1]
%   NumLayers: number of layers in the unit cell; determined by the number
%       if unique l values in the Positions matrix [no units] [1x1]
%   ZSpacing: vertical spacing of the layer in the multislice method
%       [angstroms] [1x1]
% 	NumUnitCells: number of unit cells in each direction (should be a 
%       power of 2 so that it can evenly divide the number of pixels, 
%       excluding 2^0) [no units] [1X1]
%   NumPixels: %number of total pixels in along each direction (should be 
%       a power of 2 ie 2^n to increase speed of calculation) [no units]
%       [1x1]. Not currently implemented.
%   UseGPU: specfiy whether to use GPU or CPU
%   ReflList: lists the reflections as they are appearin the H,K index
%       matrices [cx2]
%   Ind: gives the linear indices in the wave matrix that are associated
%       with each reflection (can be greater than 1 when partial coherence
%       is included) [cx1...d]
%   K2reflection: gives the inverse space squared value for each
%       reflection, which is useful for applying a profile in the case of
%       partial coherence.
%   TwoBeamCondition: if this is nonempty, then it contains the hkl miller
%       indices of the reflection that will be used in applying the two
%       beam condition - changes the propogator function [2 or 3 x 1]

    properties
        Lattice  
        Positions
        BravaisLattice = 'sc'
        AtomicNum             
        BeamEnergy
        CrystalThickness
        RotationCrystal = 0;
        Temperature = 295;
        Trans
        Prop
        Wave
    end
    
    properties (Access = private)
        Lambda
        Absorption = false
        KMax = 3.5
        partKMax = 0
        partKExtent = 0
        KGrid
        xPixelSize
        kPixelSize
        NumLayers
        ZSpacing
        NumUnitCells = 1    
        NumPixels = 64        
        UseGPU = false     
        ReflList
        Ind
        K2reflection
        TwoBeamCondition = []
    end
    
    methods
        function obj = MUSL(objLatticeInfo,varargin)
        %MUSL setups up the parameters of the multislice simulation. 
        %   Calculates the projected potential as well as the transmission 
        %    function. Vargin allows the user to pass in pairs consisting
        %    of a string that is the property name to be changed and the
        %    new value, i.e. 
        %       'Absorption',false
        %       'NumUnitCells',2
        %       'NumPixels', 128
        %

            %makes changes based on varargin
            if(mod(length(varargin),2)==0)
                for ii = 1:2:(length(varargin))
                    obj.propertyChange(varargin{ii},varargin{ii+1})
                end
            else
                error('Varargin does not have a matching pair: one of the property specfiers is missing a new value')
            end
            
            obj.Lattice = objLatticeInfo.Lattice;
            obj.Positions = objLatticeInfo.Positions;
            obj.AtomicNum = objLatticeInfo.AtomicNum;
            obj.BeamEnergy = objLatticeInfo.BeamEnergy;
            obj.CrystalThickness = objLatticeInfo.crystalThickness;
            
            obj.NumLayers = length(unique(obj.Positions(:,3))); %number of unique z positions
            obj.ZSpacing = obj.Lattice(3)/obj.NumLayers;
            obj.Lambda = 12.264306/(sqrt(obj.BeamEnergy)*sqrt(1+0.97846707*10^(-6)*obj.BeamEnergy));
            
            %setup parameters
            obj.kPixelSize = 1./obj.Lattice(1:2)./obj.NumUnitCells;
            obj.xPixelSize = obj.Lattice(1:2)*obj.NumUnitCells./obj.NumPixels;
            
            %Calculate the projected potential for each layer, which will be used to
            %calculate the transmission function. Only beams that fall within kMax are
            %used in calculating the final real space potential
            V = nan(obj.NumPixels,obj.NumPixels,obj.NumLayers);
            uniqueLayerVal = unique(obj.Positions(:,3));
            for ii = 1:obj.NumLayers
                condLayer = obj.Positions(:,3) == uniqueLayerVal(ii);
                V(:,:,ii) = MUSL.CrystalPotentialMax(obj.Positions(condLayer,:),obj.Lattice,obj.NumPixels,obj.NumUnitCells,obj.AtomicNum(condLayer,:),obj.BeamEnergy,obj.KMax,obj.Absorption,obj.Temperature);
            end
            
            %find the matrix indices that correspond to each reflection,
            %with or without an initial angular spread.
            obj.ParCoIndFind();      
                         
            %Calculate the transmission function
            obj.Trans = MUSL.transmission(V,obj.BeamEnergy,obj.ZSpacing);
            
            %check that the mesh is fine enough by checking the magnitude
            %of the transmission functions for each layer
            obj.checkFinenessMesh();
        end
        
        function [intVal,reflList] = intensity(obj,angleX,angleY,varargin)
            % Calculates the intensity of the diffraction pattern at exit
            % from the sample. angleX and angleY are in mrad
            % set varargin{1} = 'depth' to record the intensity of all spots
            % at each step
            if(isempty(varargin))
                    intVsDepth = false;
            else
                if(strcmp(varargin{1},'depth'))
                    intVsDepth = true;
                end
            end
                
            %check that angleX and angleY are less than or equal to 1
            %degree where the propogator approximation holds
            %convert 
            if(rad2deg(angleX*10^-3)>1 || rad2deg(angleY*10^-3)>1)
                error('Tilt angles are larger than 1 degree: approximation fails')
            end
            
            %initialize wave
            obj.initWave();
            
            %perform rotation of the crystal plane relative to x/y axes
            [angleX,angleY] = obj.rotateAngles(angleX,angleY);
            
            depthValues = 0:floor(obj.CrystalThickness/obj.ZSpacing);
            
            if(intVsDepth)
                intVal = zeros(size(obj.Ind,1),size(obj.Ind,2),length(depthValues));
                intVal(1,1,1) = 1;
            else
                intVal = nan(size(obj.Ind,1),size(obj.Ind,2),1);
            end
            %Calculate the propagator
            obj.Prop = MUSL.PropMesh(obj.Lattice, obj.NumUnitCells, obj.NumPixels, obj.ZSpacing, obj.BeamEnergy, angleX, angleY);
            
            %apply the two beam condition if spots are specified
            if(~isempty(obj.TwoBeamCondition))
                obj.Prop = obj.TwoBeamReduction(obj.Prop);
            end
            
            if(obj.UseGPU)
                gpuWave  = gpuArray(complex(obj.Wave));   %pass array to gpu, making sure to indicate the values may be complex
                gpuTrans = gpuArray(complex(obj.Trans));  %pass arrays to gpu to save time later
                gpuProp  = gpuArray(complex(obj.Prop));
                
                %multislice calculation with partial coherence - weight will
                %be applied later
                for slice = depthValues(2:end)
                    if(obj.NumLayers ==1)   %only 1 layer
                        gpuWave = pagefun(@times, gpuTrans, gpuWave); %multiply real space transmission function and wave
                    else
                        gpuWave = pagefun(@times, gpuTrans(:,:,rem(slice-1,obj.NumLayers)+1), gpuWave);
                    end
                    gpuWave = fft2(gpuWave);
                    gpuWave = pagefun(@times,gpuWave,gpuProp);
                    gpuWave = ifft2(gpuWave);
                    if(intVsDepth)
                        intVal(:,:,slice+1) = obj.calcInt(gpuWave);
                    end
                end
                if(~intVsDepth)
                    obj.Wave = gather(gpuWave);
                    intVal(:,:) = obj.calcInt(gpuWave); %find intensity for listed reflections
                end
            else    % No gpu       
                for slice = depthValues(2:end)
                    if(obj.NumLayers ==1)   %only 1 layer
                        obj.Wave = bsxfun(@times, obj.Trans, obj.Wave);
                    else
                        obj.Wave = bsxfun(@times, obj.Trans(:,:,rem(slice-1,obj.NumLayers)+1), obj.Wave);
                    end
                    obj.Wave = fft2(obj.Wave);
                    obj.Wave = bsxfun(@times,obj.Wave,obj.Prop);
                    obj.Wave = ifft2(obj.Wave);
                    if(intVsDepth)
                        intVal(:,:,slice+1) = obj.calcInt(obj.Wave);
                    end
                end
                if(~intVsDepth)
                    intVal(:,:) = obj.calcInt(obj.Wave);
                end
            end
            reflList = obj.ReflList;
        end
        
        function initWave(obj)
            [Kx,Ky] = obj.gridInvSpace;
            [X,Y] = obj.gridRealSpace;
            K2 = Kx.^2+Ky.^2;
            %create a 2 by ~ vector containing the momentum space coordinates of
            %the k values less than kmax in the momentum space mesh            
            obj.KGrid = [ Kx(K2 <= (obj.partKExtent*obj.partKMax)^2), Ky(K2 <= (obj.partKExtent*obj.partKMax)^2)];
            
            %Create numPixels by numPixels by size(kScan,1) matrix to be filled with
            %the initialized wave functions
            obj.Wave = zeros(obj.NumPixels,obj.NumPixels,size(obj.KGrid,1));
            
            %initialize wave function in real space, making sure to normalize it
            for kValues = 1:size(obj.KGrid,1)
                obj.Wave(:,:,kValues) = 1/(obj.NumPixels^2).*exp(2*pi*1j*(X.*obj.KGrid(kValues,1) + Y.*obj.KGrid(kValues,2)));
            end
        end
        
        %---------------set methods----------
        function set.BravaisLattice(obj,val)
            MUSL.mustBeMember(val,{'sc','fcc','bcc','diamond'});
            obj.BravaisLattice = val;
        end
        
        function set.Lattice(obj,val)
            MUSL.mustBePositive(val);
            obj.Lattice = val;
        end
        
        function set.AtomicNum(obj,val)
            MUSL.mustBePositive(val);
            MUSL.mustBeInteger(val);
            obj.AtomicNum = val;
        end
        
        function set.BeamEnergy(obj,val)
            MUSL.mustBePositive(val);
            MUSL.mustBeScalar(val);
            obj.BeamEnergy = val;
        end
        
        function set.CrystalThickness(obj,val)
            MUSL.mustBeNonnegative(val);
            obj.CrystalThickness = val;
        end
        
        function set.Lambda(obj,val)
            MUSL.mustBeNonnegative(val);
            MUSL.mustBeScalar(val);
            obj.Lambda = val;
        end
        
        function set.Absorption(obj,val)
            MUSL.mustBeScalar(val);
            MUSL.mustbeLogical(val);
            obj.Absorption = val;
        end
        
        function set.RotationCrystal(obj,val)
            MUSL.mustBeScalar(val);
            MUSL.mustBeReal(val);
            obj.RotationCrystal = val;
        end
        
        function set.KMax(obj,val)
            MUSL.mustBeReal(val);
            MUSL.mustBeScalar(val);
            obj.KMax = val;
        end
        
        function set.partKMax(obj,val)
            MUSL.mustBeReal(val);
            MUSL.mustBeScalar(val);
            obj.partKMax = val;
        end
        
        function set.xPixelSize(obj,val)
            MUSL.mustBePositive(val);
            obj.xPixelSize = val;
        end
        
        function set.kPixelSize(obj,val)
            MUSL.mustBePositive(val);
            obj.kPixelSize = val;
        end
        
        function set.NumLayers(obj,val)
            MUSL.mustBeInteger(val);
            MUSL.mustBePositive(val);
            MUSL.mustBeScalar(val);
            obj.NumLayers = val;
        end
        
        function set.ZSpacing(obj,val)
            MUSL.mustBePositive(val);
            MUSL.mustBeScalar(val);
            obj.ZSpacing = val;
        end
        
        function set.NumUnitCells(obj,val)
            MUSL.mustBePositive(val);
            MUSL.mustBeInteger(val);
            MUSL.mustBeScalar(val);
            obj.NumUnitCells = val;
        end
        
        function set.NumPixels(obj,val)
            MUSL.mustBePositive(val);
            MUSL.mustBeInteger(val);
            MUSL.mustBeScalar(val);
            obj.NumPixels = val;
        end
        
        function set.UseGPU(obj,val)
            MUSL.mustBeLogical(val)
            if(val==true)
                obj.UseGPU = true;
                %check to see if GPU is present
                gpuPresent = (gpuDeviceCount >= 1);
                if(gpuPresent == 1)
                    gpuArray(ones(1));
                    %deals with odd error that apears for large matrices - turns off memory pool
                    feature('GpuAllocPoolSizeKb', 0);   
                    %https://www.mathworks.com/matlabcentral/answers/167007-any-way-to-configure-gpu-based-fft-library-fft-ifft-to-free-memory-after-call
                end
            end
        end
        
        function set.TwoBeamCondition(obj,val)
            if(~isempty(val))
                MUSL.mustBeInteger(val)
                obj.TwoBeamCondition = val;
            else
                obj.TwoBeamCondition = val;
            end
        end
        
        %--------get methods----------------
        function lambda = getLambda(obj)
            lambda = obj.Lambda;
        end
        
        function rl = getReflList(obj)
            rl = obj.ReflList;
        end
        
        function zSpacing = getZSpacing(obj)
            zSpacing = obj.ZSpacing;
        end
        
        %--------calculate mesh details-------
        function x = GridPoints(obj)
            %vector used to assemble final real and momentum space grids with zero at center
            x = ((0:(obj.NumPixels-1)) - obj.NumPixels/2); 
        end
        
        function [H,K] = gridIndices(obj)
            %grid of the indices
            [H,K] = meshgrid(obj.GridPoints/obj.NumUnitCells);
        end
        
        function [X,Y] = gridRealSpace(obj)
            %grid in real space [angstroms]
            [X,Y] = meshgrid(obj.xPixelSize(1).*obj.GridPoints,obj.xPixelSize(2).*obj.GridPoints);
        end
        
        function [Kx,Ky] = gridInvSpace(obj)
            %grid in inverse space [angstroms^-1]
            [Kx,Ky] = meshgrid(obj.kPixelSize(1).*obj.GridPoints,obj.kPixelSize(2).*obj.GridPoints);
        end
        
        function idx = findSpot(obj,indices)
        % returns the spots linear index in the reflList given the 
        % miller indices
            [q, idx] = ismember(indices, obj.ReflList,'rows');
            if(any(q==0))
                error('The specified Miller indices could not be found in the list of reflections')
            end
        end
             
        function plotMap(obj,xAngleScale,yAngleScale,intVal,millerInd,fontSize,varargin)
        %Plots intensity values for matrices of x and y tilt combinations. 
        % intVal is the output of the MUSL.intensity command, with the 
        % reflection we wish to plot for being the millerInd.  
        % varargin is an array of the spots to normalize the simulations
        
            KX = xAngleScale*1e-3/obj.Lambda;
            KY = yAngleScale*1e-3/obj.Lambda;
            if(isa(intVal,'cell'))
                IntMap = obj.getIntMap(intVal,millerInd,varargin{:});
            elseif(isa(intVal,'double'))
                IntMap = intVal;
            else
                error('intVal is an unexpected data type.\n')
            end
            xlimitValues = [min(KX(1,:)) max(KX(1,:))];
            ylimitValues = [min(KY(:,1)) max(KY(:,1))];
            set(gcf,'defaultTextInterpreter','latex');
            
            surf(KX,KY,100*IntMap,'edgecolor','none');
            set(gca,'FontSize',fontSize,'Layer','top')
            grid off
            axis('equal')
            xlabel('$K_x (ang^{-1})$','FontSize', fontSize)
            ylabel('$K_y (ang^{-1})$','FontSize', fontSize)
            xlim(xlimitValues);
            ylim(ylimitValues);
            view(0,90)
        end
        
        function plotMapEPS(obj,xAngleScale,yAngleScale,intVal,millerInd,fontSize,varargin)
        % The same as plotMap, however, the output is optimized to look
        % correct when the figure is saved as a vector image, i.e. as an eps
        % image.
            KX = xAngleScale*1e-3/obj.Lambda;
            KY = yAngleScale*1e-3/obj.Lambda;
            if(isa(intVal,'cell'))
                IntMap = obj.getIntMap(intVal,millerInd,varargin{:});
            elseif(isa(intVal,'double'))
                IntMap = intVal;
            else
                error('intVal is an unexpected data type.\n')
            end
            xlimitValues = [min(KX(1,:)) max(KX(1,:))];
            ylimitValues = [min(KY(:,1)) max(KY(:,1))];
            set(gcf,'defaultTextInterpreter','latex');
            %----changes to allow for eps------
            imshow(IntMap, 'XData', unique(KX), 'YData', unique(KY))
            set(gca,'FontSize',fontSize,'defaultTextInterpreter','latex');
            set(gca,'YDir','normal')
            %-----------------------------------
            set(gca,'FontSize',fontSize,'Layer','top')
            grid off
            axis('equal')
            xlabel('$K_x (ang^{-1})$','FontSize', fontSize)
            ylabel('$K_y (ang^{-1})$','FontSize', fontSize)
            xlim(xlimitValues);
            ylim(ylimitValues);
            view(0,90)
        end
        
        function plotMapRealSpace(obj,xAngleScale,yAngleScale,intVal,millerInd,fontSize,varargin)
        %The same as plotMap but the printed axes correspond to the real 
        % space tilt angles. 
            KX = xAngleScale;
            KY = yAngleScale;
            if(isa(intVal,'cell'))
                IntMap = obj.getIntMap(intVal,millerInd,varargin{:});
            elseif(isa(intVal,'double'))
                IntMap = intVal;
            else
                error('intVal is an unexpected data type.\n')
            end
            xlimitValues = [min(KX(1,:)) max(KX(1,:))];
            ylimitValues = [min(KY(:,1)) max(KY(:,1))];
            set(gcf,'defaultTextInterpreter','latex');
            
            surf(KX,KY,100*IntMap,'edgecolor','none');
            set(gca,'FontSize',fontSize,'Layer','top')
            grid off
            axis('equal')
            xlabel('$\theta_x (mrad)$','FontSize', fontSize)
            ylabel('$\theta_y (mrad)$','FontSize', fontSize)
            xlim(xlimitValues);
            ylim(ylimitValues);
            view(0,90)
        end
        
        function plotIntVsDepth(obj,intVal,millerInd,fontSize)
        % Plots the insensity vs the thickness of the crystal. intVal is
        % the output of MUSL.intensity when the intVsDepth option is set to
        % true.
            depthValues = (0:floor(obj.CrystalThickness/obj.ZSpacing))*obj.ZSpacing/10;
            intValFiltered = obj.getIntVsDepth(intVal,millerInd);
            set(gcf,'defaultTextInterpreter','latex');
            plot(depthValues,intValFiltered)
            set(gca,'FontSize',fontSize,'Layer','top','TickLabelInterpreter','latex')
            xlabel('Depth (nm)','FontSize', fontSize)
            ylabel('Intensity','FontSize', fontSize)
        end
        
        function normConst = getNormConst(obj,intVal,varargin)
        % Returns the normalization constant to scale the intensity to 
        % include the specified reflections. intVal is an mxn cell array 
        % where each cell contains the output of MUSL.intensity. If no 
        % reflections are specified, then the normalization constant is 1 
        % for each orientation.    
            if(~isempty(varargin))
                spots = varargin{:};
                tempIntMap = nan(size(intVal,1),size(intVal,2),length(spots));
                for ii = 1:length(spots)
                    tempIntMap(:,:,ii) = obj.getIntMap(intVal,spots(ii,:))';
                end

                normConst = sum(tempIntMap,3)';
            else
                normConst = ones(size(intVal'));
            end
        end
        
        function intMap = getIntMap(obj,intVal,millerInd,varargin)
        % Returns the intensity map (the intensity matrix for a single 
        % reflection for all the crystallographic orientations found in 
        % intVal). Note that intVal is an mxn cell array where each cell 
        % contains the output of MUSL.intensity. varargin is an array of 
        % the spots to normalize intVal.
            
            normConst = getNormConst(obj,intVal,varargin{:});    

            ind = obj.findSpot(millerInd);
            intMap = nan(size(intVal'));
            if(~isvector(intVal{1}))
                for ii=1:numel(intVal)
                    intMap(ii) = sum(intVal{ii}(:,ind));
                end
            else
                for ii=1:numel(intVal)
                    intMap(ii) = sum(intVal{ii}(ind));
                end
            end
            intMap = intMap./normConst;
        end
        
        function intValFiltered = getIntVsDepth(obj,intVal,millerInd)
        % Returns the intensity vs depth for a single reflection from the
        % intVal produced by MUSL.intensity with the 'depth' option 
        % enabled.
        
            ind = obj.findSpot(millerInd);
            if(size(intVal,2)==1) %No partial coherence
                intValFiltered = squeeze(intVal(ind,:,:));
            else %partial coherence
                intValFiltered = squeeze(sum(intVal(:,ind,:),1));
            end
        end
    end

    methods (Access = private)
        function propertyChange(obj,inputProp,newVal)
            % Checks to see if the specified property is part of the MUSL
            % object, then changes it to the new value
            mco = ?MUSL;
            plist = mco.PropertyList;
            propFound = false;
            ii = 1;
            while(propFound~=true && ii<=length(plist))
                if(strcmp(inputProp,plist(ii).Name))
                    obj.(plist(ii).Name) = newVal;
                    propFound = true;
                end
                ii = ii + 1;
            end
        end
        
        function checkFinenessMesh(obj)
            % Checks to see if the mesh resolution needs to be increased by
            % checking to see if the absolute value of the transmission
            % function is sufficiently close to 1.
            if(obj.Absorption ~= true && isempty(obj.TwoBeamCondition))
                for i = 1:size(obj.Trans,3)
                    diff = max(abs(max(max(abs(obj.Trans(:,:,i))))^2-1),abs(min(min(abs(obj.Trans(:,:,i))))^2-1));
                    if( diff > 10^-5)
                        error(['\nThe mesh is not fine enough - add more pixels\n', ...
                            'or decrease the number of unit cells. Difference\n', ...
                            'of transmission function from 1 is %e \n'], diff);
                    end
                end
            end
        end
        
        function [reflList] = ParCoIndFind(obj)
            % PARCOINDFIND Calculates the matrices containing the H and K indices associated
            % with the momentum space (as arranged in the fft2 output along with the
            % maximum deflection vector (in inverse angstroms) of the beam,
            % so that it will determine the number of linear indices m that will be
            % associated with given reflections (2 by n of reflections matrix and return
            % an m by n matrix. lattice is in angstroms
                            
            %H and K space ordered to suit FFT2 output
            indSpace(1:(obj.NumPixels/2)) = ((1:(obj.NumPixels/2))-1);
            indSpace((obj.NumPixels/2+1):(obj.NumPixels)) = (((obj.NumPixels/2):(obj.NumPixels-1)) - obj.NumPixels);
            [H,K] = meshgrid(indSpace/obj.NumUnitCells);
            kPix = obj.kPixelSize;
            
            %Find the nonzero reflections
            condAllowedRefl = obj.AllowedReflections(H,K);
            reflList = [H(condAllowedRefl),K(condAllowedRefl)];
            
            if(obj.partKMax ~=0)
                indMax2 = (obj.partKExtent.*obj.partKMax./obj.kPixelSize'./obj.NumUnitCells).^2;   %maximum deflection index allowed

                %Find position of the reflections without deflection in momentum
                %space rearranged for FFT 2
                if(obj.partKMax ~=0)
                    numPosCond = H.^2./indMax2(1) + K.^2./indMax2(2) <= 1;
                else
                    ind2 = H.^2 + K.^2;
                    numPosCond = ind2 <= 0;
                end
                %spy(numPosCond)
                numPos = length(find(numPosCond));
                obj.Ind = NaN(numPos,size(reflList,1));
                obj.K2reflection = NaN(numPos,size(reflList,1));
                edgeValue = max(abs(reflList(:,1))); %magnitude of maximum which appears at edge
                clipMask = (abs(reflList(:,1))+sqrt(indMax2(1))>= edgeValue);
                clipVals = unique(abs(reflList(clipMask,1)));
                for refl = 1:size(reflList,1)
                    %dealing with edge cases where values come from repeat pattern
                    %(due to nature of fft)
                    if(any(abs(reflList(refl,1)) == clipVals))
                        indShiftH = fftshift(H) + sign(reflList(refl,1))*(abs(max(abs(clipVals))-abs(reflList(refl,1))));
                        if(reflList(refl,1)+max(abs(clipVals)) < 0)
                            indShiftH(indShiftH<-reflList(refl,1)) = indShiftH(indShiftH<-reflList(refl,1))+2*abs(reflList(refl,1));
                        end
                    else
                        indShiftH = H - reflList(refl,1);
                    end
                    if(any(abs(reflList(refl,2)) == clipVals))
                        indShiftK = fftshift(K) + sign(reflList(refl,2))*(abs(max(abs(clipVals))-abs(reflList(refl,2))));
                        if(reflList(refl,2)+max(abs(clipVals)) < 0)
                            indShiftK(indShiftK<-reflList(refl,2)) = indShiftK(indShiftK<-reflList(refl,2))+2*abs(reflList(refl,2));
                        end
                    else
                        indShiftK = K - reflList(refl,2);
                    end

                    ShiftedCond = indShiftH.^2./indMax2(1) + indShiftK.^2./indMax2(2) <= 1;
                    %spy(ShiftedCond)
                    obj.Ind(:,refl) = find(ShiftedCond);
                    obj.K2reflection(:,refl) =  obj.NumUnitCells^2*(indShiftH(obj.Ind(:,refl)).^2.*kPix(1)^2 + indShiftK(obj.Ind(:,refl)).^2.*kPix(2)^2);
                end
            else
                obj.Ind = find(condAllowedRefl);
                obj.K2reflection = zeros(size(obj.Ind));
            end
            obj.ReflList = reflList;
        end
        
        function [angleX,angleY] = rotateAngles(obj,angleX,angleY)
            % converts the x and y angles in mrad to their equivalent x and
            % y angles if the coordinate system was rotated by some amount
            % counter-clockwise
            
            H = obj.Lattice(1)/obj.Lambda*sin(angleX*10^-3);
            K = obj.Lattice(2)/obj.Lambda*sin(angleY*10^-3);
            A= [ H;K];
            rot = deg2rad(obj.RotationCrystal);
            A = [cos(rot), -sin(rot); sin(rot), cos(rot)]*A;
            H = reshape(A(1,:),size(H));
            K = reshape(A(2,:),size(K));
            
            angleX = 10^3*asin(obj.Lambda/obj.Lattice(1)*H); %tilt of the speciment along the x direction in mrad 0.81mrad for 2,2 spot
            angleY = 10^3*asin(obj.Lambda/obj.Lattice(2)*K); %same as above but along y direction        end
        end   
        
        function intVal = calcInt(obj,wavefunct)
             %Calculates the intensity from the wave function.
            temp = abs(fft2(wavefunct)).^2;   %back to reciprocal space
            % add ability to make beam profile gaussian before summing
            % intensity of layers
            if(obj.partKExtent ~=0)
                weight = exp(-obj.K2reflection./2./(obj.partKMax).^2); %normalization is handled later
                weight = permute(weight,[3 2 1]);
                temp = bsxfun(@times,temp,weight);
                %temp = bsxfun(@rdivide,temp,sum(sum(temp,1),2));
            end
            temp = sum(temp,3);             %compress layers for a single 2d image
            temp = temp/sum(sum(temp)); %2d normalization
            %fprintf('%f\n',sum(sum(temp)))
            intVal = gather(temp(obj.Ind)); %find intensity for listed reflections
        end
        
        function condAllowedRefl = AllowedReflections(obj,H,K)
        % Determines which indices are allowed reflections depending on the
        % Bravais lattice type. H,K are the miller indices we are checking
            switch obj.BravaisLattice
                case 'sc'
                    %only integers
                    condAllowedRefl = floor(H)==H & floor(K)==K;
                case 'fcc'
                    %h+k+l=even
                    condAllowedRefl = floor(H)==H & floor(K)==K & mod(H+K,2)==0;
                case 'bcc'
                    %h,k,l all odd or all even
                    isOdd = mod(H,2)==1 & mod(K,2)==1;
                    isEven = mod(H,2)==0 & mod(K,2)==0;
                    condAllowedRefl = isOdd | isEven;
                case 'diamond'
                    %all odd, or all even with h+k+l=4*n
                    %since we are only looking at the zeroth order Laue
                    %Zone (l=0) we get cannot get any odd
                    condAllowedRefl = mod(H,2)==0 & mod(K,2)==0 & mod(H+K,4)==0;            
                otherwise
                    condAllowedRefl = true(numel(H));
            end
        end
        
        function Prop = TwoBeamReduction(obj,Prop)
        % Removes all but the (000) and the specified spot from the
        % propogator function to mimic the two beam case
            spot1 = obj.findSpot(zeros(1,length(obj.TwoBeamCondition)));
            spot2 = obj.findSpot(obj.TwoBeamCondition);
            Vspot1temp = Prop(obj.Ind(spot1));
            Vspot2temp = Prop(obj.Ind(spot2));
            Prop = zeros(size(Prop));
            Prop(obj.Ind(spot1)) = Vspot1temp;
            Prop(obj.Ind(spot2)) = Vspot2temp;
        end
    end
    
    methods (Static)
        Vmesh = CrystalPotentialMax( unitCellSites, lattice, numPixels, numUnitCells,element, E, kMax, absorp, temperature)
        [allSites, allAtoms] = populateUnitCells( sites,atoms, numUnitCells, lattice )
        fe = ElectronScatFactor( q2, element )
        trans = transmission( V, E, zSpacing )
        prop = PropMesh( lattice, numUnitCells, numPixels, zSpacing, E, angleX, angleY  )
        fi = AbsorptivePotential( q2, E, temperature)
        B = DebyeWaller(element,temperature);
        function lambda = electronWavelength(accelVolt)
        % Returns the electron wavelength in angstroms. accelVolt is in
        % Volts.
            lambda = 12.264306./(sqrt(accelVolt).*sqrt(1+0.97846707*10^(-6).*accelVolt));
        end
        
        % The following functions act similar to the property validator
        % functions that became available with matlab 2017a
        function mustBePositive(A)
            if(any(~isreal(A)))
                error('Input must be real')
            else
                if(any(A <= 0))
                    error('Input must be greater than 0')
                end
            end      
        end
        function mustBeInteger(A)
            if(any(~isreal(A)))
                error('Input must be real')
            else
                if(any(floor(A)~=A))
                    error('Input must be integer')
                end
            end
        end
        function mustBeMember(A,B)
            if(~ismember(A,B))
                error('Input must be a member of given set');
            end
        end
        function mustBeNonnegative(A)
            if(any(~isreal(A)))
                error('Input must be real')
            else
                if(any(A<0))
                    error('Input must be non-negative')
                end
            end
        end
        function mustBeLogical(A)
            if(any(~islogical(A)))
                error('Error: must be logical')
            end
        end
        function mustBeScalar(A)
            if(length(A)>1)
                error('Input must be a scalar')
            end
        end
        function mustBeReal(A)
            if(any(~isreal(A)))
                error('Input must be real')
            end
        end
    end
end

