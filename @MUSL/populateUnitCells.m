function [ allSites, allAtoms ] = populateUnitCells( sites,atoms, numUnitCells, lattice )
%POPULATEUNITCELLS Given initial atomic sites in a single unit cell in the
%form of an array of dimension #sites by 2, outputs all possible atomic sites 
%within the number of unit cells given as an array scaled to the particular
%lattice. Note that by using the Discrete Fourier Transform in the form of
%FFT, we are assuming that the patterns in momentum space and real space
%are periodic, repeating forever in the x and y directions. This means we
%need to only include the atoms appearing along two of the outer edges,
%i.e. including atoms on the left and bottom edges (the corner but not the
%two outer extents of these sides). Also changes the atoms matrix to reflect
%the change from sites to allSites.
    if(numUnitCells > 1)
        tempSites = NaN(size(sites,1)*numUnitCells^2,3); %maximum dimension
        tempAtom = NaN(size(sites,1)*numUnitCells^2,1); %maximum dimension
        tempCounter = 1;
        for atom = 1:size(sites,1)
            for i = -(numUnitCells/2):(numUnitCells/2-1)
                for j = -(numUnitCells/2):(numUnitCells/2-1)
                    if((i~=0) || (j~=0))
                        %if((abs(sites(atom,1)+i)<=(numUnitCells/2)) && (abs(sites(atom,2)+j)<=(numUnitCells/2)))
                        tempSites(tempCounter,1) = sites(atom,1)+i;
                        tempSites(tempCounter,2) = sites(atom,2)+j;
                        tempSites(tempCounter,3) = sites(atom,3);
                        tempAtom(tempCounter) = atoms(atom);
                        tempCounter = tempCounter+1;
                    end
                end
            end
        end

        allSites = [sites;tempSites(all(~isnan(tempSites),2),:);]; %only add none NaN values to sites array
        allAtoms = [atoms;tempAtom(all(~isnan(tempSites),2),:);];
        allSites = lattice'.*allSites;
    else
        allSites = lattice'.*sites;
        allAtoms = atoms;
    end
end

