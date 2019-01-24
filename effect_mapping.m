%owleywhere compares significant effects from a multivariate analysis in GIFT (e.g. age) to (thresholded)
%neurosynth reverse inference maps or cogmaps (cogmaps do not differentiate positive/negative) voxels
%cogmaps must be resampled to same space as t-maps
%please cite our paper in Human Brain Mapping: de Lacy et al 2018

%%things to prep the script
%1. set the tmaps dir and workingmapname including the correct analysis prefix and effect type
%2. set the owleyoo save file name
%3. cd to correct sm_stats dir
%4. consider changing t-map threshold 

clear all; clc;

%set analysis parameters for sig fx covariate and cogmap threshold
tmapth = 1;
cogmapth = 0;

%%
%set basefiles: cogmap 
cogmaps = dir('/yourdirectory_of_cogmaps/*.nii');
nummaps = length(cogmaps);

%set basefiles: tmaps of significant effects
%tmaps = dir('match_low_mancovan_sm_gend_(1) - (0)_sig_effects_comp_0*.img');
tmaps = dir('effect_of_interest_*.img');
numfiles = length(tmaps);

%%
%create structs for each component of interest
for n=(1:numfiles);
    owleyoo(n).comp = struct();
    workingmapname = strrep(tmaps(n).name,'effect_of_interest','');
    owleyoo(n).comp.name = workingmapname;
end
    
%%
%extracts nonzero voxelwise z scores and location indices for each
%thresholded mancovan t map (positive and negative)
for n=(1:numfiles);
    A=spm_vol(tmaps(n).name);
    [A,ABC]=spm_read_vols(A);
    owleyoo(n).comp.nnz = nnz(A);
    u = find(A>0);
    owleyoo(n).comp.pos.voxelvalues = A(u);
    owleyoo(n).comp.pos.voxelindices = u;
    owleyoo(n).comp.pos.voxellocations = ABC(:,owleyoo(n).comp.pos.voxelindices);
    R=owleyoo(n).comp.pos.voxellocations;
    owleyoo(n).comp.pos.voxellocationscolumn = R';
    t = find(A<0);
    owleyoo(n).comp.neg.voxelvalues = A(t);
    owleyoo(n).comp.neg.voxelindices = t;
    owleyoo(n).comp.neg.voxellocations = ABC(:,owleyoo(n).comp.neg.voxelindices);
    Q=owleyoo(n).comp.neg.voxellocations;
    owleyoo(n).comp.neg.voxellocationscolumn = Q';
end
 

%%
%creates structs for each cogmap
for j=(1:nummaps);
    owleyoo(j).cogmap = struct();
    owleyoo(j).cogmap.name = cogmaps(j).name;
end
%%   
%extracts nonzero voxelwise z scores and location indices for each
%thresholded cogmap (positive and negative are not differentiated)
for j=(1:nummaps);
    fullfilepath = strcat('/yourdirectory_of_cogmaps/',cogmaps(j).name);
    V=spm_vol(fullfilepath);
    [V,XYZ]=spm_read_vols(V);
    V(V<cogmapth)=0;
    owleyoo(j).cogmap.nnz = nnz(V);
    w=find(V>0);
    owleyoo(j).cogmap.voxelvalues=V(w);
    owleyoo(j).cogmap.voxelindices=w;
end

%%
%computes indices overlap for every combination of tmap pos indices and cogmap
for n=(1:numfiles);
   for j=(1:nummaps);
   owleyoo(n).comp(j).overlappos(:) = intersect(owleyoo(j).cogmap.voxelindices,owleyoo(n).comp(1).pos.voxelindices);
   end
end

%%
%computes indices overlap for every combination of tmap neg indices and
%cogmaps
for n=(1:numfiles);
   for j=(1:nummaps);
   owleyoo(n).comp(j).overlapneg(:) = intersect(owleyoo(j).cogmap.voxelindices,owleyoo(n).comp(1).neg.voxelindices);
   end
end
  
%%
%computes mni locations in mm corresponding to pos overlap indices
for n=(1:numfiles);
    A=spm_vol(tmaps(n).name);
    [A,ABC]=spm_read_vols(A);
    for j=(1:nummaps);
    owleyoo(n).comp(j).locationspos = ABC(:,owleyoo(n).comp(j).overlappos);
    end
end

%%
%computes mni locations in mm corresponding to neg overlap indices
for n=(1:numfiles);
    A=spm_vol(tmaps(n).name);
    [A,ABC]=spm_read_vols(A);
    for j=(1:nummaps);
    owleyoo(n).comp(j).locationsneg = ABC(:,owleyoo(n).comp(j).overlapneg);
    end
end

%%
%computes common (intersection) voxel cluster for each cogmap of this t-map
%set pos locations
for j=(1:nummaps);
    F=[]; 
    for n=(1:numfiles); 
          F=[F,owleyoo(n).comp.overlappos];
    end
    owleyoo(j).cogmap.uniqueindicespos=unique(F);
end

%%
%computes common (intersection) voxel cluster for each cogmap of this t-map
%set neg locations
for j=(1:nummaps);
    F=[]; 
    for n=(1:numfiles); 
          F=[F,owleyoo(n).comp.overlapneg];
    end
    owleyoo(j).cogmap.uniqueindicesneg=unique(F);
end

%%
%computes mni locations in mm corresponding to unique pos locations (concatenated tmaps) for each cogmap
for j=(1:nummaps);
    fullfilepath = strcat('/yourdirectory_of_cogmaps/',cogmaps(j).name);
    V=spm_vol(fullfilepath);
    [V,XYZ]=spm_read_vols(V);
    owleyoo(j).cogmap.locationsuniquepos = XYZ(:,owleyoo(j).cogmap.uniqueindicespos);
    G=owleyoo(j).cogmap.locationsuniquepos;
    owleyoo(j).cogmap.locationsuniqueposcolumns = G';
end    

%%
%computes mni locations in mm corresponding to unique neg locations (concatenated tmaps) for each cogmap
for j=(1:nummaps);
    fullfilepath = strcat('/yourdirectory_of_cogmaps/',cogmaps(j).name);
    V=spm_vol(fullfilepath);
    [V,XYZ]=spm_read_vols(V);
    owleyoo(j).cogmap.locationsuniqueneg = XYZ(:,owleyoo(j).cogmap.uniqueindicesneg);
    H=owleyoo(j).cogmap.locationsuniqueneg;
    owleyoo(j).cogmap.locationsuniquenegcolumns=H';
end    



%%
%computes intersection between concatenated tmap unique pos indices and cogmap
%indices (this is the overlap between all the tmaps and the individual
%cogmap
for j=(1:nummaps);
    owleyoo(j).cogmap.overlapindicespos = intersect(owleyoo(j).cogmap.voxelindices,owleyoo(j).cogmap.uniqueindicespos)
end

%%
%computes intersection between concatenated tmap unique neg indices and cogmap
%indices (this is the overlap between all the tmaps and the individual
%cogmap
for j=(1:nummaps);
    owleyoo(j).cogmap.overlapindicesneg = intersect(owleyoo(j).cogmap.voxelindices,owleyoo(j).cogmap.uniqueindicesneg)
end

%%
%computes mni locations in mm corresponding to overlaplocations between unique pos concatenated tmaps and each cogmap
for j=(1:nummaps);
    fullfilepath = strcat('/yourdirectory_of_cogmaps/',cogmaps(j).name);
    V=spm_vol(fullfilepath);
    [V,XYZ]=spm_read_vols(V);
    owleyoo(j).cogmap.locationsoverlappos = XYZ(:,owleyoo(j).cogmap.overlapindicespos);
    J=owleyoo(j).cogmap.locationsoverlappos;
    owleyoo(j).cogmap.locationsoverlapposcolumns=J';
end 

%%
%computes mni locations in mm corresponding to overlaplocations between unique neg concatenated tmaps and each cogmap
for j=(1:nummaps);
    fullfilepath = strcat('/yourdirectory_of_cogmaps/',cogmaps(j).name);
    V=spm_vol(fullfilepath);
    [V,XYZ]=spm_read_vols(V);
    owleyoo(j).cogmap.locationsoverlapneg = XYZ(:,owleyoo(j).cogmap.overlapindicesneg);
    K=owleyoo(j).cogmap.locationsoverlapneg;
    owleyoo(j).cogmap.locationsoverlapnegcolumns=K';
end 

%%
%save output 
save('yourfilename.mat','owleyoo');cd ..
cd
