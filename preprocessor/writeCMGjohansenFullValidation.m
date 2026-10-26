
mrstModule add deckformat;

% Send to output folder NorneforCMG
outpth = fullfile('C:\Users\cigweo1\Desktop\mrst-2024a_Geomech_v3\examples\data', 'JohansenOUT2CMG') ; %will output to OUT2CMG folder in \mrst-2021b\examples\data\ folder
% outpth = fullfile('C:\Users\cigweo1\Desktop\mrst-2024a_Geomech_v3\examples\data\NorneforCMGOMO', 'NorneforCMGOMO') ;
%% Read deck for Johansen
jdir = fullfile(mrstPath('co2lab-common'), 'data', 'johansen');
sector = 'NPD5';
fprintf('    Reading %s\n', sector);
sector = fullfile(jdir, sector);
grdecl = readGRDECL([sector '.grdecl']);

G = processGRDECL(grdecl);
G = computeGeometry(G);
plotGrid(G);


% Load permeability and porosity
K = reshape(load([sector, '_Permeability.txt'])', [], 1);
p = reshape(load([sector, '_Porosity.txt'])',     [], 1);
grdecl.PERMX=K;
grdecl.PERMY=K;
grdecl.PERMZ=0.1*K;
grdecl.PORO=p;


%% Write out data for CMG
fid = fopen([outpth 'ZCORN.GRID'], 'w'); 
fprintf(fid, '*ZCORN ** %d float\n',  numel(grdecl.ZCORN));
fprintf(fid, '%0.8f %0.8f %0.8f %0.8f %0.8f %0.8f\n', grdecl.ZCORN(:));fprintf(fid, '\n');
fclose(fid);

fid = fopen([outpth 'COORD.GRID'], 'w'); 
fprintf(fid, '*COORD ** %d float\n', numel(grdecl.COORD));
fprintf(fid, '%0.8f %0.8f %0.8f %0.8f %0.8f %0.8f\n', grdecl.COORD(:));fprintf(fid, '\n');
fclose(fid);

fid = fopen([outpth 'ACTNUM.GRID'], 'w'); 
fprintf(fid, '*NULL *ALL ** %d float\n', numel(grdecl.ACTNUM));
fprintf(fid, '%d %d %d %d %d %d\n', grdecl.ACTNUM(:));fprintf(fid, '\n');
fclose(fid);

fid = fopen([outpth 'PERMI.txt'], 'w'); 
fprintf(fid, '*PERMI *ALL ** %d float\n', numel(grdecl.PERMX));
fprintf(fid, '%d %d %d %d %d %d\n', grdecl.PERMX(:));fprintf(fid, '\n');
fclose(fid);

fid = fopen([outpth 'PERMJ.txt'], 'w'); 
fprintf(fid, '*PERMI *ALL ** %d float\n', numel(grdecl.PERMY));
fprintf(fid, '%d %d %d %d %d %d\n', grdecl.PERMX(:));fprintf(fid, '\n');
fclose(fid);

fid = fopen([outpth 'PERMK.txt'], 'w'); 
fprintf(fid, '*PERMI *ALL ** %d float\n', numel(grdecl.PERMZ));
fprintf(fid, '%d %d %d %d %d %d\n', grdecl.PERMX(:));fprintf(fid, '\n');
fclose(fid);

fid = fopen([outpth 'POR.txt'], 'w'); 
fprintf(fid, '*POR *ALL ** %d float\n', numel(grdecl.PORO));
fprintf(fid, '%d %d %d %d %d %d\n', grdecl.PORO(:));fprintf(fid, '\n');
fclose(fid);


topLayer_global = false(G.cartDims);
topLayer_global(1:end, 1:end, 1) = true;
topLayer = find(topLayer_global(G.cells.indexMap));

rock.perm = [grdecl.PERMX, grdecl.PERMY, grdecl.PERMZ];

rock.perm(topLayer, 1:2) = rock.perm(topLayer, 1:2)*1e4;
k_mD = rock.perm(topLayer, :) ./ 9.869233e-16;




