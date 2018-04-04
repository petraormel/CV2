fname = 'Data/data/0000000001.pcd';

% data = readPcd(fname);


load('Data/source.mat');
load('Data/target.mat');

icp(source,target);
