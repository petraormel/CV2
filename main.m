% fname = 'Data/data/0000000001.pcd';

% data = readPcd(fname);


load('source.mat');
load('target.mat');

icp(source,target);
