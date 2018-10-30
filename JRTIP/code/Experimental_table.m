%%%%%%%% PVO�ĸĽ�汾��Xiaolong Li��2018-07-22

tic
clear all;
clc

addpath(genpath('Origin Images')); 
addpath(genpath('tools')); 
addpath(genpath('result')); 
Imgs = {'Lena', 'Baboon', 'Airplane', 'Lake', 'Peppers', 'Boat', 'Barbara', 'Elaine'};

maps = {};
map = cell(3);
mapECs = {};
mapEDs = {};

[mapEDs, mapECs, maps] = getOneMap(mapEDs, mapECs, maps,map,[0,0]);
%%
sumED = 0;
for tt = 1:8
    Iname = Imgs{tt};
    
    istr = ['Proposed_2019_',Iname,'.mat'];
    
    r = load(istr);
    res = r.res;
    
    fprintf("10000: %s:\t%.2f (%d %d) \tT:%d \tMap:%d\n", Iname, res(2,6), res(5,6), res(6,6), res(4,6), ...
        res(7,6));
    if tt ~= 2
        fprintf("20000: %s:\t%.2f (%d %d) \tT:%d \tMap:%d\n", Iname, res(2,16), res(5,16), res(6,16), res(6,16), ...
            res(7,16));
        sumED = sumED  + res(2,16);
    end
    fprintf("----------------------------------------------------------\n");
end
sumED = sumED/7;
toc