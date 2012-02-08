function [uvel,vvel,wvel,time] = ReadData(datadir,flag)
    tic; % enable timer
    uvel=importdata([datadir,'/',flag,'/uvel']);
    vvel=importdata([datadir,'/',flag,'/vvel']);
    wvel=importdata([datadir,'/',flag,'/wvel']);
    time = toc; % end timer
end
