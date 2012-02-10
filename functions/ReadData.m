function [uvel,vvel,wvel,time] = ReadData(datadir,flag,u_name,v_name,w_name)
    tic; % enable timer
    uvel=importdata([datadir,'/',flag,'/',u_name]);
    vvel=importdata([datadir,'/',flag,'/',v_name]);
    wvel=importdata([datadir,'/',flag,'/',w_name]);
    time = toc; % end timer
end
