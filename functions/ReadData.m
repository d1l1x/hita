function [uvel,vvel,wvel,time,dim] = ReadData(datadir,...
                                              flag,...
                                              u_name,...
                                              v_name,...
                                              w_name)
    tic; % enable timer
    uvel=importdata([datadir,'/',u_name]);
    vvel=importdata([datadir,'/',v_name]);
    wvel=importdata([datadir,'/',w_name]);
    time = toc; % end timer
    if strcmp(flag,'3D')
        dim=round((size(uvel,1))^(1/3));
    end
end
