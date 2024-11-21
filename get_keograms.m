%used to download the keogram NCDF files from the website
%Make  sure there is a file that exists and is indicated in save_folder
%where the keogram NCDFs will be saved
function [] = get_keograms(KeoDates, Years_Running, root_dir) 
%% Setting Stuff

ftobj = ftp('optics.gi.alaska.edu');
addpath(genpath(root_dir));
% File name example PKR_SMSP_STD_20160101.NC
% ftp://optics.gi.alaska.edu/PKR/DMSP/NCDF/2016/ URL
dir_beg = 'PKR/DMSP/NCDF/';

%% Downloading Keograms
% Comment in to download all of the keograms
for i = 1:length(Years_Running)
    save_folder = fullfile(root_dir, 'Data', num2str(Years_Running(i)));
    if ~exist(save_folder, 'dir')
        mkdir(save_folder)
    end
    cd(save_folder);
    ftobj = ftp('optics.gi.alaska.edu');
    cd(ftobj, ['PKR/DMSP/NCDF/' num2str(Years_Running(i))]);
    for k = 1:length(KeoDates{i})
        try
            filename = ['PKR_SMSP_STD_' num2str(KeoDates{i}(k)) '.NC'];
        catch
            disp('pause');
        end
        if isfile(filename)
            disp('Already Downloaded');
        else
            try
                  mget(ftobj, filename, save_folder);
                disp('File Downloaded');
            catch
                disp('File does not exist');
            end
            
        end
    end
    cd([root_dir]);
end

disp('All Done!');
end