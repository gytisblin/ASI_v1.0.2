%Deleting unused iamges
load('DownloadList.mat');
arch_dir = 'C:\Users\Alex\Desktop\Old Images';
addpath(genpath(pwd));
addpath(genpath(arch_dir));
filesdownloaded = dir('C:\Users\Alex\Desktop\CloudShapeRun\ASCDownloads\Download\ASCDownloadFolder');
filesmoved = 0;
for i = 3:size(filesdownloaded, 1)
    i
    filename = filesdownloaded(i).name;
    ind = find(filename == ASCScintList);
    if isempty(ind)
        filesmoved = filesmoved+1;
        movefile(['C:\Users\Alex\Desktop\CloudShapeRun\ASCDownloads\Download\ASCDownloadFolder\' filename], arch_dir);
    end
end

disp([num2str(filesmoved) ' files were moved']);