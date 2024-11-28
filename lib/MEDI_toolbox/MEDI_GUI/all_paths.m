function filelist=all_paths()
DicomFolder='F:\FatWaterSeparation\zhai\fwtoolbox\fwToolBox\fwtoolbox_v1_code\gui';    
folder_list = get_all_pathes(DicomFolder,cell(0));    
    tempfilelist = struct([]);
    filelist = cell(0);
    for i=1:length(folder_list)
        tempfilelist = dir(folder_list{i});        
        for j =1:length(tempfilelist)           
            if tempfilelist(j).isdir == 0
                pathname =[folder_list{i},'\',tempfilelist(j).name];
                filelist = [filelist;pathname];
            end
        end
    end
end

function folder_list =  get_all_paths(path,folder_list)
    filelist = dir(path);
    i=1;
    tempFolder = cell(0);
    while(i<=length(filelist))
        if(filelist(i).isdir == 1)
            tempFolder = [tempFolder,filelist(i)];
        end
        i =i +1;
    end
    tempFolder = tempFolder(3:numel(tempFolder));
    numFolderPath = numel(tempFolder);
    tempFolderPath = cell(numFolderPath,1);
    for i = 1:numFolderPath
        tempFolderPath{i} = [path '\' tempFolder{i}.name];
    end
    folder_list = [folder_list;path];
    for i = 1 : numFolderPath
        folder_list = get_all_paths(tempFolderPath{i}, folder_list);
    end
end
