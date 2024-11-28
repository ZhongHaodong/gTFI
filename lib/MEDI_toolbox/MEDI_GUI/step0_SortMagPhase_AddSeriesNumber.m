clc

if (1|[]), is_matlab = 1, else, is_matlab = 0, end
% This tests whether we are in matlab or octave
% I dont use this as I had trouble getting dicomread to work from octave
% but this is a simple way to test which is running due to a difference in matlab/octave syntax
undersample = 1;
max_n = 2000;  %ljq
filename = cell(1,max_n);
%path_source = '.';
%path_source = uigetdir
%path_source = 'E:\QSM\MarrowQSM\data_3T_20150912\monopolar_6echo'

%path_source = 'E:\QSM\data\LiCunboData'
path_source = 'F:\MRIdata\zhangmiao\phantom_2%_20171128\phantom_2%';

pathname_mag = [path_source '/' 'images_mag_all_echoes'];
pathname_phase = [path_source '/' 'images_phase_all_echoes'];
pathname_mag_phase = [path_source '/' 'images_mag_phase_all_echoes'];
                        
if(exist(pathname_mag) ~= 7)
           mkdir(pathname_mag)
 end;
 
 if(exist(pathname_phase) ~= 7)
           mkdir(pathname_phase)
 end;

 if(exist(pathname_mag_phase) ~= 7)
           mkdir(pathname_mag_phase)
 end;
 
dir_read = dir(fullfile(path_source,'0*.*'))

disp(['sorting the directoris.......']);

FirstSeriesNumber = 2;
LastSeriesNumber = 17;

for i=1:length(dir_read)
    
        if(dir_read(i).isdir == 0)
            continue;
        end
        
        i
        subpath_source = [path_source '/' dir_read(i).name];
        
        file_read = dir(fullfile(subpath_source,'*.dcm'));
        
        for j=1:length(file_read)
           
           file_name{j}=file_read(j).name;% 获取文件名的列表
           sourcename = [subpath_source '/' file_name{j}];
       
           meta =dicominfo(sourcename);
           
           if((meta.SeriesNumber >= FirstSeriesNumber) && (meta.SeriesNumber <= LastSeriesNumber))
                         
            if(strfind(meta.Private_0051_1016,'M')>0)
               destname = [pathname_mag '/' num2str(meta.SeriesNumber) file_name{j}];
            else
               destname = [pathname_phase '/' num2str(meta.SeriesNumber) file_name{j}];
            end
             copyfile(sourcename,destname,'f');
             
             destname = [pathname_mag_phase '/' num2str(meta.SeriesNumber) file_name{j}];
             copyfile(sourcename,destname,'f');
             
           end
                   
        end      
end

%clear all;

disp(['sort is done']);

return;


