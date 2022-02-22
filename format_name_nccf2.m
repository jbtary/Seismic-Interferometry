% Get common SAC data files for 2 OBS
% GOFAR experiment
% INPUTS:
% dir1,dir2: path to directories where the directories of OBS data are.
% obs1,obs2: name of directories of each OBS.
% 
% SAC names convention: 'G04_01_Feb_2008_H00.mat'
%                       name DD MMM YYYY  hh
% 
% OUPUTS
% info1,info2: list of data file names in common to both OBS
% 
% Examples
% dir1 = '/GOFAR/G04mat/';
% dir2 = '/GOFAR/G06mat/';
% 
% [info1,info2] = format_name_nccf(dir1,dir2);

function [info1,info2] = format_name_nccf(dir1,dir2)

inf1 = dir([dir1 '/*.mat']);
inf2 = dir([dir2 '/*.mat']);
name1 = extractfield(inf1,'name'); name1 = name1';
name2 = extractfield(inf2,'name'); name2 = name2';

dat1 = zeros(size(name1,1),1);
dat2 = zeros(size(name2,1),1);

for ii = 1:size(name1,1);
    
    day1 = str2num(name1{ii}(5:6));
    month1 = name1{ii}(8:10);
    yr1 = str2num(name1{ii}(12:15));
    hrm1 = str2num(name1{ii}(18:19));
    date = [num2str(day1) '-' month1 '-' num2str(yr1) ' ' num2str(hrm1) ':00:00'];
    dat1(ii,1) = datenum(date,'dd-mmm-yyyy HH:MM:SS');
    
    clear day1 month1 yr1 hrm1 date
end

for ii = 1:size(name2,1);
    
    day2 = str2num(name2{ii}(5:6));
    month2 = name2{ii}(8:10);
    yr2 = str2num(name2{ii}(12:15));
    hrm2 = str2num(name2{ii}(18:19));
    date = [num2str(day2) '-' month2 '-' num2str(yr2) ' ' num2str(hrm2) ':00:00'];
    dat2(ii,1) = datenum(date,'dd-mmm-yyyy HH:MM:SS');
    
    clear day2 month2 yr2 hrm2 date
end

[~,ia,ib] = intersect(dat1,dat2);
name1b = name1(ia);
name2b = name2(ib);

info1 = cell2struct(name1b,'name',2);
info2 = cell2struct(name2b,'name',2);
