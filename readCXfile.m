function [t,a] = readCXfile(cxFile)
%READCXFILE Summary of this function goes here
%   Detailed explanation goes here
    if isfile(cxFile)
        tb = readtable(cxFile);
        t = posixtime(datetime(string(tb(:,2).Date) + " " + string(tb(:,3).Time),'InputFormat',"MM/dd/uuuu HH:mm:ss.SSSSSS"));
        a = table2array(tb(:,4:6));
    else
        disp("Error: the input file cannot be opened.")
    end
end

