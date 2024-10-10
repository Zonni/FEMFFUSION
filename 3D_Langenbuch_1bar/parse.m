function data = parse(filename, nlines, mode, headerline)
% usage:
% data = parse('filename', nlines, 'headerline');
% where for the input variables:
%     'filename'   = String with the name of the file
%     nlines       = Integer specifying how many lines to read
%     'headerline' = String with the header we will use to find the lines to read
% and the output:
%     data         = Matrix of size (nlines, lengthoflines) with numerical data.

% read the whole file to a temporary cell array
fid = fopen(filename,'rt');
tmp = textscan(fid,'%s','Delimiter','\r\n', 'CommentStyle', '%');
fclose(fid);

%% look for the headerline
tmp = tmp{1};

%# Ignore empty lines
tmp(cellfun(@(x)isempty(x), tmp)) = [];

%# Ignore small lines

%cellfun(@disp, tmp1)

% here we look for the header
idx = cellfun(@(x) strcmp(x,headerline), tmp);


ind_vec = find(idx);

ind=ind_vec(mode);

% idx is the line we are interesting in. We copy the lines
% from idx+1:idx+nlines to the data out

data = [];
for i = 1:nlines
    data =[data str2num(tmp{ind+i})];
end

% delete temporary array (if you want)
clear tmp

end
