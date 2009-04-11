function outputStuff = mfcsvread(fileName)
%mfcsvread reads a CSV file containing both text & numeric data.  MATLAB's
%csvread function will work with all numeric, or all text data, but not
%both.  It's common to have a file with a single line of comma separated
%text headers followed by many rows of numeric data.  xlsread is limited in
%the number of rows & colums (actually, Excel is the limitation) it can
%read.
%
% The CSV file should look like:
% comma, separated, text, ...
% 1,2,3,4,5,...
% 6,7,8,9,10,...
% etc...
%
% The output is a structure with the column headers as fields in a
% structure each with a vector of data.

% Open the file
fid=fopen(fileName);
 if fid == -1 
     disp('Cannot open file');
     return;
 end

% Start reading the data
tline = fgetl(fid); % read and discard first 
tline = fgetl(fid); % and second lines
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid); % Read in the first line only (text headers)
tline = tline(tline~=' '); %Get rid of any spaces
commaLocs=findstr(',',tline); % find the commas
fieldNames=cell(1,(length(commaLocs)+1));
start=1;
for colIdx=1:length(commaLocs)
    fieldNames{colIdx}=tline(start:commaLocs(colIdx)-1);
    start=commaLocs(colIdx)+1;
end
fieldNames{colIdx+1}=tline(start:end);

%Read in the rest of the data (should be numeric)
outputStuff=csvread(fileName,5,0);

%Convert the data into a single structure.
%[cellW cellH]=size(fieldNames);
%numFields=max([cellW cellH]);
%Needs to be a 1xN or Nx1 cell-array.
%if ~iscell(fieldNames)
%    disp('Needs to be a cell-array');
%    return;
%elseif cellW ~=numFields && cellH ~=numFields
%    disp('Needs to be 1xN or Nx1');
%    return;
%end

%dataSize=size(fieldData);
%arrayDim=find(dataSize==numFields);
%if isempty(arrayDim)
%    disp('Dimensions are wrong');
%end

%outputStuff=[];
%if arrayDim==2
%    for idx=1:numFields
%        outputStuff.(fieldNames{idx})=fieldData(:,idx);
%    end
%else
%    for idx=1:numFields
%        outputStuff.(fieldNames{idx})=fieldData(idx,:);
%    end
%end
