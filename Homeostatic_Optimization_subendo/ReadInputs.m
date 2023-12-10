function outputs = ReadInputs(file)

strPath = 'Input';
strName = [file,'.txt'];
strFull = fullfile(strPath,strName);
fid = fopen(strFull);
iRow = 1;
while (~feof(fid)) 
    temp = textscan(fid,'%f\n','CommentStyle','%%');
    outputs(iRow, 1) = temp{1};
    iRow = iRow + 1;
end
outputs = [outputs];
fclose(fid);