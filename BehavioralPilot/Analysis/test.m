prompt = {'Debug mode?','Subject Initials','Subject Number','Session (1 or 2)',...
             'Part (1-6)', 'Run Number (1 or 2)', 'Training Run? (0 or 1)','Image Set'};
dlgtitle = 'Enter Run Parameters';
dims = [1 35];
definput = {'0','XX','99','1','1','1','0','3'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
p.debug = str2double(answer{1});
p.Subject = answer{2};
p.SubNum = sprintf('%02d',str2double(answer{3}));
p.Session = str2double(answer{4});
p.Part = str2double(answer{5});
p.RunNumThisPart = str2double(answer{6});
p.Training = str2double(answer{7});
p.ImageSet = str2double(answer{8});

p.SubNum
if ~mod(p.SubNum,2)
    fprintf('even, swap\n')
else
    fprintf('odd, no swap\n')
end