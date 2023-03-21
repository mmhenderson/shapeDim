function p = deg2pix(p)
% converts degrees visual angle to pixel units before rendering
% with PTB. Needs p.screenWidthCM and p.vDistCM
% js - 10.2007

% figure out pixels per degree, p.sRect(1) is x coord for upper left of
% screen, and p.sRect(3) is x coord for lower right corner of screen
if isfield(p,'screenWidthCM')
    p.ppd = pi * (p.sRect(3)-p.sRect(1)) / atan(p.screenWidthCM/p.vDistCM/2) / 360;
elseif isfield(p,'screenHeightCM')
    p.ppd = pi * (p.sRect(4)-p.sRect(2)) / atan(p.screenHeightCM/p.vDistCM/2) / 360;
else
    error('need either width or height of screen!')
end

% get name of each field in p
s = fieldnames(p);

% convert all fields with the word 'Deg' from degrees visual angle to
% pixels, then store in a renmaed field name
for i=1:length(s)
    ind = strfind(s{i}, 'Deg');
    if ind
        curVal = getfield(p,s{i});
        tmp = char(s{i});
        newfn = [tmp(1:ind-1), 'Pix'];
        p = setfield(p,newfn,curVal*p.ppd);
    end
end

end