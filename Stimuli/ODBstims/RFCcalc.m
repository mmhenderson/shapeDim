function X = RFCcalc(par,mysize,stimtype,pointslist)
% modified from Hans Op de Beeck 
% 2006-9 ddrucker@psych.upenn.edu
% copied from github.com/dmb/thesis by MMH 2019
% 07/2019 modified by MMH (mmhender@ucsd.edu)
%% 
w = par.w; 
A = par.A;
P = par.P;
step = par.step;

assert(mod(mysize,2)==0);
center = mysize/2;

% all the original parameters here are set to work with size [299x299] so
% we need to adjust based on the size we actually want here (higher
% resolution!)
size_adjustment=mysize/299; 
A=A.*size_adjustment;
step=step.*size_adjustment;

% this is an adjustment that seems to be arbitrary
new_step=5*step./(sqrt(16)-1); 

% this is the approx radius of the shape, in pixels
r=round(70*size_adjustment);
% this is the approx circumference in pixels
npixel=round(2*pi*r);


for pointnum=1:size(pointslist,1)
    dim1val = pointslist(pointnum,1);
    dim2val = pointslist(pointnum,2);

    % how much are we stepping in amplitude? This is the coordinate you
    % entered, scaled by the "step" which is set in ODBParams.m
    % note that we're always changing only dims 2 and 5, others are fixed
    amp_step = new_step.*[0,dim1val-1, 0, 0, dim2val-1, 0,0];
   
    % now define the relationship between angle and radius along the curve.
    % this is equal to whatever radius we've set up above (r) plus the
    % modulation of the curve at this point (amp*sin(freq*theta+phase), 
    % summed over all freq components)
    curve_func = @(theta) r + sum((A+amp_step).*sin(w.*theta'+P),2);
    
    % a hacky solution to correct for the fact that we have non-integer
    % frequencies here...the start and stop points will initially not line
    % up perfectly. 
    % first define the start point
    yb = curve_func(0);
    % then define the end point
    ye = curve_func(2*pi);
    % this is how much they'll be off by...
    d_1e = ye-yb;

    % initialize the coordinate list. 
    Xpoly = zeros(npixel+1,1);
    Ypoly = zeros(npixel+1,1);

    %% looping over pixels in my curve here (npixels+1 total pts)
    for ii = 0:npixel
        
        % where are we along the curve, in radians?
        % (when ii==npixels we're back to the start)
        theta = 2*pi/npixel*ii;
        rho = curve_func(theta);
        
        % what are the x and y coordinates of this point? Polar to
        % cartesian conversions. Except we actually want to flip the x and 
        % y coordinates and make x negative, just to match the original code.
        yval = rho*cos(theta);
        xval = (-1)*rho*sin(theta);
       
        % apply a correction to the y coordinate based on what we did above
        % (correct for the difference between start and stop pts)
        yval = yval + cos(theta/2)*d_1e/2;
        
        % finally add the center on and round to get pixel space coordinates
        xcoord = round(xval+center);
        ycoord = round(yval+center);
            
        Xpoly(ii+1) = xcoord;
        Ypoly(ii+1) = ycoord;
         
    end

     %% Smooth the curve a bit  
    windsize = floor(5.84*size_adjustment);
    wind = ones(1,windsize);
    Xpoly = filtfilt(wind/windsize,1,Xpoly);
    Ypoly = filtfilt(wind/windsize,1,Ypoly);
  
    %% Make the final image (filled contour)
    Im = double(roipoly(zeros(mysize),Xpoly,Ypoly));
    
    %% Get rid of the holes!    
    Im = imfill(Im, 'holes');

    Im=Im';
    X{pointnum}=Im;

end
