function X = RFCcalc(par,mysize,stimtype,pointslist)
% modified from Hans Op de Beeck 
% 2006-9 ddrucker@psych.upenn.edu
% 07/2019 modified by MMH (mmhender@ucsd.edu)
%% 
w = par.w; 
A = par.A;
P = par.P;
step = par.step;

if strcmp(stimtype,'grid')
    pointslist=points(0,4)+2.5;
end
assert(mod(mysize,2)==0);
center = mysize/2;

size_adjustment=mysize/299; %original parameters calculated for larger stimulus size (299x299 instead of 256x256 as generated now)
A=A.*size_adjustment;
step=step.*size_adjustment;
%new_step=5*step./(sqrt(length(pointslist))-1); % we REALLY want length(pointslist) to be 16 though
new_step=5*step./(sqrt(16)-1); % we REALLY want length(pointslist) to be 16 though

%names=char(65:65+length(pointslist)-1);


% smallest distance between adjacent stimuli becomes (5*step)/(Ndim-1)
%     the original step was calculated for 6 values/dimension instead of 3
% step is not necessarily equal for each manipulated dimension (maximum = 6)
%		because it was originally calibrated to generate
%     equal stimulus differences in pixel space

%number of pixels on circle outline
%amplitude expressed in radials is a function of the ratio of A to r

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

    % initialize the coordinate list. Note this is more points than the
    % npixels, we'll use interpolation so we don't know exactly how many
    % points until the end.
    Xpoly = [];Ypoly = [];

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

        if ii>0
            % if after the first point, see if we need to interpolate.
            % this happens if the modulation rate is high. 
            xrange = [Xpoly(end),xcoord];
            yrange = [Ypoly(end),ycoord];
          
            npts_interp = max(abs([diff(xrange), diff(yrange)]))+1;
            
            x_interp = linspace(xrange(1),xrange(2),npts_interp)';
            y_interp = linspace(yrange(1),yrange(2),npts_interp)';
            
            % take off the first element here becuse it's already in the
            % list.
            x_interp = round(x_interp(2:end));
            y_interp = round(y_interp(2:end));

%             for p=1:numel(x_interp)
            Xpoly = [Xpoly; x_interp];
            Ypoly = [Ypoly; y_interp];
%                 [Xpoly(end-3:end),Ypoly(end-3:end)]
%             end
        else
            Xpoly =[Xpoly, xcoord];
            Ypoly =[Ypoly, ycoord];
        end

    end

     %% Smooth the curve a bit  
    windsize = floor(5.84*size_adjustment);
    wind = ones(1,windsize);
    Xpoly = filtfilt(wind/windsize,1,Xpoly);
    Ypoly = filtfilt(wind/windsize,1,Ypoly);
  
%     Xpoly = round(Xpoly, 1);
%     Ypoly = round(Ypoly, 1);
     %% now look for intersections.
    pgon = polyshape([Xpoly,Ypoly]);
    holes = pgon.holes;
    
    maxdist = 1;
    
    pts2remove = [];
    for hh=1:length(holes)
        
        vertices = holes(hh).Vertices;
        all_vertices = interppolygon(holes(hh).Vertices,10000);
        
        all_dist = squeeze(sqrt(sum((repmat([Xpoly, Ypoly], 1, 1, 10000) - permute(repmat(all_vertices, 1,1,numel(Xpoly)),[3,2,1])).^2,2)));
        pts = find(any(all_dist<maxdist, 2));
%         pts = find(ismember([Xpoly,Ypoly],all_vertices,'rows'));
        %        pts = find(ismember([Xpoly,Ypoly],pgon,'rows'));
        pts2remove = [pts2remove;pts];
    end
    
    % can plot these areas to verify that they're the correct portion to be
    % removing...
    figure;hold all;
    plot(Xpoly,Ypoly,'.');
    for hh=1:length(holes)
           
        vertices = holes(hh).Vertices;
        all_vertices = interppolygon(holes(hh).Vertices,1000);
        
        plot(Xpoly(pts2remove),Ypoly(pts2remove),'ro');
    end
    
    %     %% now look for intersections.
%     % note there is probably a better way to do this! 
%     maxdist = 0.0001;
%     minloopsize = 10;
%     inner_loop = 0;
%     pts2remove = zeros(length(Xpoly),1);
%     
%     % looping over all points except the first and last because they're
%     % equal
%     for pp=2:length(Xpoly)-1
% 
%         this_coord = [Xpoly(pp),Ypoly(pp)];   
%         
%         % see if this point is intersecting any other points.
%         int_points = find(sqrt(sum((this_coord-[Xpoly, Ypoly]).^2,2))<maxdist);
%         
%         % doesn't count if it intersects with itself 
%         int_points(int_points==pp|int_points==pp+npixel) = [];
%         
%         % also doesn't count if the intersection is with a point very
%         % nearby in time
%         int_points(abs(int_points-pp)<minloopsize) = [];
%         
%         %finally, doesn't count if we have found an intersection very
%         %recently
%         last_int = find(diff([0;pts2remove(1:pp-1)]),1,'last');
%         
%         if ~isempty(int_points) 
%             if (isempty(last_int) || pp-last_int>=minloopsize)
%                 inner_loop = 1-inner_loop;
%             end
%         end
%         
%         pts2remove(pp) = inner_loop;
%     end
%     
%     pts2remove = boolean(pts2remove);
% 
%     % can plot these areas to verify that they're the correct portion to be
%     % removing...
%     figure;hold all;
%     plot(Xpoly,Ypoly,'.');
%     plot(Xpoly(pts2remove),Ypoly(pts2remove),'ro');
    
    % take them out!
%     Xpoly(pts2remove) = [];
%     Ypoly(pts2remove) = [];
    
     %% Smooth the curve a bit  
% %     windsize = floor(5.84*size_adjustment);
%     windsize = floor(1*size_adjustment);
%     wind = ones(1,windsize);
%     Xpoly = filtfilt(wind/windsize,1,Xpoly);
%     Ypoly = filtfilt(wind/windsize,1,Ypoly);
    %% Make the final image (filled contour)
    Im = roipoly(zeros(mysize),Xpoly,Ypoly);
    Im=Im';

    X{pointnum}=Im;

end
