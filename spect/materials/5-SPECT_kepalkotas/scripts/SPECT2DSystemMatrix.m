%
% SPECT2DSystemMatrix Creates a system matrix for SPECT
%
% A = SPECT2DSystemMatrix( imagesize, theta, [projsize] )
%
%    imagesize    is the size of the emission image
%    theta          is a vector of detector angles 
%    projsize    is the size of each projection.
%    
%    The parameter projsize is optional.
%
%   Note that this function creates the WHOLE system matrix, which consumes
%   a LOT of memory. (For a 128x128 image with 128 projection the matrix is
%   approximatley 3GB.)
%
%   Example
%   -------
%
%    A = SPECT2DSystemMatrix( 64, 0:6:359 );
%    img = readRawR32( 'phantom-ring-128x128.raw', [128 128] );
%    p = SPECT2DForwardProjSM( img, A );
%    showImage(p);
%

function A = SPECT2DSystemMatrix( imagesize, theta, projsize )

if (nargin < 2)
    error('invalid number of parameters!');
end
if (nargin < 3 || projsize <= 0)
    d=ceil(imagesize*sqrt(2));
    projsize = d + mod(d,2) + 3;
end

A = {};
A.imagesize = imagesize;
A.theta = theta;
A.projsize =  projsize;
A.matrix = zeros( [imagesize*imagesize, size(theta,2)*projsize] );

h = waitbar( 0.0, 'SPECT 2D ideal System Matrix generator');
center=imagesize/2+0.5;
detcenter=projsize/2+0.5;
for y=1:imagesize
for x=1:imagesize
    img_o = x+(y-1)*imagesize;
    for t=1:size(theta,2)
        % project pixel x, y in direction theta(t)
        r = ([x,y]-[center,center])*[-sin(theta(t)/180*pi), cos(theta(t)/180*pi)]' + detcenter;
        r1 = floor(r);
        r2 = r1+1;
        s1 = r-r1;

        if (r1 > 0 && r1 <= projsize)
            A.matrix(img_o, (t-1)*projsize+r1) = 1-s1;
        end
        if (r2 > 0 && r2 <= projsize)
            A.matrix(img_o, (t-1)*projsize+r2) = s1;
        end
    end
end
waitbar(1.0/(imagesize)*y, h);
end
close(h);

end
