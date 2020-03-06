% SPECT2DRRSystemMatrix Creates a system matrix for SPECT
%
%       S = SPECT2DRRSystemMatrix( imagesize, theta, projsize, psf_a, psf_b, sigma_intr, ror, step )
%
%   Where
%       imagesize       is the size of the emission image
%       theta           is a vector of detector angles 
%       projsize        is the size of each projection.
%       psf_a, psf_b	are the parameters of the collimator's depth 
%                   	dependent point-spread function in mm
%       sigma_intr      is the crystal's intrinsic blur in mm
%       ror             is the detector's radius of rotation in mm
%       step            is the size of each voxel in mm
%
%       If the parameter projsize is 0, then the default value is computed.
%
%   Note that this function creates the WHOLE system matrix, which consumes
%   a LOT of memory. (For a 128x128 image with 128 projection the matrix is
%   approximatley 3GB.)
%
%   Example
%   -------
%
%       S = SPECT2DRRSystemMatrix( 64, 0:3:179, 0, 0.20, 0.03, 2.0, 235, 6 );
%       img = readRawR32( 'phantom-ring-128x128.raw', [128 128] );
%       p = SPECT2DForwardProjSM( img, S );
%       showImage(p);
%
function A = SPECT2DRRSystemMatrix( imagesize, theta, projsize, psf_a, psf_b, sigma_intr, ror, step )

if (nargin < 8)
	fprintf([ ...
'Use: SPECT2DRRSystemMatrix( imagesize, theta, projsize, psf_a, psf_b, sigma_intr, ror, step )\n' ...
'	Where \n' ...
'	imagesize: dimension of the image\n' ...
'	theta: list of projection angles\n' ...
'	projsize: dimension of each projection\n' ...
'	psf_a, psf_b: collimators depth dependent point-spread function in mm\n' ...
'	sigma_intr: crystals intrinsic blur in mm\n' ...
'	ror: spect head rotation radius in mm\n' ...
'	step: size of each voxel in mm\n' ...
'\n' ...
]);
	error('Bad number of arguments');
end

if ( projsize <= 0 )
        d=ceil(imagesize*sqrt(2));
        projsize = d + mod(d,2) + 3;
end

A = {};
A.imagesize = imagesize;
A.theta = theta;
A.projsize =  projsize;
A.psf_a = psf_a;
A.psf_b = psf_b;
A.sigma_intr = sigma_intr;
A.ror = ror;
A.step = step;
A.matrix = zeros( [imagesize*imagesize, size(theta,2)*projsize] );

h = waitbar( 0.0, 'SPECT 2D System Matrix generator width depth-dependent blurring');
center=imagesize/2+0.5;
detcenter=projsize/2+0.5;

scale=8;

for y=1:imagesize
    for x=1:imagesize
        img_o = x+(y-1)*imagesize;
        for t=1:size(theta,2)
            % project pixel x, y in direction theta(t)
            r = ([x,y]-[center,center])*[-sin(theta(t)/180*pi), cos(theta(t)/180*pi)]' + detcenter;
            % distance from the collimator
            d = ror + ([x,y]-[center,center])*step*[cos(theta(t)/180*pi), sin(theta(t)/180*pi)]';
            sigma = sqrt( sigma_intr^2 + (psf_a + d*psf_b)^2 );

            if (sigma > 0)
                lr = max(1, floor(r-3*sigma/step));
                hr = min(projsize, ceil( r+3*sigma/step));
                hp = step/scale/(sqrt(2*pi)*sigma) * exp(-(((lr:(1/scale):hr) - r)*step).^2/2/sigma^2);
                for p=lr:(hr-1)
                    A.matrix(img_o,(t-1)*projsize+p) = sum(hp( ((p-lr)*scale+1):((p-lr+1)*scale) ));
                end
            else
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
    end
    if (ishandle(h))
        waitbar(img_o/(imagesize^2), h);
    end
end
close(h);
end
