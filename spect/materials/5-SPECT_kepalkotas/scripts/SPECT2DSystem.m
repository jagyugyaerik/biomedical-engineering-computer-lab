
% SPECT2DSystem Creates a system descriptor for SPECT
%
%       S = SPECT2DSystem( imagesize, theta, projsize, psf_a, psf_b, sigma_intr, ror, step )
%
%   Where
%       imagesize       is the size of the emission image
%       theta           is a vector of detector angles 
%       projsize        is the size of each projection.
%	psf_a, psf_b	are the parameters of the collimator's depth 
%		dependent point-spread function in mm
%       sigma_intr	is the crystal's intrinsic blur in mm
%	ror		is the detector's radius of rotation in mm
%	step		is the size of each voxel in mm
%
%       If the parameter projsize is 0, then the default value is computed.
%
%
%   Example
%   -------
%
%       S = SPECT2DSystem( 128, 0:3:179, 0, 0.20, 0.03, 2.0, 235, 1 );
%       img = readRawR32( 'phantom-ring-128x128.raw', [128 128] );
%       p = SPECT2DForwardProj( img, S );
%       showImage(p);
%

function S = SPECT2DSystem( imagesize, theta, projsize, psf_a, psf_b, sigma_intr, ror, step )

if (nargin < 8)
	error('Bad number of arguments: SPECT2DSystem( imagesize, theta, projsize, psf_a, psf_b, sigma_intr, ror, step)');
end
if ( projsize <= 0 )
        d=ceil(imagesize*sqrt(2));
        projsize = d + mod(d,2) + 3;
end


S = {};
S.imagesize = imagesize;
S.theta = theta;
S.projsize =  projsize;
S.psf_a = psf_a;
S.psf_b = psf_b;
S.sigma_intr = sigma_intr;
S.ror = ror;
S.step = step;

end
