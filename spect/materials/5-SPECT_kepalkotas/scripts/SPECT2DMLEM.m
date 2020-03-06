function img=SPECT2DMLEM( projection, S, iter, img, bar )

% SPECT2DMLEM	Implements the ML-EM algorithm.
%
%	img=SPECT2DMLEM( projection, S, iter, [img], [bar] )
%
% Where
%	projection	is a series of projection images
%	S		is a System descriptor (see SPECT2DSystem)
%	iter		is the number of iterations
%	img		is the initial image (optional)
%	bar		set this to false to disable the progress bar
%
% Example
%	
%	source = readRawR32('phantom-ring-128x128.raw', [128 128]);
%	mumap = readRawR32('attenuation-cyl-128x128.raw', [128 128]);
%	S = SPECT2DSystem(128, 0:3:359, 128, 0.97030, 0.017239, 2, 280, 3 );
%	proj = SPECT2DAttenuatedForwardProj( source, S, mumap );
%	proj = SPECT2DAddNoise( proj );
%	r = SPECT2DMLEM( proj, S, 10 );
%       showImage(r, 'MLEM reconstruction (10 iterations) of noisy attenuated projection series');


if (nargin < 5)
	bar = 1;
end
if (nargin < 4)
	img = ones( [S.imagesize S.imagesize] );
end
if (nargin < 3)
	error('Too few arguments. <reconstructed image> = SPECT2DMLEM( <projection>, <system descriptor>, <number of iterations>, [initial image] )');
end
if (size(projection) ~= [S.projsize, size(S.theta,2)])
	error('Projection image is not compatible with the defined SPECT system.');
end

if (bar)
	h=waitbar(0.0, 'SPECT2D ML-EM operations');
end
norm = SPECT2DBackwardProj( ones([S.projsize, size(S.theta)]), S );
for i=0:(iter-1)
	% forward projection
	fp = SPECT2DForwardProj( img, S );

	% compute the error
    warning off MATLAB:divideByZero;
	e = projection ./ fp;
	e(find(~isfinite(e)))=0;
    warning on MATLAB:divideByZero;

	% backproject on norm
	bp = SPECT2DBackwardProj( e, S ) ;
	img = img .* bp ./ norm;

	if (bar && ishandle(h))
		waitbar( 1.0/iter * (i+1), h);
	end
end
if (bar && ishandle(h))
	close(h);
end

end
