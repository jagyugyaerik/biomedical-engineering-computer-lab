function img=SPECT2DMLEM( projection, S, mumap, iter, img, bar )

% SPECT2DAttenuatedMLEM   Implements the ML-EM algorithm.
%
%       img=SPECT2DAttenuatedMLEM( projection, S, mumap, iter, [img], [bar] )
%
% Where
%       projection      is a series of projection images
%       S               is a System descriptor (see SPECT2DSystem)
%	mumap		is an image sized linear attenuation coefficient map
%       iter            is the number of iterations
%       img             is the initial image (optional)
%       bar             set this to false to disable the progress bar
%  
% Example
%       
%       source = readRawR32('phantom-ring-128x128.raw', [128 128]);
%       mumap = readRawR32('attenuation-cyl-128x128.raw', [128 128]);
%       S = SPECT2DSystem(128, 0:3:359, 128, 0.97030, 0.017239, 2, 280, 3 );
%       proj = SPECT2DAttenuatedForwardProj( source, S, mumap );
%       proj = SPECT2DAddNoise( proj );
%       r = SPECT2DAttenuatedMLEM( proj, S, mumap, 10 );
%       showImage(r, 'MLEM reconstruction (10 iterations) of noisy attenuated projection series');

if (nargin < 6)
	bar = 1;
end
if (nargin < 5)
	img = ones( [S.imagesize S.imagesize] );
end
if (nargin < 4)
	error('Too few arguments. <reconstructed image> = SPECT2DMLEM( <projection>, <system descriptor>, <number of iterations>, [initial image] )');
end
if (size(projection) ~= [S.projsize, size(S.theta,2)])
	error('Projection image is not compatible with the defined SPECT system.');
end
if (size(img) ~= size(mumap))
        error('Image has to be the same dimensions as the mumap!');
end

if (bar)
h=waitbar(0.0, 'SPECT2D ML-EM operations');
end
norm = SPECT2DAttenuatedBackwardProj( ones([S.projsize, size(S.theta)]), S, mumap );
for i=0:(iter-1)
	% forward projection
	fp = SPECT2DAttenuatedForwardProj( img, S, mumap );

	% compute the error
	e = projection ./ fp;
	e(find(~isfinite(e)))=0;

	% backproject on norm
	bp = SPECT2DAttenuatedBackwardProj( e, S, mumap ) ;
	img = img .* bp ./ norm;

	if (bar && ishandle(h))
	waitbar( 1.0/iter * (i+1), h);
end
end
if (bar && ishandle(h))
close(h);
end

end	%function
