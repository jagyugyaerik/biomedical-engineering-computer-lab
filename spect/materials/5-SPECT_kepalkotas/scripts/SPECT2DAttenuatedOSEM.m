function img=SPECT2DAttenuatedOSEM( projection, S, mumap, iter, subsets, img, bar )

% SPECT2DAttenuatedOSEM   Implements the OS-EM algorithm.
%
%       img=SPECT2DAttenuatedOSEM( projection, S, mumap, iter, subsets, [img], [bar]
% )
%
% Where
%       projection      is a series of projection images
%       S               is a System descriptor (see SPECT2DSystem)
%       mumap           is an image sized linear attenuation coefficient map
%       iter            is the number of iterations
%	subsets		is the number of subsets (1 for ML-EM algorithm)
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
%       r = SPECT2DAttenuatedOSEM( proj, S, mumap, 10, 4 );
%       showImage(r, 'OSEM reconstruction (10 iterations, 4 subsets) of noisy attenuated projection series');

if (nargin < 7)
	bar = true;
end
if (nargin < 6)
	img = ones( [S.imagesize S.imagesize] );
end
if (size(img) == [1 1])
   	img = ones( [S.imagesize S.imagesize] );
end
if (nargin < 5)
	error('Too few arguments. <reconstructed image> = SPECT2DAttenuatedOSEM( <projection>, <system descriptor>, <attenuation map>, <number of iterations>, <number of subsets>, [initial image] )');
end
if (size(projection) ~= [S.projsize, size(S.theta,2)])
	error('Projection image is not compatible with the defined SPECT system.');
end
if (size(img) ~= size(mumap))
        error('Image has to be the same dimensions as the mumap!');
end

if (bar)
h=waitbar(0.0, 'SPECT2D OS-EM operations');
end

norm = {};
sorder=SPECTSubsetOrder(subsets);
for s=1:subsets
	st = sorder(s):subsets:size(S.theta,2);
	norm{s} = SPECT2DAttenuatedBackwardProj( ones([S.projsize, size(st)]), S, mumap, st );
end

for i=0:(iter-1)
    for s=1:subsets
        st = sorder(s):subsets:size(S.theta,2);

        % forward projection
        fp = SPECT2DAttenuatedForwardProj( img, S, mumap, st );

        % compute the error
%        e = projection(:,st) ./ fp;
%        e(find(~isfinite(e)))=0;
        e = zeros(size(fp));
        for ps = 1:size(fp,2)
            for p=1:size(fp,1)
                if (fp(p,ps) ~= 0)
                    e(p, ps) = projection(p, st(ps)) / fp(p, ps);
                end
            end
        end

        % backproject on norm
        bp = SPECT2DAttenuatedBackwardProj( e, S, mumap, st ) ;
        img = img .* bp ./ norm{s};
    end
    if (bar)
    	waitbar( 1.0/iter * (i+1), h);
    end
end
if (bar)
    close(h);
end

end
