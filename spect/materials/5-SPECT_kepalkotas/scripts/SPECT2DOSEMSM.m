function img=SPECT2DOSEMSM( projection, A, iter, subsets, img, bar )

% SPECT2DOSEMSM   Implements the OS-EM algorithm.
%
%       img=SPECT2DOSEMSM( projection, A, iter, subsets, [img], [bar] )
%
% Where
%       projection      is a series of projection images
%       A               is a System matrix (see SPECT2DSystemMatrix, SPECT2DRRSystemMatrix)
%       iter            is the number of iterations
%	subsets             is the number of subsets (1 for ML-EM algorithm)
%       img             is the initial image (optional)
%       bar             set this to false to disable the progress bar
%  
% Example
%       
%       source = readRawR32('phantom-ring-128x128.raw', [128 128]);
%       A = SPECT2DRRSystemMatrix(128, 0:3:359, 128, 0.97030, 0.017239, 2, 280, 3 );
%       proj = SPECT2DForwardProjSM( source, A );
%       proj = SPECT2DAddNoise( proj, 100000 );
%       r = SPECT2DOSEMSM( proj, A, 10, 4 );
%       showImage(r, 'OSEM reconstruction (10 iterations, 4 subsets) of noisy projection series');

if (nargin < 6)
	bar = true;
end
if (nargin < 5)
	img = ones( [A.imagesize A.imagesize] );
end
if (length(img) == 1)
   	img = ones( [A.imagesize A.imagesize] );
end
if (nargin < 4)
	error('Too few arguments. <reconstructed image> = SPECT2DOSEMSM( <projection>, <system matrix>, <number of iterations>, <number of subsets>, [initial image] )');
end
if (size(projection) ~= [A.projsize, length(A.theta)])
	error('Projection image is not compatible with the defined SPECT system.');
end

if (bar)
    h=waitbar(0.0, 'SPECT2D OS-EM operations');
end

norm = {};
sorder=SPECTSubsetOrder(subsets);
for s=1:subsets
	st = sorder(s):subsets:size(A.theta,2);
	norm{s} = SPECT2DBackwardProjSM( ones([A.projsize, length(st)]), A, st );
end

for i=0:(iter-1)
    for s=1:subsets
        st = sorder(s):subsets:size(A.theta,2);

        % forward projection
        fp = SPECT2DForwardProjSM( img, A, st );

        % compute the error
        e = projection(:,st) ./ fp;
        e(find(~isfinite(e)))=0;

        % backproject on norm
        bp = SPECT2DBackwardProjSM( e, A, st ) ;
        img = img .* bp ./ norm{s};
    end
    if (bar && ishandle(h))
    	waitbar( 1.0/iter * (i+1), h);
    end
end
if (bar && ishandle(h))
    close(h);
end

end
