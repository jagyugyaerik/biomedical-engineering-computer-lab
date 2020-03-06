function img=SPECT2DMLEMSM( projection, A, iter, img, bar )

if (nargin < 3)
        error('Too few arguments: proj = img=SPECT2DMLEMSM( <projection>, <System Matrix>, <number of iterations>, [<initial image>])');
end
if (nargin < 4)
	img = ones( [A.imagesize A.imagesize] );
end
if (nargin < 5)
    bar = true;
end
if (size(img) ~= [A.imagesize A.imagesize])
        error('Invalid image dimensions!');
end

if (bar)
    h=waitbar(0.0, 'SPECT2D System matrix based ML-EM operations');
end
warnzero=warning( 'off', 'MATLAB:divideByZero' );
for i=0:(iter-1)
	% forward projection
	fp = SPECT2DForwardProjSM( img, A );

	% compute the error
	e = projection ./ fp;
	e(find(~isfinite(e)))=0;

	r = SPECT2DBackwardProjSM( e, A );
%	for y=1:A.imagesize
%	for x=1:A.imagesize
%		img(x,y) = img(x,y) / sum(A.matrix(x+(y-1)*A.imagesize,:)) * r(x,y);
%	end
%	end
	img(:) = img(:) ./ sum(A.matrix,2) .* r(:);
    if (bar && ishandle(h))
        waitbar( 1.0/iter * (i+1), h);
    end
end
warning( warnzero );
if (bar && ishandle(h))
    close(h);
end

end
