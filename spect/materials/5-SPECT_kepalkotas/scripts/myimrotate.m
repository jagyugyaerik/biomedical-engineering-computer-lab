function tmp=myimrotate( img, theta, method, bbox )
% MYIMROTATE Rotate image.
%    B = IMROTATE(A,ANGLE) rotates image A by ANGLE degrees in a 
%    counterclockwise direction around its center point. To rotate the image
%    clockwise, specify a negative value for ANGLE. IMROTATE makes the output
%    image B large enough to contain the entire rotated image. IMROTATE uses
%    nearest neighbor interpolation, setting the values of pixels in B that 
%    are outside the rotated image to 0 (zero).
% 
%    B = IMROTATE(A,ANGLE,METHOD) rotates image A, using the interpolation
%    method specified by METHOD. METHOD is a string that can have one of the
%    following values. The default value is enclosed in braces ({}).
% 
%         {'nearest'}  Nearest neighbor interpolation
% 
%         'bilinear'   Bilinear interpolation
% 
%         'bicubic'    Bicubic interpolation. Note: This interpolation
%                      method can produce pixel values outside the original
%                      range.
% 
%    B = IMROTATE(A,ANGLE,METHOD,BBOX) rotates image A, where BBOX specifies 
%    the size of the output image B. BBOX is a text string that can have 
%    either of the following values. The default value is enclosed in braces
%    ({}).
% 
%         {'loose'}    Make output image B large enough to contain the
%                      entire rotated image. B is generally larger than A.
% 
%         'crop'       Make output image B the same size as the input image
%                      A, cropping the rotated image to fit. 
% 

if (nargin < 2)
	help myimrotate;
	error('Too few arguments');
end
if (nargin < 3)
	method = 'nearest';
end;
if (nargin < 4)
	bbox = 'loose';
end;

if (exist('imrotate') >= 2)
%	% if matlab has imrotate, then we should use that ... 
%	%	it's faster and more complete :)
	tmp = imrotate( img, theta, method, bbox );
%    tmp = tmp / sum(tmp(:)) * sum(img(:));
	return
end

if (size(size(img)) ~= [1 2])
	error('Image img must be a 2D matrix!');
end

if (strcmpi( bbox, 'crop'))
	tmp = zeros(size(img));
elseif (strcmpi( bbox, 'loose'))
	% it's really loose
	ms = max(size(img));
	newsize = ms + 2*ceil( (sqrt(2)-1)*ms/2 );
	tmp = zeros(newsize, newsize);
else
	error( ['BBox type [' bbox '] is not available, choose from crop or loose!']);
end

tmpcenter = size(tmp,1)/2;
imgcenter = size(img)/2;

if (strcmpi(method, 'nearest'))
elseif (strcmpi(method, 'bilinear'))
else
	warning([ method ' interpolation is not (yet) available, using bilinear instead.']);
	method='bilinear';
end


for y=1:size(tmp,2)
for x=1:size(tmp,1)
	ix =  (x-tmpcenter) * cos(theta/180*pi) + (y-tmpcenter) * sin(theta/180*pi) + imgcenter(1);
	iy = -(x-tmpcenter) * sin(theta/180*pi) + (y-tmpcenter) * cos(theta/180*pi) + imgcenter(2);
	val = 0;
	if (strcmpi(method, 'nearest'))
		ix = round(ix);
		iy = round(iy);
		if (ix>=1 && iy >=1 && ix <= size(img,1) && iy <= size(img,2))
			val = img(ix,iy);
		end
	elseif (strcmpi(method, 'bilinear'))
		px = floor(ix);
		py = floor(iy);
		if (py >= 1 && px >= 1  && px <= size(img,1)-1 && py <= size(img,2))
			%valx1 = (1 - rem(ix, 1))* img(px, py) + rem(ix,1) * img(px+1, py);
			valx1 = img(px:(px+1), py)' * [ 1-rem(ix,1), rem(ix,1) ]';
		else
			valx1 = 0;
		end
		if (py >= 0 && px >= 1  && px <= size(img,1)-1 && py <= size(img,2)-1)
			%valx2 = (1 - rem(ix, 1))* img(px, py+1) + rem(ix,1) * img(px+1, py+1);
			valx2 = img(px:(px+1), py+1)' * [ 1-rem(ix,1), rem(ix,1) ]';
		else
			valx2 = 0;
		end
		%val = (1 - rem(iy,1))*valx1 + rem(iy,1)*valx2;
		val = [valx1 valx2] * [ 1-rem(iy,1); rem(iy,1) ];
	end
	tmp(x,y) = val;
end
end

%tmp = tmp / sum(tmp(:)) * sum(img(:));
end
