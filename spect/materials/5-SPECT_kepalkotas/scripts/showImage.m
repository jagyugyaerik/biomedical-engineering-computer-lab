function showImage( img, caption, figureID, sub )

% showImage( img[, caption, [, figureID[, sub]]] )
%	Shows an image in the figure given by it's figureID.
%	You can select a subplot too with the sub argument.
%
% Example
%	figure(10);
%	subplot(2,2,1);
%	showImage( rand(10,10), 'Random image', 10, 221);
%	showImage( rand(10,10), 'Random image', 10, 222);
%	showImage( rand(10,10), 'Random image', 10, 223);
%	showImage( rand(10,10), 'Random image', 10, 224);
%


% check input arguments and report errors
if (nargin == 0 || nargin > 4)
    error('Use: showImage( img[, caption, [, figureID[, sub]]] );');
end

% change to the given figure and subplot
if (nargin >= 3)
    figure(figureID);
end
if (nargin == 4)
	if (size(sub)==[1 3])
		subplot( sub(1), sub(2), sub(3) );
	else 
		subplot( sub );
	end
end

imagesc(img);
colormap(gray);
colorbar;
axis image;

%if (exist('impixelinfo')>=2) 
%    impixelinfo; 
%elseif (exist('pixval')>=2)
%    pixval;
%end

if (nargin >= 2)
    title( caption );
end

end
