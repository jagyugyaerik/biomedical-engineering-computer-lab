function p=SPECT2DForwardProj( img, S, theta)

% SPECT2DForwardProj   Calculates the forward projection series of an image.
%
%       p=SPECT2DForwardProj( image, S, theta )
%
% Where 
%       image  		is the source distribution to be projected
%       S               is a SPECT2D System descriptor (see SPECT2DSystem)
%       theta           is a vector of projection IDs from the descriptor
%
% Example
%       p = SPECT2DForwardProj( img, S )
%               Creates the full forward projection of an image.
%
%       p = SPECT2DBackwardProj( img, S, 1 )
%               Creates only one projection image (the first one)
%

if (nargin < 3) 
	theta=1:size(S.theta,2);
end;

p=zeros( S.projsize, size(theta,2));

imgsize= size(img);
offset = ceil((sqrt(2)-1)*max(imgsize)/2);
newsize= imgsize+2*offset;

if ( max(newsize)*S.step > 2*S.ror )
        error('Invalid SPECT System: the whole image must be in the FOV!');
end

for t=1:size(theta,2)
	bigimg = zeros( newsize );
	bigimg( (1+offset):(imgsize(1)+offset), (1+offset):(imgsize(2)+offset) ) = img;
	
	rot = myimrotate( bigimg, -S.theta(theta(t)), 'bilinear', 'crop' );
%	rot = rot * sum(img(:))/sum(rot(:));
	if (S.psf_a || S.psf_b)
		pc = size(rot,2);
		for d=1:size(rot,1)
%			dist = S.ror + ((size(rot,1)-d+0.5) - size(rot,1)/2)*S.step;
			dist = S.ror + ((d+0.5) - size(rot,1)/2)*S.step;
			sigma = S.psf_a + S.psf_b * dist;
			gauss = 1/(sqrt(2*pi)*sigma) * exp( -(((1:pc)-(pc/2)) * S.step ).^2 / 2 /sigma^2 );
			gauss = gauss(find(gauss>1e-6));
			gc = size(gauss,2);
			sp = conv( rot(d,:), gauss);
			rot(d,:) = sp( (floor(gc/2)+1):(pc+floor(gc/2)) );
		end
	end
	proj = sum(rot);
	if (S.sigma_intr)
		pc = size(proj,2);
		gauss = 1/(sqrt(2*pi)*S.sigma_intr) * exp( -(((1:pc)-floor(pc/2)) * S.step ).^2 / 2 /S.sigma_intr^2 );
		gauss = gauss(find(gauss>1e-6));
		gc = size(gauss,2);
		sp = conv( proj, gauss);
		proj = sp( (floor(gc/2)+1):(pc+floor(gc/2)) );
	end
	pc = floor(size(proj,2)/2);
	ac = floor(S.projsize / 2);
	if (size(proj,2) > S.projsize)
		p( :, t)=proj( (1+offset):(offset+S.projsize) );
	else
		p( (ac-pc+1):(ac+pc+mod(size(proj,2),2)), t)=proj( : );
	end
end

end
