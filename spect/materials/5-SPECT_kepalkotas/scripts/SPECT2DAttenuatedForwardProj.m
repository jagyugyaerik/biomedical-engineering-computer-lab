function p=SPECT2DAttenuatedForwardProj( img, S, mumap, theta)

% SPECT2DAttenuatedForwardProj   Calculates the forward projection series of an image.
%
%       p=SPECT2DAttenuatedForwardProj( image, S, mumap, theta )  
%
% Where 
%       image           is the source distribution to be projected
%       S               is a SPECT2D System descriptor (see SPECT2DSystem)
%	mumap		is an image sized linear attenuation coefficient map
%       theta           is a vector of projection IDs from the descriptor
%
% Example
%       p = SPECT2DAttenuatedForwardProj( img, S, mumap )  
%               Creates the full forward projection of an image.
%
%       p = SPECT2DAttenuatedBackwardProj( img, S, mumap, 1 )
%               Creates only one projection image (the first one)
%

if (nargin < 3)
	error('Too few arguments: proj = SPECT2DAttenuatedForwardProj( image, System, mumap, [angle IDs])');
end
if (nargin < 4) 
	theta=1:size(S.theta,2);
end;
if (max(theta) > size(S.theta,2))
	error('Invalid angle ID');
end
if (size(img) ~= size(mumap))
	error('Image has to be the same dimensions as the mumap!');
end

imgsize = size(img);
offset = ceil((sqrt(2)-1)*max(imgsize)/2); 
newsize = imgsize+2*offset;

if ( max(newsize)*S.step > 2*S.ror )
	error('Invalid SPECT System: the whole image must be in the FOV!');
end


p=zeros( S.projsize, size(theta,2));

	bigimg  = zeros( newsize );
	bigimg( (1+offset):(imgsize(1)+offset), (1+offset):(imgsize(2)+offset) ) = img;

	bigmumap  = zeros( newsize );
	bigmumap( (1+offset):(imgsize(1)+offset), (1+offset):(imgsize(2)+offset) ) = mumap;
	

for t=1:size(theta,2)
	rot = myimrotate( bigimg, -S.theta(theta(t)), 'bilinear', 'crop' );
%    rot(find (rot < 0)) = 0;
%?	rot = rot * sum(img(:))/sum(rot(:));

	murot = myimrotate( bigmumap, -S.theta(theta(t)), 'bilinear', 'crop' );
%    murot(find(murot < 0)) = 0; % make sure all values are above 0 ... bicubic interpolation suxx on edges
%?	murot = murot * sum(mumap(:))/sum(murot(:));

	% compute the integral mu-map
	for d=2:size(murot,1)
		murot(d,:)=murot(d,:)+murot(d-1,:);
	end
	% add the attenuation-effect to the rotated map
	for d=1:size(murot,1)
		rot(d,:)=rot(d,:).*exp(-S.step*murot(d,:)/10);
	end

	if (S.psf_a || S.psf_b)
		pc = size(rot,2);
		for d=1:size(rot,1)
%			dist = S.ror + ((size(rot,1)-d+0.5) - size(rot,1)/2)*S.step;
			dist = S.ror + ((d+0.5) - size(rot,1)/2)*S.step;
            if (dist < 0)
                error('Detector crossed the volume!');
            end
			sigma = S.psf_a + S.psf_b * dist;
			gauss = 1/(sqrt(2*pi)*sigma) * exp( -(((1:pc)-(pc/2)) * S.step ).^2 / 2 /sigma^2 );
			gauss = gauss(find(gauss>1e-12));
			gc = size(gauss,2);
			sp = conv( rot(d,:), gauss);
			rot(d,:) = sp( (floor(gc/2)+1):(pc+floor(gc/2)) );
		end
	end
	proj = sum(rot);
	if (S.sigma_intr)
		pc = size(proj,2);
		gauss = 1/(sqrt(2*pi)*S.sigma_intr) * exp( -(((1:pc)-floor(pc/2)) * S.step ).^2 / 2 /S.sigma_intr^2 );
		gauss = gauss(find(gauss>1e-12));
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
