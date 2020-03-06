function img=SPECT2DAttenuatedBackwardProj( projection, S, mumap, theta )

% SPECT2DAttenuatedBackwardProj   Calculates the backward projection of a series.
%
%       img=SPECT2DAttenuatedBackwardProj( projection, S, mumap, theta )
%
% Where 
%       projection      is a projection series
%       S               is a SPECT2D System descriptor (see SPECT2DSystem)
%	mumap		is an image sized linear attenuation coefficient map
%       theta           is a vector of projection IDs from the descriptor
%
% Example
%       img = SPECT2DAttenuatedBackwardProj( projection, S, mumap )
%               Creates a backprojection from a full series
%
%       img = SPECT2DAttenuatedBackwardProj( projection, S, mumap, 1 )
%               Creates a backprojection from the first projection image only
%

if (nargin < 3)
        error('Too few arguments: proj = SPECT2DAttenuatedBackwardProj( projection, System, mumap, [angle IDs])');
end
if (nargin < 4)
        theta=1:size(S.theta,2);
end;
if (max(theta) > size(S.theta,2))
        error('Invalid angle ID');
end

img=zeros( S.imagesize, S.imagesize );
if (size(img) ~= size(mumap))
        error('Image has to be the same dimensions as the mumap!');
end

imgsize= size(img);
offset = ceil((sqrt(2)-1)*max(imgsize)/2); 
newsize= imgsize+2*offset;

if ( max(newsize)*S.step > 2*S.ror )
        error('Invalid SPECT System: the whole image must be in the FOV!');
end

for t=1:size(theta,2)
	bigimg = zeros( newsize );

	proj = projection(:,t)';
	if (S.sigma_intr)
		pc = size(proj,2);
		gauss = 1/(sqrt(2*pi)*S.sigma_intr) * exp( -(((1:pc)-floor(pc/2)) * S.step ).^2 / 2 /S.sigma_intr^2 );
		gauss = gauss(find(gauss>1e-12));
		gc = size(gauss,2);
		sp = conv( proj, gauss);
% ezt le kell itt vágni?
%		proj = sp( (floor(gc/2)+1):(pc+floor(gc/2)) );
		proj = sp;
	end

	pc = floor(size(proj,2)/2);
	ac = floor(newsize(1) / 2);
	if (size(proj,2) > newsize(1))
		bigimg( 1, :)=proj( (pc-ac+1):(pc+ac) );
	else
		bigimg( 1, (ac-pc+1):(ac-pc+size(proj,2)))=proj( : );
	end
	for y=(2:newsize(2))
		bigimg(y,:)=bigimg(1,:);
	end

	if (S.psf_a || S.psf_b)
		pc = size(bigimg,2);
		for d=1:size(bigimg,1)
%			dist = S.ror + ((size(rot,1)-d+0.5) - size(rot,1)/2)*S.step;
			dist = S.ror + ((d+0.5) - size(bigimg,1)/2)*S.step;
			sigma = S.psf_a + S.psf_b * dist;
			gauss = 1/(sqrt(2*pi)*sigma) * exp( -(((1:pc)-(pc/2)) * S.step ).^2 / 2 /sigma^2 );
			gauss = gauss(find(gauss>1e-12));
			gc = size(gauss,2);
			sp = conv( bigimg(d,:), gauss);
			bigimg(d,:) = sp( (floor(gc/2)+1):(pc+floor(gc/2)) );
		end
	end

        bigmumap  = zeros( newsize );
        bigmumap( (1+offset):(imgsize(1)+offset), (1+offset):(imgsize(2)+offset) ) = mumap;
        murot = myimrotate( bigmumap, -S.theta(theta(t)), 'bilinear', 'crop' );
        % compute the integral mu-map
        for d=2:size(murot,1)
                murot(d,:)=murot(d,:)+murot(d-1,:);
        end
        % add the attenuation-effect to the rotated map
        for d=1:size(murot,1)
                bigimg(d,:)=bigimg(d,:).*exp(-S.step*murot(d,:)/10);
        end

	rot = myimrotate( bigimg, S.theta(theta(t)), 'bilinear', 'crop' );
%	rot = rot * sum(bigimg(:))/sum(rot(:));

	img = img + rot( (1+offset):(imgsize(1)+offset), (1+offset):(imgsize(2)+offset) );
end

end
