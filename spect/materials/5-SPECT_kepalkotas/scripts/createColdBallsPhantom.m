
function low = createColdBallsPhantom( dim, step, scale )
if (nargin < 2)
	error('Invalid number of arguments, use: createColdBallsPhantom( <image dimension>, <image pixel size>, [scale=4])');
end
if (nargin < 3)
	scale=4;
end

img=zeros([ dim*scale dim*scale ]);
step=step/scale;

balls=[4.75 6.35 7.95 9.55 12.7 15.9];
centerdist=54.8;
cyl=120;

cylr1=centerdist+max(balls(:));
cylr2=centerdist-max(balls(:));

origo=dim*scale*step/2-0.5;

n=0;
h=waitbar(0, 'create phantom');
for iy=1:size(img,2)
for ix=1:size(img,1)
	x=ix*step-origo;
	y=iy*step-origo;
	d = sqrt((x^2 + y^2));
	if ( d < cyl )
		ok=1;
		if ( d <= cylr1 && d >= cylr2)
			for j=1:size(balls,2)
				cx = centerdist * cos(j*2*pi/6);
				cy = centerdist * sin(j*2*pi/6);
				if ( (x-cx)^2 + (y-cy)^2 < balls(j)^2 ) 
					ok = 0; 
					break;
				end
			end
		end
		if (ok == 1)
			img(ix,iy) = 1;
		end
	end
end;
waitbar(iy/(1+size(img,2)) );
end;
close(h);

low=zeros([ dim dim ]);
for iy=1:dim
for ix=1:dim
	low(ix,iy)=sum(sum( img( ((ix-1)*scale+1):((ix)*scale), ((iy-1)*scale+1):((iy)*scale) )));
end;
end;


end
