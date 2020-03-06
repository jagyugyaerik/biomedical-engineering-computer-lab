%% -*- texinfo -*-
%% @defun @var{recon} = iradon (@var{proj}, @var{theta}, @var{interp}, @
%%                              @var{filter}, @var{scaling}, @var{output_size})
%%
%% Performs filtered back-projection on the projections in @var{proj}
%% to reconstruct an approximation of the original image.
%%
%% @var{proj} should be a matrix whose columns are projections of an
%% image (or slice).  Each element of @var{theta} is used as the angle
%% (in degrees) that the corresponding column of @var{proj} was
%% projected at.  If @var{theta} is omitted, it is assumed that
%% projections were taken at evenly spaced angles between 0 and 180 degrees.
%% @var{theta} can also be a scalar, in which case it is taken as the
%% angle between projections if more than one projection is provided.
%% 
%% @var{interp} determines the type of interpolation that is used
%% in the back-projection.  It must be one of the types accepted by
%% @command{interp1}, and defaults to 'Linear' if it is omitted.
%%
%% @var{filter} and @var{scaling} determine the type of rho filter 
%% to apply.  See the help for @command{rho_filter} for their use.
%%
%% @var{output_size} sets the edge length of the output image (it
%% is always square).  This argument does not scale the image.  If it
%% is omitted, the length is taken to be
%% @group
%% 2 * floor (size (proj, 1) / (2 * sqrt (2))).
%% @end group
%% 
%% If @var{proj} was obtained using @command{radon}, there is no
%% guarantee that the reconstructed image will be exactly the same
%% size as the original.
%% 
%% @defunx [@var{recon}, @var{filt}] = iradon (...)
%%
%% This form also returns the filter frequency response in the vector
%% @var{filt}.
%% @end defun
%%
%% Performs filtered back-projection in order to reconstruct an
%% image based on its projections.
%%
%% Filtered back-projection is the most common means of reconstructing
%% images from CT scans.  It is a two step process: First, each of 
%% the projections is filtered with a `rho filter', so named due
%% to its frequency domain definition, which is simply |rho|, where
%% rho is the radial axis in a polar coordinate system.  Second, 
%% the filtered projections are each `smeared' across the image
%% space.  This is the back-projection part.
%%
%% Usage example:
%% @example
%%   P = phantom ();
%%   projections = radon (P, 1:179);
%%   reconstruction = iradon (filtered_projections, 1:179, 'Spline', 'Hann');
%%   figure, imshow (reconstruction, [])
%% @end example

%% Copyright (C) 2010 Alex Opie <lx_op@orcon.net.nz>
%%
%%
%% This program is free software; you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 3 of the License, or (at
%% your option) any later version.
%%
%% This program is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%% General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; see the file COPYING.  If not, see
%% <http://www.gnu.org/licenses/>.


function [recon, filt] = iradon (proj, theta, interp, filter, scaling, output_size)
  
  if (nargin == 0)
    error ('No projections provided to iradon');
  end
  
  if (nargin < 6)
    output_size = 2 * floor (size (proj, 1) / (2 * sqrt (2)));
  end
  if (nargin < 5) || (length (scaling) == 0)
    scaling = 1;
  end
  if (nargin < 4) || (length (filter) == 0)
    filter = 'Ram-Lak';
  end
  if (nargin < 3) || (length (interp) == 0)
    interp = 'linear';
  end
  if (nargin < 2) || (length (theta) == 0)
    theta = 180 * (0:1:size (proj, 2) - 1) / size (proj, 2);
  end
  
  if (isscalar (theta)) && (size (proj, 2) ~= 1)
    theta = (0:size (proj, 2) - 1) * theta;
  end
  
  if (length (theta) ~= size (proj, 2))
    error ('iradon: Number of projections does not match number of angles')
  end
  if (isscalar (scaling)==0)
    error ('iradon: Frequency scaling value must be a scalar');
  end
  if (length (find (strcmpi (interp, {'nearest', 'linear', 'spline', ...
                                       'pchip', 'cubic'})))==0)
    error ('iradon: Invalid interpolation method specified');
  end
  
  %% Convert angles to radians
  theta = theta * pi / 180;
  
  %% First, filter the projections
  [filtered, filt] = rho_filter (proj, filter, scaling);
  
  %% Next, back-project
  recon = back_project (filtered, theta, interp, output_size);
  
end


function recon = back_project (proj, theta, interpolation, dim)
  %% Make an empty image
  recon = zeros (dim, dim);
  
  %% Zero pad the projections if the requested image
  %% has a diagonal longer than the projections
  diagonal = ceil (dim * sqrt (2)) + 1;
  if (size (proj, 1) < diagonal)
    diff = 2 * ceil ((diagonal - size (proj, 1)) / 2);
    proj = padarray (proj, diff / 2);
  end
  
  %% Create the x & y values for each pixel
  centre = floor ((dim + 1) / 2);
  x = (0:dim - 1) - centre + 1;
  x = repmat (x, dim, 1);
   
  y = (dim - 1: -1 : 0)' - centre;
  y = repmat (y, 1, dim);
  
  %% s axis for projections, needed by interp1
  s = (0:size (proj, 1) - 1) - floor (size (proj, 1) / 2);
  
  %% Sum each projection's contribution
  for i = 1:length (theta)
    s_dash = (x * cos (theta (i)) + y * sin (theta (i)));
    interpolated = interp1 (s, proj (:, i), s_dash (:), ['*', interpolation]);
    recon = recon + reshape (interpolated, dim, dim);
  end
  
  %% Scale the reconstructed values to their original size
  recon = recon * pi / (2 * length (theta));
  
end


function [filtered_proj, filt] = rho_filter (proj, type, scaling)

  filtered_proj = proj;
  
  if (nargin < 3)
    scaling = 1;
  end
  if (nargin < 2) || (size (type, 2) == 0)
    type = 'ram-lak';
  end

  %% Extend the projections to a power of 2
  new_len = 2 * 2^nextpow2 (size (filtered_proj, 1));
  
  %% Scale the frequency response
  int_len = (new_len * scaling);
  if (mod (floor (int_len), 2))
    int_len = ceil (int_len);
  else
    int_len = floor (int_len);
  end

  if (strcmpi (type, 'none'))
      filtered_proj = proj;
      filt=zeros([1 int_len]);
      filt(1)=1;
    return;
  end
  
  if (scaling > 1) || (scaling < 0)
    error ('Scaling factor must be in [0,1]');
  end

  filtered_proj (new_len, 1) = 0;
  
  %% Create the basic filter response
  rho = scaling * (0:1 / (int_len ):1);
  rho = rho(1:int_len)';
  
  %% Create the window to apply to the filter response

  if (strcmpi (type, 'ram-lak'))
    filt = ones(size(rho));
  elseif (strcmpi (type, 'hamming'))
    filt = fftshift (hamming (length (rho)));
  elseif (strcmpi (type, 'hann'))
    filt = fftshift (hann (length (rho)));
  elseif (strcmpi (type, 'cosine'))
    f = 0.5 * (0:length (rho) - 1)' / length (rho);
    filt = fftshift (sin (2 * pi * f));
  elseif (strcmpi (type, 'shepp-logan'))
    f = (0:length (rho) / 2)' / length (rho);
    filt = sin (pi * f) ./ (pi * f);
    filt (1) = 1;
    filt = [filt; filt((size (filtered_proj, 1) - 1):-1:2)];
  else
    error ('rho_filter: Unknown window type');
  end

  %% Apply the window
  filt = rho.*filt;
  
  %% Pad the response to the correct length
  len_diff = new_len - int_len;
  if (len_diff ~= 0)
    pad = len_diff / 2;
    filt = padarray (fftshift (filt), pad);
    filt = fftshift (filt);
  end
  
  filtered_proj = fft (filtered_proj);
  %% Perform the filtering
  for i = 1:size (filtered_proj, 2)
    filtered_proj (:, i) = filtered_proj (:, i) .* filt;
  end
  
  %% Finally bring the projections back to the spatial domain
  filtered_proj = real (ifft (filtered_proj));
  
  %% Chop the projections back to their original size
  filtered_proj = filtered_proj ( 1 : size (proj, 1), :);
  
end


function yi = interp1(x, y, xi, method, extrap)

  if nargin<3 || nargin>5
    usage('yi = interp1(x, y, xi [, "method" [, "extrap"]])');
  end

  if nargin < 4, 
    method = 'linear';
  else
    method = lower(method); 
  end

  if nargin < 5
    extrap = NaN;
  end

  %% reshape matrices for convenience
  x = x(:);
  if size(y,1)==1, y=y(:); end
  transposed = (size(xi,1)==1);
  xi = xi(:);

  %% determine sizes
  nx = size(x,1);
  [ny, nc] = size(y);
  if (nx < 2 || ny < 2)
     error ('interp1: table too short');
  end

  %% determine which values are out of range and set them to extrap,
  %% unless extrap=='extrap' in which case, extrapolate them like we
  %% should be doing in the first place.
  minx = x(1);
  if (method(1) == '*')
     dx = x(2) - x(1);
     maxx = minx + (ny-1)*dx;
  else
     maxx = x(nx);
  end
  if strcmp(extrap,'extrap')
    range=1:size(xi,1);
    yi = zeros(size(xi,1), size(y,2));
  else
    range = find(xi >= minx & xi <= maxx);
    yi = extrap*ones(size(xi,1), size(y,2));
    if isempty(range), 
      if transposed, yi = yi.'; end
      return; 
    end
    xi = xi(range);
  end

  if strcmp(method, 'nearest')
    idx = lookup(0.5*(x(1:nx-1)+x(2:nx)), xi)+1;
    yi(range,:) = y(idx,:);

  elseif strcmp(method, '*nearest')
    idx = floor((xi-minx)/dx+1.5);
    yi(range,:) = y(idx,:);

  elseif strcmp(method, 'linear')
    %% find the interval containing the test point
    idx = lookup (x(2:nx-1), xi)+1; 
                                % 2:n-1 so that anything beyond the ends
                                % gets dumped into an interval
    %% use the endpoints of the interval to define a line
    dy = y(2:ny,:) - y(1:ny-1,:);
    dx = x(2:nx) - x(1:nx-1);
    s = (xi - x(idx))./dx(idx);
    yi(range,:) = s(:,ones(1,nc)).*dy(idx,:) + y(idx,:);

  elseif strcmp(method, '*linear')
    %% find the interval containing the test point
    t = (xi - minx)/dx + 1;
    idx = floor(t);

    %% use the endpoints of the interval to define a line
    dy = [y(2:ny,:) - y(1:ny-1,:); zeros(1,nc)];
    s = (t - idx)./dx;
    yi(range,:) = s(:,ones(1,nc)).*dy(idx,:) + y(idx,:); 

  elseif strcmp(method, 'pchip') || strcmp(method, '*pchip')
    if (nx == 2) x = linspace(minx, maxx, ny); end
    yi(range,:) = pchip(x, y, xi);

  elseif strcmp(method, 'cubic')
    if (nx < 4 || ny < 4)
      error ('interp1: table too short');
    end
    idx = lookup(x(3:nx-2), xi) + 1;

    %% Construct cubic equations for each interval using divided
    %% differences (computation of c and d don't use divided differences
    %% but instead solve 2 equations for 2 unknowns). Perhaps
    %% reformulating this as a lagrange polynomial would be more efficient.
    i=1:nx-3;
    J = ones(1,nc);
    dx = diff(x);
    dx2 = x(i+1).^2 - x(i).^2;
    dx3 = x(i+1).^3 - x(i).^3;
    a=diff(y,3)./dx(i,J).^3/6;
    b=(diff(y(1:nx-1,:),2)./dx(i,J).^2 - 6*a.*x(i+1,J))/2;
    c=(diff(y(1:nx-2,:),1) - a.*dx3(:,J) - b.*dx2(:,J))./dx(i,J);
    d=y(i,:) - ((a.*x(i,J) + b).*x(i,J) + c).*x(i,J);
    yi(range,:) = ((a(idx,:).*xi(:,J) + b(idx,:)).*xi(:,J) ...
                   + c(idx,:)).*xi(:,J) + d(idx,:);

  elseif strcmp(method, '*cubic')
    if (nx < 4 || ny < 4)
      error ('interp1: table too short');
    end

    %% From: Miloje Makivic 
    %% http://www.npac.syr.edu/projects/nasa/MILOJE/final/node36.html
    t = (xi - minx)/dx + 1;
    idx = max(min(floor(t), ny-2), 2);
    t = t - idx;
    t2 = t.*t;
    tp = 1 - 0.5*t;
    a = (1 - t2).*tp;
    b = (t2 + t).*tp;
    c = (t2 - t).*tp/3;
    d = (t2 - 1).*t/6;
    J = ones(1,nc);
    yi(range,:) = a(:,J) .* y(idx,:) + b(:,J) .* y(idx+1,:) ...
                  + c(:,J) .* y(idx-1,:) + d(:,J) .* y(idx+2,:);

  elseif strcmp(method, 'spline') || strcmp(method, '*spline')
    if (nx == 2) x = linspace(minx, maxx, ny); end
    yi(range,:) = spline(x, y, xi);

  else
    error(['interp1 does not understand method "', method, '"']);
  end
  if transposed, yi=yi.'; end

end


%!demo
%! P = phantom ();
%! figure, imshow (P, []), title ("Original image")
%! projections = radon (P, 0:179);
%! reconstruction = iradon (projections, 0:179, 'Spline', 'Hann');
%! figure, imshow (reconstruction, []), title ("Reconstructed image")

% $Log$
% 2010-30-03 lxop
% First submitted to Octave-Forge
