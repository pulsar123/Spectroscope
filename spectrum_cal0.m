# Calibration script for a DIY spectrometer. The input spectrum file spectrum.dat has two columns - 
# x pixel coordinate, and pixel brightness. It can be produced from an image of the spectrum using
# ImageMacick package and bash scripting, like this (assuming the image width is 1280 pixels):

#   convert input_image  -resize 1280x1\! __a.txt
#   cat __a.txt  |cut -d\) -f1 |sed 's/,0: (/ /g' |grep -v \# > spectrum.dat


pkg load optim
pkg load specfun
format long;

# Different stages:
# 1: Initial stage. First calibration. Will read target spectrum.dat, detect main spectral lines,
#    write file calibr_measured.txt which contains pixel coordinates of the detected lines. It also plots
#    the raw spectrum. You use the plot to identify some/all lines in the spectrum, and add the actual
#    wavelengths on the corresponding lines in calibr_measured.txt as a second column. 
#    Delete all the lines without an identification.
# 2: Second stage. Reads the file calibr_measured.txt, computes polynomial fitting parameters, saves them
#    in params0.txt file.
# 3: Third (recalibration, optional) stage. Reads initial parameters guess params0.txt, current calibration target image,
#    finds all spectral lines, matches them to a list of known spectral lines (second column in calibr_measured.txt)
#    computes and writes improved calibration parameters in params.txt.
stage = 1;

dirname = '' # Working directory with the trailing slash, e.g. 'E:\ASTRO\Spectra\'. Leave empty for current folder
xsign = -1;  # Set to -1 if your red spectrum side is to the left; 1 otherwise
# Spectral line identification parameters:
del0 = 0.2;  # fractional decrease vs the peak value to draw line borders
dx = 20; # how far from the peak in x direction the line (within del0 brightness) can stretch 
noise =5e-4; # no lines with peak value below this times the max brightness in the plot
sgm0 = 2; # Initial sigma value in pixels for Gaussian peak fitting
dl_max = 20; # If a nearest detected line is farther than this many nm from an actual spectral line, skip this actual line
dlk_max = 2.6; # Rejecting lines deviating by more than this many nm in the new fit

n_poly = 4;  # Polynome order (>=2) for calibration
debug = 1;
skewed_gauss = 0; # If 1, use skewed gaussians for line top fitting

clf;
hold on;

# The entire stage 2 stuff:
if (stage == 2)
  # Reading pixel coordinates and corresponding actual wavelengths for calibration target 
  # spectral lines from this file (output of stage 1):
  filename = [dirname 'calibr_measured.txt'];
  separator = '';
  skipped_rows = 0;
  skipped_columns = 0;
  m = dlmread(filename, separator, skipped_rows, skipped_columns);
  x = m(:,1);  # Pixel coordinates
  lambda = (m(:, 2));  #  Actual wavelengths, nm
  scatter(x, lambda, "b", "filled"); # Blue dots are the data
  xlabel('Pixel coordinate')
  ylabel('lambda, nm')
  [p s] = polyfit(x,lambda,n_poly);
  printf("Fitting parabola parameters:\n")
  p
  sd2 = 0;
  evalc('y0 = polyval (p, x)');
  for k=1:rows(x)
    sd2 = sd2 + (lambda(k)-y0(k))^2;
  endfor
  std = sqrt(sd2/(rows(x) - 1));
  printf("std (nm): %f\n", std);
  # Writing calibration parameters to a file:
  printf("Uncertainties of the parameters:\n");
  sqrt (diag (s.C)/s.df)*s.normr
  r=xlim();
  x2=r(1):1:r(2);
  evalc('y2 = polyval (p, x2)');
  plot(x2,y2,"r");  # Red line is the parabolic fit to the data
  printf("Writing parameters to file params0.txt...\n");
  dlmwrite([dirname 'params.txt'], p,' ');
  dlmwrite([dirname 'params0.txt'], p,' ');
  hold off
  return;
endif


# Raw calibration spectrum (two columns - pixel coordinate, intensity):
filename = [dirname 'spectrum.dat'];
# empty separator means 'automatic'
separator = '';
skipped_rows = 0;
skipped_columns = 0;
m = dlmread(filename, separator, skipped_rows, skipped_columns);

lambda = xsign*m(:,1);
I = log10(m(:, 2));

noise_floor = noise * max(m(:,2));
line_id = 0;
w=[]; s=[]; a=[];
dirs=[-1, 1];
edges=[0,0];

# Finding all spectral lines:
for i=3:rows(m)-2
  x0=m(i,1);
  y0=m(i,2);
  if (m(i,2)>=m(i-1,2) && m(i,2)>m(i+1,2) && m(i,2)>noise_floor)
    # We found a local maximum
    # Finding the extent of the line in both (positive and negative) directions
    good = 1;
    edges=[i,i];
    for j=1:2
      dir = dirs(j);
      i1 = i;
      # Covering the allowed ranges of x:
      while (i1>2 && i1<rows(m)-1)
        i1 = i1 + dir;
        if (m(i1,2) > m(i,2) || abs(i1-i)>dx)
          # Bad line: found brighter pixels then the peak, or hit dx limits
          good = 0;
          break;
          endif

        edges(j) = i1;

        if (m(i1,2)<m(i,2)*(1-del0))
          # We hit the edge of the line
          break;
          endif 
        endwhile
      
      if (good == 0)
        break;
        endif
      endfor  #  j/dir loop
  
    if (good == 1)  
      line_id = line_id + 1;
      # The line is good, and we found both edges at this point
      xx=xsign*m(edges(1):edges(2),1);
      yy=m(edges(1):edges(2),2);
      if (skewed_gauss == 0)
      # Nonlinear Gaussian fitting to the peak:
      f = @ (p) p(1) * exp (-(p(2)-xx).^2/2/p(3).^2) - yy;
      init = [y0; xsign*x0; sgm0];
      # Using evalc to suppress output to screen:
      evalc('[p, residuals, cvg, outp] = nonlin_residmin (f, init)');
      w(line_id) = xsign*p(2);  # coordinate of the center of the line
      s(line_id) = abs(p(3));  # sigma (half-width) of the line
      a(line_id) = log10(p(1));  # amplitude of the line
      xx1=xsign*(m(edges(1),1)-1:0.1:m(edges(2),1)+1);
      yy1 = log10(p(1) * exp (-(p(2)-xx1).^2/2/p(3).^2));
        else
      # Nonlinear two-sided Gaussian fitting to the peak:
      f = @ (p) p(1) * (exp (-(p(2)-xx).^2/2/p(3).^2).*heaviside(xx-p(2)) + exp (-(p(2)-xx).^2/2/p(4).^2).*heaviside(p(2)-xx)  ) - yy;
      init = [y0; xsign*x0; sgm0; sgm0];
      # Using evalc to suppress output to screen:
      evalc('[p, residuals, cvg, outp] = nonlin_residmin (f, init)');
      w(line_id) = xsign*p(2);  # coordinate of the center of the line
      s(line_id) = abs(p(3));  # left sigma (half-width) of the line
      sr(line_id) = abs(p(4));  # right sigma (half-width) of the line
      a(line_id) = log10(p(1));  # amplitude of the line
      xx1=xsign*(m(edges(1),1)-1:0.1:m(edges(2),1)+1);
      yy1 = log10(p(1) * ( exp (-(p(2)-xx1).^2/2/p(3).^2).*heaviside(xx1-p(2)) + exp (-(p(2)-xx1).^2/2/p(4).^2).*heaviside(p(2)-xx1) ));
      
      endif
#    plot(xx1, yy1, "b","linewidth", 2);
#    scatter(xx, log10(yy), "r","linewidth", 2);
#    pause();
      endif
    endif
  endfor  # i loop

#evalc('lambda1 = polyval (p0, w)');
#scatter(lambda1, a, "g", "filled","linewidth", 2)
  
N = line_id;  # Number of detected lines
printf("Detected %d spectral lines\n", N);
printf("Range of Gaussian fit sigmas (pixels): %f ... %f\n",min(s), max(s));

if (stage == 1)
  xlabel('x');
  ylabel('log10(I)');
  plot(lambda, I, "b");  # Plotting the raw spectrum in blue
  scatter(xsign*w, a, "r", "filled")  # Centers of found spectral lines are shown as red dots
  # Writing the raw coordinates of the detected lines:
  printf("Writing raw lines to calibr_measured.txt ...\n");
  dlmwrite([dirname 'calibr_measured.txt'], [transpose(w) transpose(a)],' ');
#  dlmwrite([dirname 'calibr_measured.txt'], transpose(w),' ');
  hold off;
  return;
endif

    
# Reading the list of actual lines wavelengths for the calibration target (second column)
filename2 = [dirname 'calibr_measured.txt'];
skipped_columns = 1;
lines = dlmread(filename2, separator, skipped_rows, skipped_columns);
N_actual = numel(lines);  # Number of actual lines

# Reading initial guess for the calibration parameters:
filenamep = [dirname 'params0.txt'];
separator = '';
skipped_rows = 0;
skipped_columns = 0;
p0 = dlmread(filenamep, separator, skipped_rows, skipped_columns);
  
# Searching for detected lines matching the actual lines list  
N_matches = 0;
#  Initial guess for the wavelengths of the detected lines:
evalc('lambda = polyval (p0, w)');
#lambda = p0(1)*w.^2 + p0(2)*w + p0(3);
used = zeros(N,1);
used_dl = zeros(N,1);
for i=1:N_actual
  # Finding the nearest detected line
  dl_min = 1e6;
  for j=1:N
    dl = abs(lambda(j)-lines(i)); # distance between detected and actual lines in nm
    if (dl < dl_min)
      dl_min = dl;
      j_min = j;
      endif
    endfor
  if (dl_min > dl_max)
    # Skipping this actual line as it looks like we couldn't detect a good match
    continue;
    endif;
  if (used(j_min) == 0)
    # If the line hasn't been claimed yet
    used(j_min) = i; # Storing the matching actual line id here
    used_dl(j_min) = dl_min;
  else
    # If the line has already been claimed, checking if the new claim is a better one (smaller dl), and replacing it if so
    if (dl_min < used_dl(j_min))
      # Replacing the old detected line with the current, better one
      used(j_min) = i;
      used_dl(j_min) = dl_min;
      endif
    endif
  endfor  #  i , actual lines

# Now that all matches have been found, running a polynomial linfit to get the improved calibration parameters  
# Only lines with matched pairs (used >0) are used
iter = 0;
while (1)
  iter = iter + 1;
printf("===============  Iteration # %d ================\n", iter);
x=[];  y=[];
k = 0; 
for j=1:N
  if (used(j) > 0)
    k = k + 1;
    x(k) = w(j);  # Raw pixel coordinate for the detected line center
    y(k) = lines(used(j));  # Wavelength of the nearest actual line, nm
    endif
  endfor
N_match = k;
printf("N_match = %d\n", N_match);
[p s] = polyfit(x,y,n_poly);
evalc('y0 = polyval (p, x)');
printf("Fitting parameters:\n")  
p
sd2 = 0;
for k=1:N_match
  sd2 = sd2 + (y(k)-y0(k))^2;
endfor
std = sqrt(sd2/(N_match - 1));
printf("std (nm): %f\n", std);
# If some individual lines have too large dl, we remove them, and redo the fitting:
k = 0; dk = 0;
for j=1:N
  if (used(j) > 0)
    k = k + 1;
    dlk = abs(y(k) - y0(k));
    if (dlk > dlk_max)
      # This line deviates too much, removing it
      used(j) = 0;
      used_dl(j) = 0;
      dk = dk + 1;
      endif
    endif
  endfor
if (dk == 0)
  break;
  endif
endwhile  

printf("Uncertainties of the parameters:\n");
sqrt (diag (s.C)/s.df)*s.normr

# Saving the updated calibration parameters:
printf("Writing parameters to file params.txt...\n");
dlmwrite([dirname 'params.txt'], p,' ');

# Polynomic fitting:
evalc('lambda = polyval (p, m(:,1))');
I = log10(m(:, 2));
x1=min(lambda);
x2=max(lambda);
xlim([x1-25 x2+25]);
xlabel('lambda, nm');
ylabel('log10(I)');
plot(lambda, I, "b");

yy=ylim();
y1=yy(1);  y2=yy(2);
for k=1:N_match
  lam = y(k);
  line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "r")  
endfor
# All lines including skipped lines:
for i=1:N_actual
  lam = lines(i);
  line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "r", "linestyle", "--")  
endfor

if (debug == 1)
  evalc('lambda1 = polyval (p, w)');
  scatter(lambda1, a, "g", "filled","linewidth", 2)
endif



# Writing out the wavelengths of the lines found:
fid=fopen([dirname 'tmp.txt'],'w');
for i=1:N
  fprintf(fid, "%f\n", w(i));
  endfor
fclose(fid);

lam=576.959; # 
#line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "cyan")
lam=579.065; # 
#line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "cyan")
lam=546.5; # 
#line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "cyan")
lam=436.6; # 
#line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "cyan")
lam=611.6 ; # 
#line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "green")


hold off;
