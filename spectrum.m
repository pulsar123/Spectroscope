# Script to plot a (calibrated) spectrum

dirname = '' # Working directory with the trailing slash, e.g. 'E:\ASTRO\Spectra\'. Leave empty for current folder

use_response_curve = 0;  # Whether to use response curve, tr_final.dat

remove_continuum = 1; # Remove continuum spectrum using a polynome
n_polynome = 6;  # polynome degree
show_all = 0; # show both the original spectrum, polynome fitting, and the corrected spectrum

clf;
hold on;

N = 1;
dy = [];


for i=1:N

filename = [dirname 'spectrum.dat'];  # Input spectrum file 
#filename = [dirname 'spectrum.dat_d' num2str(i)];
separator = '';
skipped_rows = 0;
skipped_columns = 0;
m = dlmread(filename, separator, skipped_rows, skipped_columns);
filenamep = [dirname 'params.txt'];  # Calibration file
p = dlmread(filenamep, separator, skipped_rows, skipped_columns);
# Converting pixel coordinates to nanometers using calibration parameters p:
evalc('lambda = polyval (p, m(:,1))');

#y = m(:, 2);
y = log10(m(:, 2));
y = y - max(y);
if (i==1)
  y0 = y;
endif
dy(i) = sum(y-y0)/rows(y);  # log difference between different plots (compared to plot1)

  y_copy = y;
#plot(lambda, y_copy, "r");

if (use_response_curve)
  filenamec = [dirname 'tr_final.dat'];  # Camera response curve polynome file
  pc = dlmread(filenamec, separator, skipped_rows, skipped_columns);
  lam_tr = pc(:,1);
  tr = pc(:,2);
  tr = tr - max(tr);
  tr_interp = interp1 (lam_tr, tr, lambda, "extrap");
  y = y - tr_interp;  # Multiplying by inverse of the camera sensitivity, to remove the camera factor
endif

I = y;
#I = log10(y);
#I = y - dy(i);

xlim([375 725]);
x1=min(lambda);
x2=max(lambda);
xlim([x1-25 x2+25]);
xlabel('lambda, nm');
#ylabel('I');
ylabel('log10(I)');

if (show_all == 1 || remove_continuum == 0)
#  plot(lambda, I, "g");  # Plotting the spectrum
  plot(lambda, I, 'color', [(i-1.0)/(N-1.0)  0.5  1-(i-1.0)/(N-1.0)]);
endif


if (remove_continuum == 1)
  [a s] = polyfit(lambda,I,n_polynome);
  evalc('y = polyval (a, lambda)');
  if (show_all == 1)
    plot(lambda, y, "r");
    endif
  plot(lambda, I-y, "b");
#  plot(lambda, 10.^(I), "b");
  endif

endfor  

  
yy=ylim();
y1=yy(1);  y2=yy(2);

lam = 440;
line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "y")
lam = 680;
line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "y")


# Plotting main nebula spectral lines
lam=687.55; # O2
line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "r")
lam=656.3; # H
line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "g")
lam=628.15; # O2
line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "r")
lam=589.3; # Na
line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "y")
lam=527; # F
line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1)
lam=486.1; # H
line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "g")
lam=518.4; # Mg
line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "magenta")
lam=517.3; # Mg
line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "magenta")
lam=430.8; # F
line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1)
lam=393.4; # Ca
line ("xdata",[lam,lam], "ydata",[y1,y2], "linewidth", 1, "color", "cyan")

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
