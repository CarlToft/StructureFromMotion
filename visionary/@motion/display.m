function display(m)
%MOTION/DISPLAY displays motion object
disp(' ');
disp([inputname(1),' = ']);
disp(' ');
nbr = length(m.cam);
  disp(['  Motion:']);
if nbr == 1,
  disp(['   1 camera']);
else
  disp(['   ' num2str(nbr),' cameras']);
end

disp(' ');

