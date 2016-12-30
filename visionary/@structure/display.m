function display(s)
%STRUCTURE/DISPLAY displays structure object

disp(' ');
disp([inputname(1),' = ']);
disp(' ');
disp('  Structure:');
np=size(s.points,2);
nl=size(s.lines,2);
nq=size(s.quadrics,2);
if np==0 & nl==0 & nq==0,
  disp('   No features');
else
  if np>0,
   if np==1,
     disp(['   1 point']);
   else
     disp(['   ',num2str(np),' points']);
   end
  end
  if nl>0,
   if nl==1,
     disp(['   1 line']);
   else
     disp(['   ',num2str(nl),' lines']);
   end
  end
  if nq>0,
   if nq==1,
     disp(['   1 quadric']);
   else
     disp(['   ',num2str(nq),' quadrics']);
   end
  end
end

disp(' ');

