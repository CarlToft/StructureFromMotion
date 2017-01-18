function display(i)
%IMAGEDATA/DISPLAY displays imagedata object

disp(' ');
disp([inputname(1),' = ']);
disp(' ');

disp('  Imagedata:');

if ~isempty(i.filename),
  disp(['   Image file: ',i.filename]);
end
if ~isempty(i.im)
  disp(['   Size of image: ',num2str(size(i.im,1)),'x',num2str(size(i.im,2))]);
end

np=size(i.points,2);
nl=size(i.lines,2);
nc=size(i.conics,2);
if np==0 & nl==0 & nc==0,
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
  if nc>0,
   if nc==1,
     disp(['   1 conic']);
   else
     disp(['   ',num2str(nc),' conics']);
   end
  end
end

disp(' ');


