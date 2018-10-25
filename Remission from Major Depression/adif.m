function c=adif(a,b)
 if abs(a-b)>pi
  c=2*pi-abs(a-b);
 else
  c=abs(a-b);
 end
