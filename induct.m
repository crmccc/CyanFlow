
a23=1.3761-4.582i;
a21=1.25-3.75i;
a14=3.33-10i;
res=zeros(4,4);
res(1,1)=res(1,1)+(a21+a14);
res(1,2)=res(1,2)-a21;
res(1,4)=res(1,4)-a14;
res(2,2)=res(2,2)+a23+a21;
res(2,3)=res(2,3)-a23;
res(2,1)=res(2,1)-a21;
res(3,3)=res(3,3)+a23;
res(4,1)=res(4,1)-a14;
res(4,4)=res(4,4)+a14;
