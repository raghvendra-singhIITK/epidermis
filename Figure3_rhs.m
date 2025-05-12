function ydot=Figure3_rhs(t, y)
ydot=zeros(4, 1);

ANotch=y(1);
AWntSC=y(2);
AWntTA=y(3);
AYAP=y(4);

kNotch=1; 
kdNotch=0.2;
kpNotch=0.05;

kWntSC=1.4;
kdWntSC=0.2;
kpWntSC=0.1;

kWntTA=1;
kdWntTA=0.2;
kpWntTA=0.1;

kYAP=0.188; 
kdYAP=0.02;
kpYAP=0.1;

n=2; 

ydot(1)=kNotch*AWntSC^n/(1+AWntSC^n+AYAP^n+AWntTA^n) - kdNotch*ANotch+ kpNotch;
ydot(2)=kWntSC*ANotch^n/(1+ANotch^n)- kdWntSC*AWntSC+kpWntSC;         
ydot(3)=kWntTA/(1+ANotch^n)- kdWntTA*AWntTA+kpWntTA;
ydot(4)=kYAP*ANotch^n/(1+ANotch^n)- kdYAP*AYAP+kpYAP;




                                                                            