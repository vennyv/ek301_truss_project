%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input C matrix
C= [1 1 0 0 0 0 0 0 0 0 0 0 0;1 0 1 1 0 0 0 0 0 0 0 0 0;0 1 1 0 1 1 0 1 0 0 0 0 0;0 0 0 1 1 0 1 0 0 0 0 0 0;0 0 0 0 0 1 1 0 1 1 0 0 0; 0 0 0 0 0 0 0 1 1 0 1 1 0;0 0 0 0 0 0 0 0 0 1 1 0 1; 0 0 0 0 0 0 0 0 0 0 0 1 1];

sizes=size(C);
joints=sizes(1,1);
members=sizes(1,2);

%Input Sx and Sy matrices
Sx=zeros(8,3);
Sx(1,1)=1;

Sy=zeros(8,3);
Sy(1,2)=1;
Sy(8,3)=1;

%Input X and Y vectors
X=[0;6.25;12.5;12.5;19.5;19.5;25.25;31];
Y=[0;6;0;12;12;0;6;0];

%Input L vector
L=zeros(16,1);
L(11)=32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setting up A vector
A=zeros(joints*2,members+3);
A(1:joints,1:members)=C;
A(joints+1:joints*2,1:members)=C;

[a b]=find(C==1);


%Calculate total distance
distvec=zeros(13,1);
for x = 1:2:members*2
    xdist=(X(a(x+1))-X(a(x)));
    xdistsq=xdist^2;
    ydist=(Y(a(x+1))-Y(a(x)));
    ydistsq=ydist^2;
    if x==1
       distvec(x)=sqrt(xdistsq+ydistsq);
    else 
        distvec(x-1)=sqrt(xdistsq+ydistsq);
end
end

distvec=nonzeros(distvec');
tdist=sum(distvec);

%Calculate total cost
cost=10*joints+tdist;

%Calculating A vector
for x = 1:2:members*2
    xdist=(X(a(x+1))-X(a(x)));
    xdistsq=xdist^2;
    ydist=(Y(a(x+1))-Y(a(x)));
    ydistsq=ydist^2;
    dist=sqrt(xdistsq+ydistsq);
    A(a(x),b(x))=(X(a(x+1))-X(a(x)))/dist;
    A(a(x+1),b(x+1))=(X(a(x))-X(a(x+1)))/dist;
end

for x = 1:2:members*2
    xdist=(X(a(x+1))-X(a(x)));
    xdistsq=xdist^2;
    ydist=(Y(a(x+1))-Y(a(x)));
    ydistsq=ydist^2;
    dist=sqrt(xdistsq+ydistsq);
    A(a(x)+8,b(x))=(Y(a(x+1))-Y(a(x)))/dist;
    A(a(x+1)+8,b(x+1))=(Y(a(x))-Y(a(x+1)))/dist;
end

A(1:joints,members+1:members+3)=Sx;
A(joints+1:joints*2,members+1:members+3)=Sy;

%Calculate T vector
T=A^-1 * L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Print Results
fprintf('EK301, Section A2, Group 20: Darren S., Venessa M., Vikram B. 4/2/2024.\n');
fprintf('Load: %.1f oz\n',sum(L));
fprintf('Member forces in oz:\n');
for z= 1:members
    if T(z)>0
        fprintf('m%d: %.3f (T)\n',z,abs(T(z)));
    elseif T(z)<0
        fprintf('m%d: %.3f (C)\n',z,abs(T(z)));
    elseif T(z)==0
        fprintf('m%d: %.3f\n',z,abs(T(z)));
    end
end

fprintf('Reaction forces in oz:\n');
fprintf('Sx1: %.2f\n',T(members+1));
fprintf('Sy1: %.2f\n',T(members+2));
fprintf('Sy2: %.2f\n',T(members+3));

fprintf('Cost of truss: $%.2f\n',cost)
fprintf('Theoretical max load/cost ratio in oz/$: 0.241\n')
