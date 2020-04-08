%clear the screen
clc
clear
%introduce name and class
fprintf(1,"Jackson Helm\n")
fprintf(1,"1282.02 PAC, 2020\n\n")

%introduce program
fprintf(1,"This program will calculate various quantities of a\n")
fprintf(1,"rectangular channel containing flowing water. See the\n");
fprintf(1,"attached equation sheet for details on calculation.\n\n")

%lay out instructions
fprintf(1,"\nInput the parameter to solve for:\n")
fprintf(1,"1 = Volumetric Flow Rate Q\n")
fprintf(1,"2 = Width W\n")
fprintf(1,"3 = Height H\n")
fprintf(1,"4 = Length L\n")
fprintf(1,"5 = Pressure Differential P\n")
fprintf(1,"6 = Viscosity V\n\n")

%prompt user for choice of solve variable
c=input("Choice?");

%utilize switch case structure to solve for selected variable

switch c
    case 1
        fprintf(1,"You have chosen: Volumetric Flow Rate Q\n\n");
        
        %input known values
        
        w=input("Input Width W:");
        h=input("Input Height H:");
        l=input("Input Length L:");
        p=input("Input Pressure Differential P:");
        v=input("Input Viscosity V:");
        
        %solve for quantity
        
        q=(w*h^3*p)/(12*v*l);
        fprintf(1,"\n\nVolumetric Flow Rate Q = %f cm^3/s",q)
    case 2
        fprintf(1,"You have chosen: Width W\n\n");
        
        %input known values
        
        q=input("Input Volumetric Flow Rate Q:");
        h=input("Input Height H:");
        l=input("Input Length L:");
        p=input("Input Pressure Differential P:");
        v=input("Input Viscosity V:");
        
        %solve for quantity
        
        w = ((12*v*l*q)/(h^3*p));
        fprintf(1,"\n\nWidth W = %f cm",w)
        
    case 3
        fprintf(1,"You have chosen: Height H\n\n");
        
        %input known values
        
        w=input("Input Width W:");
        q=input("Volumetric Flow Rate Q:");
        l=input("Input Length L:");
        p=input("Input Pressure Differential P:");
        v=input("Input Viscosity V:");
        
        %solve for quantity
        
        h=((12*v*q*l)/(w*p))^(1/3);
        fprintf(1,"\n\nHeight H = %f cm",q)
        
    case 4
        
        fprintf(1,"You have chosen: Length L\n\n");
        
        %input known values
        
        w=input("Input Width W:");
        h=input("Input Height H:");
        q=input("Input Volumetric Flow Rate Q:");
        p=input("Input Pressure Differential P:");
        v=input("Input Viscosity V:");
        
        %solve for quantity
        
        l=(w*h^3*p)/(12*v*q);
        fprintf(1,"\n\nLength L = %f cm",l)
        
    case 5
        
        fprintf(1,"You have chosen: Pressure Differential P\n\n");
        
        %input known values
        
        w=input("Input Width W:");
        h=input("Input Height H:");
        l=input("Input Length L:");
        q=input("Input Volumetric Flow Rate Q:");
        v=input("Input Viscosity V:");
        
        %solve for quantity 
        
        p=(12*v*l*q)/(w*h^3);
        fprintf(1,"\n\nPressure Differential = %f dyne/cm^2",p)
    case 6
        fprintf(1,"You have chosen: Viscosity V\n\n");
        
        %input known values
        
        w=input("Input Width W:");
        h=input("Input Height H:");
        l=input("Input Length L:");
        p=input("Input Pressure Differential P:");
        q=input("Input Volumeetric Flow Rate Q:");
        
        %solve for quantity
        
        v=(w*h^3*p)/(12*l*q);
        fprintf(1,"\n\nViscosity = %f g/cm-s",v)
end
        
%calculate average velocity

a = (h^2*p)/(12*v*l);
fprintf(1,"\n\nAverage Velocity Vavg = %f m/s",a)

%calculate shear stress at wall

t = ((h/2)*p)/l;
fprintf(1,"\n\nShear Stress at Wall T = %f dyne/cm^2",t)

%calculate Effective Diameter

d = (4*w*h)/((2*w)+(2*h));
fprintf(1,"\n\nEffective Diameter = %f cm",d)

%calculate Reynolds number

r = (1*a*d)/(v);
fprintf(1,"\n\nReynolds Number = %f",r)

%calculate Entrance Length

e = 0.06*r*d;
fprintf(1,"\n\nEntrance Length = %f cm",e)

%set up table 
Parameter = {'Volumetric Flow Rate Q';'Width W';'Height H';'Length L';'Pressure Differential';'Viscosity';'Average Velocity';'Wall Shear Stress';'Reynolds Number';'Entrance Length'};
Value = {q;w;h;l;p;v;a;t;r;e};
Unit = {'cm^3/s';'cm';'cm';'cm';'dyn/cm^2';'g/cm*s';'cm/s';'dyn/cm^2';'n/a';'cm';};

%create table
X=table(Parameter,Value,Unit)

%set up y values for location inside channel

y = linspace(-h/2,h/2,1000);

%set up velocity function

vinside = (p/8*v*l)*(h^2-(4*y.^2));

%set up shear stress function

ssinside = (p*(abs(y)/l));

%set up plot 1 for velocity
figure(1);
plot(y,vinside)
hold on
title('Velocity at distance away from center');
xlabel('y distance (cm)');
ylabel('velocity (cm/s)');
hold off

figure(2);
plot(y,ssinside)
hold on
title('Shear stress at distance away from center');
xlabel('y distance (cm)');
ylabel('shear stress (dyn/cm^2)');
hold off










