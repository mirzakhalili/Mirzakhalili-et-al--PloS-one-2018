% Code by Ehsan Mirzakhalili mirzakh@umich.edu
% https://doi.org/10.1371/journal.pone.0201302
clc;
clear;
close all;
clf;

Cut=0;
TEnd=60;
cd Data
load('Day1-Control.mat');
cd ..
x=Data(1,:);
Data(1,:)=[];
y=mean(Data,1)';
s=std(Data,0,1)';
y(x<Cut)=[];
s(x<Cut)=[];
x(x<Cut)=[];
y(x>TEnd)=[];
s(x>TEnd)=[];
x(x>TEnd)=[];
cd Main
COLORS;
Par=Init;
M=Run(x,Par);
cd ..
MinY=-5;
MaxY=12;

h1=plot(x,M);
ax=gca;

set(ax, 'NextPlot', 'add');
h2=plot(x,y);

set(ax, 'NextPlot', 'add');
syms Stim
Scale=(MaxY-MinY)/20;Offset=MinY+1;
S=fplot((heaviside(Stim-10)-heaviside(Stim-40))*Scale+Offset,[0,60]);
text(21,-2.5,'Stimulus','Color',ST,'Fontsize',14);
box off;

set(h1,'color',MO,'linewidth',2);
set(h2,'color',F2,'linestyle','none','marker','.','markersize',12);
set(S,'color',ST,'linewidth',3);

set(ax,'fontsize',14);
set(ax,'linewidth',2);

ylim([MinY,MaxY]);
yticks(0:5:MaxY);
xlim([0,55]);

xlabel('Time [s]');
ylabel('FRET ratio change [%]');
L=legend('Model','Experiment');L.Box='off';
% legendmarkeradjust(18);

MyT=title('A', 'Units', 'normalized');
set(MyT,'Position', [-0.1 1.00], 'HorizontalAlignment', 'right');