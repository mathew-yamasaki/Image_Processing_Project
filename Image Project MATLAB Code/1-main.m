%==========================================================================
%      Course: Image and Signal Processing
%  Instructor: Professor Therese Smith
%     Student: Mathew Yamasaki
%     Project: Final Image Project
%        Date: December 18, 2012
% 
% Description: This program performs image boudary, Freeman's chain code,
%              MPP, boundary length, and boundary diameter calculations.
%              Selected images are diplayed as well as the MPPs, image
%              signatures, and pattern vector plots.
% 
% The following functions Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, 
% & S. L. Eddins are used:
%                           bound2im.m
%                           connectpoly.m
%                           fchcode.m
%                           minperpoly.m
%						    signature.m
%==========================================================================

% read in images from class A
a1=imread('a1.bmp');
a2=imread('a2.bmp');
a3=imread('a3.bmp');
a4=imread('a4.bmp');
a5=imread('a5.bmp');

% read in images from class B
b_21=imread('b21.bmp');
b_22=imread('b22.bmp');
b_23=imread('b23.bmp');
b_24=imread('b24.bmp');
b_25=imread('b25.bmp'); %k5(end+1,:)=0; k5

% read in images from class R
r1=imread('r1.bmp');
r2=imread('r2.bmp');
r3=imread('r3.bmp');
r4=imread('r4.bmp');
r5=imread('r5.bmp');

% class arrays
al={a1,a2,a3,a4,a5};
bl={b_21,b_22,b_23,b_24,b_25};
rl={r1,r2,r3,r4,r5};

% used for figure titles
at={'a1.bmp','a2.bmp','a3.bmp','a4.bmp','a5.bmp'};
bt={'b21.bmp','b22.bmp','b23.bmp','b24.bmp','b25.bmp'};
rt={'r1.bmp','r2.bmp','r3.bmp','r4.bmp','r5.bmp'};

fileNames={'a1.bmp','a2.bmp','a3.bmp','a4.bmp','a5.bmp',...
           'b21.bmp','b22.bmp','b23.bmp','b24.bmp','b25.bmp',...
           'r1.bmp','r2.bmp','r3.bmp','r4.bmp','r5.bmp'};

% array for performing image boundary, Freeman's chain code, boundary
% length, and boundary diameter calculations
letters={a1,a2,a3,a4,a5,b_21,b_22,b_23,b_24,b_25,r1,r2,r3,r4,r5};

%==========================================================================
%                       Part b: show selected images 
%==========================================================================
figure(1);
for i=1:1:length(letters)
    subplot(3,5,i);
    imshow(letters{i});
    title(sprintf('%s',fileNames{i}));
end

%==========================================================================
%                  Part c.1: calculate image boundaries  
%==========================================================================
for i=1:1:length(letters)
    boundaries=bwboundaries(imcomplement(rgb2gray(letters{i})));
end

%==========================================================================
%               Part c.2: calculate Freeman chain-code using fchcode() 
% Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins, Digital 
% Image Processing Using MATLAB, Prentice-Hall
%==========================================================================
for i=1:1:length(letters)
    boundaries=bwboundaries(imcomplement(rgb2gray(letters{i})));
    fchcode(cell2mat(boundaries(1,1)));
end

%==========================================================================
%           Part c.3 & d: calculate and display MPP for letter A  
%==========================================================================
figure(2);
for i=1:length(al)
    b=bwboundaries(imcomplement(rgb2gray(al{i})));
    b=b{1};
    [M,N]=size(al{i});
    xmin=min(b(:,1));
    ymin=min(b(:,2));
    A=imcomplement(rgb2gray(al{i}));
    [x,y]=minperpoly(A,2);
    b2=connectpoly(x,y);
    b2=bound2im(b2,M,N,xmin,ymin);
    b2(:,31:150)=[];
    b2(30:50,:)=[];
    subplot(1,5,i);
    imshow(b2);
    title(sprintf('%s MPP',at{i}));
end

%==========================================================================
%           Part c.3 & d: calculate and display MPP for letter B  
%==========================================================================
figure(3);
for i=1:length(bl)
    b=bwboundaries(imcomplement(rgb2gray(bl{i})));
    b=b{1};
    [M,N]=size(bl{i});
    xmin=min(b(:,1));
    ymin=min(b(:,2));
    A=imcomplement(rgb2gray(bl{i}));
    [x,y]=minperpoly(A,2);
    b2=connectpoly(x,y);
    b2=bound2im(b2,M,N,xmin,ymin);
    b2(:,40:150)=[];
    subplot(1,5,i);
    imshow(b2);
    title(sprintf('%s MPP',bt{i}));
end

%==========================================================================
%           Part c.3 & d: calculate and display MPP for letter R  
%==========================================================================

figure(4);
for i=1:length(rl)
    b=bwboundaries(imcomplement(rgb2gray(rl{i})));
    b=b{1};
    [M,N]=size(rl{i});
    xmin=min(b(:,1));
    ymin=min(b(:,2));
    A=imcomplement(rgb2gray(rl{i}));
    [x,y]=minperpoly(A,2);
    b2=connectpoly(x,y);
    b2=bound2im(b2,M,N,xmin,ymin);
    b2(:,31:150)=[];
    b2(30:50,:)=[];
    subplot(1,5,i);
    imshow(b2);
    title(sprintf('%s MPP',rt{i}));
end

%==========================================================================
%            Part c.4: calculate length of the image boundaries  
%==========================================================================
for i=1:1:length(letters)
    boundaries=bwboundaries(imcomplement(rgb2gray(letters{i})));
    matrix=cell2mat(boundaries);
    [boundary_length dim]=size(matrix);
    
    boundary_length;                         % return image boundary length
    boundary_diameter=max(pdist(matrix));    % return max Euclidean distance
end

%==========================================================================
%                 calculating the pattern vector for A   
%==========================================================================
A1=imcomplement(im2bw(a1,0.5));
A2=imcomplement(im2bw(a2,0.5));
A3=imcomplement(im2bw(a3,0.5));
A4=imcomplement(im2bw(a4,0.5));
A5=imcomplement(im2bw(a5,0.5));

% calculate the boundary
A1B=bwboundaries(A1);
A2B=bwboundaries(A2);
A3B=bwboundaries(A3);
A4B=bwboundaries(A4);
A5B=bwboundaries(A5);

% convert to matrix
A1M=cell2mat(A1B);
A2M=cell2mat(A2B);
A3M=cell2mat(A3B);
A4M=cell2mat(A4B);
A5M=cell2mat(A5B);

[blengthA1 dim]=size(A1M);
[blengthA2 dim]=size(A2M);
[blengthA3 dim]=size(A3M);
[blengthA4 dim]=size(A4M);
[blengthA5 dim]=size(A5M);

diamA1=max(pdist(A1M));
diamA2=max(pdist(A2M));
diamA3=max(pdist(A3M));
diamA4=max(pdist(A4M));
diamA5=max(pdist(A5M));

[d,a]=signature(A1B{1});
a11=moment(d,1);
a12=moment(d,2);
a13=moment(d,3);
a14=moment(d,4);
a15=moment(d,5);

[d,a]=signature(A2B{1});
a21=moment(d,1);
a22=moment(d,2);
a23=moment(d,3);
a24=moment(d,4);
a25=moment(d,5);

[d,a]=signature(A3B{1});
a31=moment(d,1);
a32=moment(d,2);
a33=moment(d,3);
a34=moment(d,4);
a35=moment(d,5);

[d,a]=signature(A4B{1});
a41=moment(d,1);
a42=moment(d,2);
a43=moment(d,3);
a44=moment(d,4);
a45=moment(d,5);

[d,a]=signature(A5B{1});
a51=moment(d,1);
a52=moment(d,2);
a53=moment(d,3);
a54=moment(d,4);
a55=moment(d,5);

% letter B pattern vector
%{
APV=[blengthA1, distA1;
    blengthA2, distA2;
    blengthA3, distA3;
    blengthA4, distA4;
    blengthA5, distA5];
%}
% letter R pattern vector
 
APV=[blengthA1, diamA1, a11, a12, a13, a14, a15;
     blengthA2, diamA2, a21, a22, a23, a24, a25;
     blengthA3, diamA3, a31, a32, a33, a34, a35;
     blengthA4, diamA4, a41, a42, a43, a44, a45;
     blengthA5, diamA5, a51, a52, a53, a54, a55];

% prepare pattern vectors for plotting
ay=[distA1;distA2;distA3;distA4;distA5]; 
ax=[blengthA1;blengthA2;blengthA3;blengthA4;blengthA5];
 
%==========================================================================
%                 calculating the pattern vector for B   
%==========================================================================
B1=imcomplement(im2bw(b_21,0.5));
B2=imcomplement(im2bw(b_22,0.5));
B3=imcomplement(im2bw(b_23,0.5));
B4=imcomplement(im2bw(b_24,0.5));
B5=imcomplement(im2bw(b_25,0.5));

B1B=bwboundaries(B1);
B2B=bwboundaries(B2);
B3B=bwboundaries(B3);
B4B=bwboundaries(B4);
B5B=bwboundaries(B5);

B1M=cell2mat(B1B);
B2M=cell2mat(B2B);
B3M=cell2mat(B3B);
B4M=cell2mat(B4B);
B5M=cell2mat(B5B);

[blengthB1 dim]=size(B1M);
[blengthB2 dim]=size(B2M);
[blengthB3 dim]=size(B3M);
[blengthB4 dim]=size(B4M);
[blengthB5 dim]=size(B5M);

diamB1=max(pdist(B1M));
diamB2=max(pdist(B2M));
diamB3=max(pdist(B3M));
diamB4=max(pdist(B4M));
diamB5=max(pdist(B5M));

[d,a]=signature(B1B{1});
b11=moment(d,1);
b12=moment(d,2);
b13=moment(d,3);
b14=moment(d,4);
b15=moment(d,5);

[d,a]=signature(B2B{1});
b21=moment(d,1);
b22=moment(d,2);
b23=moment(d,3);
b24=moment(d,4);
b25=moment(d,5);

[d,a]=signature(B3B{1});
b31=moment(d,1);
b32=moment(d,2);
b33=moment(d,3);
b34=moment(d,4);
b35=moment(d,5);

[d,a]=signature(B4B{1});
b41=moment(d,1);
b42=moment(d,2);
b43=moment(d,3);
b44=moment(d,4);
b45=moment(d,5);

[d,a]=signature(B5B{1});
b51=moment(d,1);
b52=moment(d,2);
b53=moment(d,3);
b54=moment(d,4);
b55=moment(d,5);

% letter B pattern vector
%{
BPV=[blengthB1, distB1;
    blengthB2, distB2;
    blengthB3, distB3;
    blengthB4, distB4;
    blengthB5, distB5];
%}

BPV=[blengthB1, diamB1, b11, b12, b13, b14, b15;
     blengthB2, diamB2, b21, b22, b23, b24, b25;
     blengthB3, diamB3, b31, b32, b33, b34, b35;
     blengthB4, diamB4, b41, b42, b43, b44, b45;
     blengthB5, diamB5, b51, b52, b53, b54, b55];

% prepare pattern vectors for plotting
by=[distB1;distB2;distB3;distB4;distB5]; 
bx=[blengthB1;blengthB2;blengthB3;blengthB4;blengthB5];

%==========================================================================
%                 calculating the pattern vector for R   
%==========================================================================

R1=imcomplement(im2bw(r1,0.5));
R2=imcomplement(im2bw(r2,0.5));
R3=imcomplement(im2bw(r3,0.5));
R4=imcomplement(im2bw(r4,0.5));
R5=imcomplement(im2bw(r5,0.5));

R1B=bwboundaries(R1);
R2B=bwboundaries(R2);
R3B=bwboundaries(R3);
R4B=bwboundaries(R4);
R5B=bwboundaries(R5);

R1M=cell2mat(R1B);
R2M=cell2mat(R2B);
R3M=cell2mat(R3B);
R4M=cell2mat(R4B);
R5M=cell2mat(R5B);

[blengthR1 dim]=size(R1M);
[blengthR2 dim]=size(R2M);
[blengthR3 dim]=size(R3M);
[blengthR4 dim]=size(R4M);
[blengthR5 dim]=size(R5M);

diamR1=max(pdist(R1M));
diamR2=max(pdist(R2M));
diamR3=max(pdist(R3M));
diamR4=max(pdist(R4M));
diamR5=max(pdist(R5M));

[d,a]=signature(R1B{1});
r11=moment(d,1);
r21=moment(d,2);
r31=moment(d,3);
r41=moment(d,4);
r51=moment(d,5);

[d,a]=signature(R2B{1});
r12=moment(d,1);
r22=moment(d,2);
r32=moment(d,3);
r42=moment(d,4);
r52=moment(d,5);

[d,a]=signature(R3B{1});
r13=moment(d,1);
r23=moment(d,2);
r33=moment(d,3);
r43=moment(d,4);
r53=moment(d,5);

[d,a]=signature(R4B{1});
r14=moment(d,1);
r24=moment(d,2);
r34=moment(d,3);
r44=moment(d,4);
r54=moment(d,5);

[d,a]=signature(R5B{1});
r15=moment(d,1);
r25=moment(d,2);
r35=moment(d,3);
r45=moment(d,4);
r55=moment(d,5);

% letter R pattern vector

%{
RPV=[blengthR1, distR1;
    blengthR2, distR2;
    blengthR3, distR3;
    blengthR4, distR4;
    blengthR5, distR5];
%}
RPV=[blengthR1, diamR1, r11, r21, r31, r41, r51;
     blengthR2, diamR2, r12, r22, r32, r42, r52;
     blengthR3, diamR3, r13, r23, r33, r43, r53;
     blengthR4, diamR4, r14, r24, r34, r44, r54;
     blengthR5, diamR5, r15, r25, r35, r45, r55];
 
    
% rearrange pattern vector for plotting
ry=[distR1;distR2;distR3;distR4;distR5]; 
rx=[blengthR1;blengthR2;blengthR3;blengthR4;blengthR5];

%==========================================================================
%                         plot pattern vectors  
%==========================================================================

figure(5);
hold all; grid on;
plot(ax,ay,'or','MarkerSize',10,'LineWidth',1.25);
plot(bx,by,'xb','MarkerSize',10,'LineWidth',1.25);
plot(rx,ry,'*k','MarkerSize',10,'LineWidth',1.25);
title('Plots of pattern vectors for letters A, B, and R');
xlabel('Boundary length');
ylabel('Boundary diameter');
legend('A','B','R','location','NorthWest');

%==========================================================================
%              plot the signatures of the image boundaries 
%==========================================================================
figure(6);
for i=1:length(letters)
    bSq=bwboundaries(imcomplement(rgb2gray(letters{i})),'noholes');
    [distSq, angleSq]=signature(bSq{1});
    
    subplot(3,5,i);
    plot(angleSq, distSq)
    title(sprintf('Signature of %s', fileNames{i}));
    xlabel('Angular distance','FontSize',7);
    ylabel('Pixel distance from centroid','FontSize',7);
    
end

%==========================================================================
%                calculate the means of the pattern vectors  
%==========================================================================
MAPV=mean(APV);
MBPV=mean(BPV);
MRPV=mean(RPV);

%--------------------------------------------------------------------------
% calculate APV1 to all means
d11=APV(1,:)*(MAPV')-(1/2)*(MAPV)*(MAPV'); % A
d12=APV(1,:)*(MBPV')-(1/2)*(MBPV)*(MBPV'); % B
d13=APV(1,:)*(MRPV')-(1/2)*(MRPV)*(MRPV'); % R
[maxAPV1, indAPV1]=max([d11, d12, d13]);

% calculate APV2 to all means
d11=APV(2,:)*(MAPV')-(1/2)*(MAPV)*(MAPV');
d12=APV(2,:)*(MBPV')-(1/2)*(MBPV)*(MBPV');
d13=APV(2,:)*(MRPV')-(1/2)*(MRPV)*(MRPV');
[maxAPV2, indAPV2]=max([d11, d12, d13]);

% calculate APV3 to all means
d11=APV(3,:)*(MAPV')-(1/2)*(MAPV)*(MAPV');
d12=APV(3,:)*(MBPV')-(1/2)*(MBPV)*(MBPV');
d13=APV(3,:)*(MRPV')-(1/2)*(MRPV)*(MRPV');
[maxAPV3, indAPV3]=max([d11, d12, d13]);

% calculate APV4 to all means
d11=APV(4,:)*(MAPV')-(1/2)*(MAPV)*(MAPV');
d12=APV(4,:)*(MBPV')-(1/2)*(MBPV)*(MBPV');
d13=APV(4,:)*(MRPV')-(1/2)*(MRPV)*(MRPV');
[maxAPV4, indAPV4]=max([d11, d12, d13]);

% calculate APV5 to all means
d11=APV(5,:)*(MAPV')-(1/2)*(MAPV)*(MAPV');
d12=APV(5,:)*(MBPV')-(1/2)*(MBPV)*(MBPV');
d13=APV(5,:)*(MRPV')-(1/2)*(MRPV)*(MRPV');
[maxAPV5, indAPV5]=max([d11, d12, d13]);
%--------------------------------------------------------------------------

% calculate BPV1 to all means
d11=BPV(1,:)*(MAPV')-(1/2)*(MAPV)*(MAPV'); % A
d12=BPV(1,:)*(MBPV')-(1/2)*(MBPV)*(MBPV'); % B
d13=BPV(1,:)*(MRPV')-(1/2)*(MRPV)*(MRPV'); % R
[maxBPV1, indBPV1]=max([d11, d12, d13]);

% calculate BPV2 to all means
d11=BPV(2,:)*(MAPV')-(1/2)*(MAPV)*(MAPV');
d12=BPV(2,:)*(MBPV')-(1/2)*(MBPV)*(MBPV');
d13=BPV(2,:)*(MRPV')-(1/2)*(MRPV)*(MRPV');
[maxBPV2, indBPV2]=max([d11, d12, d13]);

% calculate BPV3 to all means
d11=BPV(3,:)*(MAPV')-(1/2)*(MAPV)*(MAPV');
d12=BPV(3,:)*(MBPV')-(1/2)*(MBPV)*(MBPV');
d13=BPV(3,:)*(MRPV')-(1/2)*(MRPV)*(MRPV');
[maxBPV3, indBPV3]=max([d11, d12, d13]);

% calculate BPV4 to all means
d11=BPV(4,:)*(MAPV')-(1/2)*(MAPV)*(MAPV');
d12=BPV(4,:)*(MBPV')-(1/2)*(MBPV)*(MBPV');
d13=BPV(4,:)*(MRPV')-(1/2)*(MRPV)*(MRPV');
[maxBPV4, indBPV4]=max([d11, d12, d13]);

% calculate BPV5 to all means
d11=BPV(5,:)*(MAPV')-(1/2)*(MAPV)*(MAPV');
d12=BPV(5,:)*(MBPV')-(1/2)*(MBPV)*(MBPV');
d13=BPV(5,:)*(MRPV')-(1/2)*(MRPV)*(MRPV');
[maxBPV5, indBPV5]=max([d11, d12, d13]);
%--------------------------------------------------------------------------

% calculate RPV1 to all means
d11=RPV(1,:)*(MAPV')-(1/2)*(MAPV)*(MAPV'); % A
d12=RPV(1,:)*(MBPV')-(1/2)*(MBPV)*(MBPV'); % B
d13=RPV(1,:)*(MRPV')-(1/2)*(MRPV)*(MRPV'); % R
[maxRPV1, indRPV1]=max([d11, d12, d13]);

% calculate RPV2 to all means
d11=RPV(2,:)*(MAPV')-(1/2)*(MAPV)*(MAPV');
d12=RPV(2,:)*(MBPV')-(1/2)*(MBPV)*(MBPV');
d13=RPV(2,:)*(MRPV')-(1/2)*(MRPV)*(MRPV');
[maxRPV2, indRPV2]=max([d11, d12, d13]);

% calculate RPV3 to all means
d11=RPV(3,:)*(MAPV')-(1/2)*(MAPV)*(MAPV');
d12=RPV(3,:)*(MBPV')-(1/2)*(MBPV)*(MBPV');
d13=RPV(3,:)*(MRPV')-(1/2)*(MRPV)*(MRPV');
[maxRPV3, indRPV3]=max([d11, d12, d13]);

% calculate RPV4 to all means
d11=RPV(4,:)*(MAPV')-(1/2)*(MAPV)*(MAPV');
d12=RPV(4,:)*(MBPV')-(1/2)*(MBPV)*(MBPV');
d13=RPV(4,:)*(MRPV')-(1/2)*(MRPV)*(MRPV');
[maxRPV4, indRPV4]=max([d11, d12, d13]);

% calculate RPV5 to all means
d11=RPV(5,:)*(MAPV')-(1/2)*(MAPV)*(MAPV');
d12=RPV(5,:)*(MBPV')-(1/2)*(MBPV)*(MBPV');
d13=RPV(5,:)*(MRPV')-(1/2)*(MRPV)*(MRPV');
[maxRPV5, indRPV5]=max([d11, d12, d13]);

