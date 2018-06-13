clc;
clear all;
close all
sift1 = importdata('./obj02_001.png.harhes.sift',' ',2); %x,y,a,b,c,desc
sift2 = importdata('./obj02_002.png.harhes.sift',' ',2);

desc1 = sift1.data(:,6:133);
coord1 = sift1.data(:,1:2);
coord2 = sift2.data(:,1:2);
desc2 = sift2.data(:,6:133);

[matches, scores] = vl_ubcmatch(desc1' , desc2',5);
nrc=99; %number of correspondences 

matchedcoord1=coord1(matches(1,:),:);
matchedcoord2=coord2(matches(2,:),:);

Dist=(matchedcoord1(:,1).^2+matchedcoord1(:,2).^2).^0.5;
[C I]=sort(Dist); %I=matchindices sorted by distance to origin from first set of coordinates

id=floor(linspace(1,size(I,1),nrc));
Isel=I(floor(linspace(1,size(I,1),nrc)));
for i=1:size(Isel,1)
    if (Dist(Isel(i))/mean(Dist))>1.5
        if i==1
            Isel(i)=I(id(i+1));
        else
            Isel(i)=I(id(i-1));
        end
    end
end

matches_selected=matches(:,Isel);
p1=coord1(matches_selected(1,:)',:);
p2=coord2(matches_selected(2,:)',:);

%Normalizing
%for image 1
mx1=sum(p1(:,1))/size(p1,1);
my1=sum(p1(:,2))/size(p1,1);
d1=sum(sqrt((p1(:,1)-mx1).^2+(p1(:,2)-my1).^2))/size(p1,1);
T1=[(sqrt(2))/d1 0 -mx1*sqrt(2)/d1; 0 (sqrt(2))/d1  -my1*sqrt(2)/d1; 0 0 1];

%for image 2
mx2=sum(p2(:,1))/size(p2,1);
my2=sum(p2(:,2))/size(p2,1);
d2=sum(sqrt((p2(:,1)-mx2).^2+(p2(:,2)-my2).^2))/size(p2,1);
T2=[(sqrt(2))/d2 0 -mx2*sqrt(2)/d2; 0 (sqrt(2))/d2  -my2*sqrt(2)/d2; 0 0 1];

%prepare coordinates for transformation
p1_norm=[p1 ones(size(p1,1),1)]'; 
p2_norm=[p2 ones(size(p2,1),1)]';

for i=1:size(p1_norm,2)
    p1_norm(:,i)=T1*p1_norm(:,i);
    p2_norm(:,i)=T2*p2_norm(:,i);
end

p1_norm=p1_norm';
p2_norm=p2_norm';

A_norm = [p1_norm(:,1).*p2_norm(:,1) p1_norm(:,1).*p2_norm(:,2) p1_norm(:,1) p1_norm(:,2).*p2_norm(:,1) p1_norm(:,2).*p2_norm(:,2) p1_norm(:,2) p2_norm(:,1) p2_norm(:,2) ones(size(matches_selected,2),1)]; 
[U,D,V] = svd(A_norm);

f_norm=V(:,end);
F_norm=reshape(f_norm,3,3);

% Decompose F
[U_f,D_f,V_f] = svd(F_norm);
D_f(end, end) = 0;
F_norm = U_f*D_f*V_f.'; 

p1=[p1 ones(size(p1,1),1)];
p2=[p2 ones(size(p2,1),1)];

mostInliers = 0;
N = 1000;
for n = 1:N
    % Get 8 random points
    perm = randperm(length(p1_norm));
    P    = 8;
    seed = perm(1:P);
    p1_ransac = p1_norm(seed,:);
    p2_ransac = p2_norm(seed,:);
    
    A_ransac = [p1_ransac(:,1).*p2_ransac(:,1) p1_ransac(:,1).*p2_ransac(:,2) p1_ransac(:,1) p1_ransac(:,2).*p2_ransac(:,1) p1_ransac(:,2).*p2_ransac(:,2) p1_ransac(:,2) p2_ransac(:,1) p2_ransac(:,2) ones(P,1)]; 
    [U,D,V] = svd(A_ransac);
    f_ransac=V(:,end);
    F_ransac=reshape(f_ransac,3,3);

    % Decompose F and enforce singularity
    [U_f,D_f,V_f] = svd(F_ransac);
    D_f(end, end) = 0;
    F_ransac = U_f*D_f*V_f.'; 
    
    %calculate the distance for each point
    d = zeros(size(p1_norm,1),1);
    for i=1:size(p1_norm,1)
        term1 = (p2_norm(i,:)*F_ransac*p1_norm(i,:)')^2;
        term2 = (F_ransac(1,:)*p1_norm(i,:)')^2;
        term3 = (F_ransac(2,:)*p1_norm(i,:)')^2;
        term4 = (F_ransac(:,1)'*p1_norm(i,:)')^2;
        term5 = (F_ransac(:,2)'*p1_norm(i,:)')^2;
        d(i) = (term1)/( term2 + term3 + term4 + term5);
    end
    treshold = 1e-5;
    inliers = length(find(d < treshold));
    disp(d)
    if inliers > mostInliers
        mostInliers = inliers
        F_ransac_best = F_ransac;
    end
end

%Denormalization of ransac 
F1=T1'*F_ransac_best*T1;
F2=T2'*F_ransac_best*T2;
F_ransac=F1;

%denormalize the 'normal' matrices and return as sanity check
F1_denorm=T1'*F_norm*T1;
F2_denorm=T2'*F_norm*T2;
F_denorm = F1_denorm

p1=p1(:,1:2);
p2=p2(:,1:2);

%plotting, plot the F_ransac

I1=imread('TeddyBearPNG/obj02_001.png');
I2=imread('TeddyBearPNG/obj02_002.png');

figure;
subplot(121)
imshow(I1)
title('Inliers and Epipolar Lines in First Image RANSAC'); hold on;
plot(p1(:,1),p2(:,2),'go');

lines1=epipolarLine(F_ransac',p2); %Ax+By+C
epipoint1=lineToBorderPoints(lines1,size(I1));
line(epipoint1(:,[1,3])',epipoint1(:,[2,4])');

subplot(122); 
imshow(I2)
title('Epipolar lines in second image RANSAC'); hold on; 
plot(p2(:,1),p2(:,2),'go');  

lines2=epipolarLine(F_ransac,p1);
epipoint2=lineToBorderPoints(lines2,size(I2));
line(epipoint2(:,[1,3])',epipoint2(:,[2,4])');
truesize; 

%plot the figure with 'normal F'
figure;
subplot(121)
imshow(I1)
title('Inliers and Epipolar Lines in First Image F_normal'); hold on;
plot(p1(:,1),p2(:,2),'go');

lines1=epipolarLine(F_denorm',p2); %Ax+By+C
epipoint1=lineToBorderPoints(lines1,size(I1));
line(epipoint1(:,[1,3])',epipoint1(:,[2,4])');

subplot(122); 
imshow(I2)
title('Epipolar lines in second image v'); hold on; 
plot(p2(:,1),p2(:,2),'go');  

lines2=epipolarLine(F_denorm,p1);
epipoint2=lineToBorderPoints(lines2,size(I2));
line(epipoint2(:,[1,3])',epipoint2(:,[2,4])');
truesize; 



