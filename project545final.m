close all; clear all; clc
%PBIO 545 Project: Projecting hyperspectral data onto cone space and modified cone
%space and PCA of hyperspectral images
%How does the overlap between the L and M spectral sensitivity affect
%correlation between L and M cone activity? What happens when we change the
%sensitivity spectra of the L cone? Can this result be explained by PCA of
%the hyperspectral data? 

%Part I. Projecting our hyperspectral data onto cone space and modified cone 
%space. 
%Found data on http://www.cvrl.org/cones.htm from Stockman and Sharpe, 2000
%on the cone spectral sensitivity functions. I used the linear energy model
%provided on the website. 
%First we pull out the cone spectral sensitivity functions we already have from
%csv file and separating it into separate variables. 
Conedata=csvread('linss2_10e_5.csv'); 
wavelength=Conedata(:,1); 
Lcone=Conedata(:,2); 
Mcone=Conedata(:,3); 
Scone=Conedata(:,4); 

%Let's graph this information. As you can see, there is a lot of overlap
%between the L and M cone spectral sensitivity indicating similar firing rates 
%between these cones at middle wavelengths.The S cone sensitivity is in blue, 
%the M cone is in green, and the L cone is in red. 
figure(1); 
plot(wavelength, Scone); 
hold on; 
plot(wavelength, Mcone, 'g'); 
hold on; 
plot(wavelength, Lcone, 'r'); 
xlabel('Wavelength (nm)'); 
ylabel('Cone spectral sensitivity'); 
title('Cone spectral Sensitivity'); 


%Now we need to load in the hyperspectral data
%This data comes from http://personalpages.manchester.ac.uk/staff/david.foster/Hyperspectral_images_of_natural_scenes_04.html
%which has figures from Foster et al., 2006. I used image 2. 
%This hyperspectral data has three dimensions: pixels in x and y
%directions and intensity across spectra in the third dimension. The spectra 
%or wavelengths vary from 400 to 720 in steps of 10 nm. 
%The data is all saved in a variable called reflectances.
load ref_ruivaes1bb_reg1

%To pull out individual pixels of reflectances, we can make a for loop. I
%chose a random row of pixels to analyze in the image. 
Spectraldata=zeros(200, 33); 
for i=1:1017
    Spectraldata(i,:)=reflectances(i,100,:); 
end 


%
%Let's transform our data into cone space! We will transform into L, S and
%M cone space. We can do this by projecting our hyperspectral data
%into the L space and the M space as described by our Sharpe and Stockman data.
%This is done using a dot product. 

 
%Calculate the norm of each sensitivity of the cones so you can project
%onto the unit vector of each of these spaces.
Lnorm=norm(Lcone); 
Snorm=norm(Scone); 
Mnorm=norm(Mcone); 
%Projecting onto cone space
for i=1:1017 
    Lconeproj(i)=Spectraldata(i,:)*(Lcone/Lnorm); 
    Mconeproj(i)=Spectraldata(i,:)*(Mcone/Mnorm); 
    Sconeproj(i)=Spectraldata(i,:)*(Scone/Snorm); 
end 

%Let's graph our results so we can visualize our hyperspectral data in L
%and M cone space. 
figure(2); 
plot3(Lconeproj, Mconeproj, Sconeproj, 'or'); 
xlabel('L cone Projection'); 
ylabel('M cone Projection'); 
zlabel('S cone Projection'); 
grid on 
rotate3d on
title('Projections of Spectral data in cone space'); 

%If you look at this graph, the L and M cone activation by this image looks
%correlated. Let's calculate the correlational coefficients. 

rLM=corrcoef(Lconeproj, Mconeproj); 
rMS=corrcoef(Mconeproj, Sconeproj); 
rSL=corrcoef(Sconeproj, Lconeproj); 

%We can see that the correlational coefficient was .9985, indicating a very
%strong positive correlation between activity from L and M cone space. 

%What happens if we shift over the L cone spectral sensitivity so it peaks 
%at a larger wavelength? Will that decrease the correlation between L and M
%activity? 

%Make new L cone spectral sensitivity data by shifting our data over by 10.
newLcone=zeros(33,1); 
for i=11:33
    newLcone(i)=Lcone(i-10);
end 


%Let's graph this new cone spectra. We can see that there is now much less 
%overlap between the L and M cone spectra. 
%Note: the information in the hyperspectral images so our spectral sensitivities 
%also only go that high. This leads to the L cone spectra being cut off in 
%this graph. 
figure(3); 
plot(wavelength, Scone); 
hold on; 
plot(wavelength, Mcone, 'g'); 
hold on; 
plot(wavelength, newLcone, 'r'); 
xlabel('Wavelength (nm)'); 
ylabel('Cone spectral sensitivity'); 
title('Cone spectral Sensitivity'); 

%Let's calculate the new projection onto the new L cone axis. We do not
%need to project the hyperspectral data again onto the M and S cone axes 
%because they were not changed. 
for i=1:1017
    Lconeprojnew(i)=Spectraldata(i,:)*(newLcone/norm(newLcone)); 
end 

%Let's graph the hyperspectral data projected into this new cone space.
%Again we can see that there appears to be correlation between L and M cone
%responses. 
figure(4); 
plot3(Lconeprojnew, Mconeproj, Sconeproj, 'or'); 
xlabel('L cone Projection'); 
ylabel('M cone Projection'); 
zlabel('S cone Projection'); 
grid on 
rotate3d on
title('Projections of Spectral Data in Modified Cone Space'); 

%Let's calculate the new correlational coefficients. We do not need to
%recalculate the correlation between M and S projections because these
%projections did not change with the current manipulation. 
rnewLM=corrcoef(Lconeprojnew, Mconeproj); 
rnewSL=corrcoef(Sconeproj, Lconeprojnew);

%We notice that the correlational coefficient for L and M projections was 
%slightly lower but still quite high at are still quite high at .9601 despite
%the shift over in the spectral sensitivity of the L cones. Why might this 
%be. 

%%
%Part 2: PCA of hyperspectral data 
%Can the PCA or principle component analysis of the hyperspectral data
%explain the high correlation seen between the modified L cone space and
%the M cone space? 
 
%We begin by plotting the  mean and wavelength-dependent variance of our
%hyperspectral data. If we look at the mean of the spectral data in a graph
%form, we can see that there are peaks in the peak near 550 nm and at the
%end of the wavelength spectra we have information for, around 700 nm. We
%see a similar shape in the variance graph. The 550nm small peak lines up well with
%the spectral sensivitiy of our M cones and 700 nm peak lines up well with
%the modified L cone spectral sensivity. However, we do not know if the
%peaks near 550 nm and 700 nm are changing together
figure(1);
clf;
subplot(1, 2, 1)
plot(wavelength, mean(Spectraldata))
xlabel('Wavelength (nm)')
ylabel('Intensity')
title('Average Spectral Intensity'); 
axis tight
subplot(1, 2, 2)
plot(wavelength, var(Spectraldata))
xlabel('Wavelength(nm)')
ylabel('variance (intensity)')
title('Variance in Spectral Intensity'); 
axis tight 

% Blue is the mean and variance of the spectral data for the first row of
% pixels in image 2. If we look at the variance and the mean plots, we can
% see that the variance looks similar to the means accross these pixels. 
% Look at the variance plot. We see the variance at each wavelength. However, 
% we do not have information about how variance between wavelengths are related
% to each other.   
% First step to computing a PCA is computing the covariance matrix of my 
% row of hyperspectral data with the mean subtracted.
SpectraldataCovar = cov(Spectraldata - repmat(mean(Spectraldata), 1017, 1));

%Now we can start breaking the data into principle components. The first
%thing we want to do is examine the variance accross our components to see
%how much of the varaince is accounted for by each of the components. This 
%will tell us which components are best able to explain our data.
SpectraldataEigVal = eig(SpectraldataCovar)

%Let's graph this so we better visualize the size of each of the
%eigenvalues which represent the variance found in my hyperspectral data.
%We can see that the first component is quite large and the second
%component may be significant but there is a significant drop off in
%eigenvalue size after these first two components.
figure(2)
clf
plot(SpectraldataEigVal, 'o')
xlabel('component')
ylabel('variance')
title('Variance accounted by EigenVectors')

%Let's see what percentage of the total variance the fist two components of
%the PCA account for. 
%I summed the total variance by summing the eigenvalues and taking the
%proportion of that sum for the first and second eigen values of that
%whole.
eigvalues=diag(EigVal);
totalVar=sum(eigvalues); 
percentlargest=eigvalues(33)/totalVar; 
percentsecondlargest=eigvalues(32)/totalVar; 
percentfirstandsecond=percentlargest+percentsecondlargest;
%We can see that the first PCA component accounts for 92.52 percent of the
%data. With the second component with the first, we can account for 96.56 
%percent.

%Now let's actually compute those principle directions by completing an 
%eigenvalue decomposition. Once we do that, we can project the Spectraldata
%onto the eigenvectors and calculate the variance. 

[EigVec, EigVal] = eig(SpectraldataCovar);
VarSpectralData = var(Spectraldata * EigVec);

% To double check the variance of the data matches the spread accross eigenvectors 
% we plot variance of SpectralData onto the eigenvectors, or possible
% principle components. We see the same spread of variance we saw earlier
% for the eigenvalues. 

figure(1)
clf
plot([1:length(VarSpectralData)], VarSpectralData, 'or')
xlabel('component')
ylabel('variance')

% We can again see from this graph that most of the structure lies in two 
% of the eigenvectors. Let's do an eigenvalue decomposition that only uses 
% the two largest components. We can then examine these two components and 
% see what we learn about the structure of the data. 

% We begin by getting a few selected eigenvectors based on our eigenvalue 
% analysis. (by default eigs gives us those with the largest eigenvalues 
% when we only ask for two):

[EigVec, EigVal] = eigs(SpectraldataCovar, 2);

% Let's look at these. The blue line refers to the first PCA component that 
% accounts for most of the data and the green line refers to the second PCA 
% component that accounts for much less of the data. 
figure(1)
clf
plot(wavelength,EigVec);
xlabel('wavelength'); 
ylabel('intensity'); 
title('PCA Components'); 
%We can see these components look different from each other only at larger
%wavelengths. However, we also see bumps in both of these components at
%around 550 nm and 700 nm. In particular, if we look at the largest
%component we can see that there is a negative bump around 550 and 700.
%This suggests that this component accounts for changing in both of these
%areas at the same time. This may indicate that a shifted over L cone
%projection will still be very correlated with the M cone projection
%because even though the L cone projection that is shifted over can pick up
%more variance at higher wavlengths, when the L cone projection is picking
%up more variance and firing more, so is the M cone projection. Therefore,
%even though the new L cone spectra is picking up new information, its projection
%is still correlated with M cone production just due to the row in the
%hyperspectral data I chose. 

%Therefore, the new L cone is not getting access to independent information.
%Thus, there is not much pressure in this context to push the L cone
%spectra to sensitivity to larger wavelengths because the pieces of the 
%spectrum are correlated so you are not getting new information. 

%If we look at this image, we can see that the image is mostly green with 
%differing levels of brightness. Therefore, the spectral properties are likely 
%very similar across a single row of the image. 

%We can look at the spectral data projected onto the two eigenvectors if we
%would like, which confirms that much of the variance in the data is on the
%axis accounting for the first component. 
figure(2)
clf
plot(Spectraldata * EigVec(:, 1), Spectraldata * EigVec(:, 2), 'or')
hold on
xlabel('component 1')
ylabel('component 2')
grid on

%To better examine this question, one would analyze more of this image and
%other hyperspectral images containing different information or pieces of
%the world, such as the sky. 