%this code reproduces figure 2 in Field and Wood, mean  water vapor path vs
%mean wind speed for cyclones in the north atlantic

load('erai_MCMS_NH_special_prw_datacyc_2001_boxlen_degrees_45.mat')

% (1) GRAB THE NECESSARY DATA
[datamat lonmat latmat]=datacyc_data_to_mat(datacyc_2001);                                         
[fulllon fulltt fullcc]=func_datacyc_var_to_mat(datacyc_2001,'fulllon');                           
[fulllat fulltt fullcc]=func_datacyc_var_to_mat(datacyc_2001,'fulllat');     

% (2) SET THE BOX REGION
lonW=310; 
lonE=350;
latN=60;
latS=30;
% (3) FIND THE POINTS IN THE BOX REGION
ilon=find_between(fulllon,lonW,lonE);
sublat=fulllat(ilon);
isublat=find_between(sublat,latS,latN);
% (4) HERE IS Index of cyclones we want numbered by their position in the
%datamat matrix.
iuse=ilon(isublat);

%find which cyclones we are interested in
cycloneindex=unique(fulltt(iuse));

%create a dictionary of those cyclones

cyclonedictionary=datacyc_2001(cycloneindex);

[newdatamat lonmat latmat]=datacyc_data_to_mat(cyclonedictionary);

%average data across time, lat and lon
prwaverage=nanmean(nanmean(newdatamat,2),3); 





%now do again with windspeed

load('erai_MCMS_NH_special_ws_datacyc_2001_boxlen_degrees_45.mat')

% (1) GRAB THE NECESSARY DATA
[wsdatamat lonmat latmat]=datacyc_data_to_mat(datacyc_2001);                                         
[fulllon fulltt fullcc]=func_datacyc_var_to_mat(datacyc_2001,'fulllon');                           
[fulllat fulltt fullcc]=func_datacyc_var_to_mat(datacyc_2001,'fulllat');     

% (2) SET THE BOX REGION
lonW=310; 
lonE=350;
latN=60;
latS=30;
% (3) FIND THE POINTS IN THE BOX REGION
ilon=find_between(fulllon,lonW,lonE);
sublat=fulllat(ilon);
isublat=find_between(sublat,latS,latN);
% (4) HERE IS Index of cyclones we want numbered by their position in the
%datamat matrix.
iuse=ilon(isublat);

%find which cyclones we are interested in
wscycloneindex=unique(fulltt(iuse));

%create a dictionary of those cyclones

cyclonedictionary=datacyc_2001(wscycloneindex);

[newdatamat lonmat latmat]=datacyc_data_to_mat(cyclonedictionary);

%average data across time, lat and lon
wsaverage=nanmean(nanmean(newdatamat,2),3); 





%%make a scatter plot
figure
plot(wsaverage,prwaverage,'k*')

title('Mean Water Vapor Path vs Mean Wind speed')