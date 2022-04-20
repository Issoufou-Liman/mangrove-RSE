<center> <h1>This repository contains a reproducible workflow for spatio-temporal analysis of mangroves vegetation in Ghana based on LANDSAT time series in Google Earth Engine.</h1> </center>

<div class="figure" style="text-align: left">
<img src="./Mangrove_Ghana.png" alt="a plot of mtcars"  />
<p class="caption"><i><b><font color="blue">A Google Earth view of mangroves in Ghana</font></b></i></p>
</div>

Mangroves are multi-purpose trees and shrubs that develop adaptation mechanisms allowing them to survive in saline water. They provide many ecosystems service to coastal communities all over the world. Evidences have shown that the coverage of these special type of forests have recently been increasing in some areas while decreasing in others. Here, we used time series of remotely-sensed data from Landsat to study the dynamics of mangroves vegetation in Ghana.

This repository contains Google Earth Engine and R codes for a spatio-temporal analysis of mangrove vegetation for Ghana

To replicate this work, you need to make sure that you have a working installation of `Google Earth Engine` `python API` for use in `rgee` package. Then, look for the following line of code under the chunk `Setting-the-Scene` to update you credentials.

```
 ee_Initialize(user = 'your user name', drive = TRUE, gcs = FALSE)
           
```
