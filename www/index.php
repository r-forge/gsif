
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='https://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
    <style type="text/css">
<!--
.style1 {font-size: x-small}
.R_code {font-family:"Courier New", Courier, monospace;
font-style:italic;
font-size: x-small;
}
-->
    </style>
</head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="<?php echo $themeroot; ?>imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->

<div>
  <p>This package contains tools and procedures to handle soil data and produce gridded soil property maps to support the global soil data initivatives such as the GlobalSoilMap.net project. This package was developed as a support to the <a href="http://africasoils.net" target="_blank">Africa Soil Information Service project</a>. To access <strong><a href="http://soilgrids1km.isric.org">SoilGrids1km</a></strong> data sets and/or <a href="http://soilinfo.isric.org" target="_blank">download the App for mobile phones</a>, please refer to <a href="http://www.soilgrids.org"><strong>www.soilgrids.org</strong></a> </p>
</div>
<table border="0" cellspacing="0" cellpadding="10">
  <tr>
    <td><div align="center"><a href="http://globalsoilmap.net/content/1-degree-tiles-world"><img src="Fig_1degree_tiles_world.thumbnail.png" alt="One degree tiles (world map)" width="220" height="120" border="0" /></a></div></td>
    <td><a href="http://soilprofiles.org"><img src="opensoilprofiles.gif" alt="Open Soil Profiles" width="320" height="100" border="0" longdesc="http://soilprofiles.org" /></a></td>
  </tr>
</table>
<p class="style1">Contact: <a href="http://www.wewur.wur.nl/popups/vcard.aspx?id=HENGL001" target="_blank">Tomislav Hengl</a></p>
<p class="style1">Contributions by: Bas Kempen, Gerard B.M. Heuvelink, Dylan Beaudette, Reuter I. Hannes, Brendan Malone, Pierre Roudier ... </p>
<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. See the complete list of <strong><a href="00Index.html">functions</a></strong> available in this package. GSIF provides access to several case studies. A tutorial to analyze soil property and soil class data with the Ebergotzen data is available <strong><a href="http://gsif.isric.org/doku.php?id=wiki:tutorial_eberg">here</a></strong>. </p>
<p><iframe src="https://docs.google.com/presentation/d/1rhK4LZWbEVPX2mwgPdeKqSKvHjJOz-phQ5py1yj-1H0/embed?start=false&loop=false&delayms=3000" frameborder="0" width="533" height="429" allowfullscreen="true" mozallowfullscreen="true" webkitallowfullscreen="true"></iframe></p>
<p><strong>Installation:</strong></p>
<p>To install this package from R-forge use (works only on<strong> &gt;= R 2.15!</strong>):</p>
<p class="R_code">&gt; install.packages(c(&quot;RCurl&quot;, &quot;XML&quot;, &quot;rgdal&quot;, &quot;raster&quot;, &quot;sp&quot;, &quot;aqp&quot;,  &quot;mda&quot;, &quot;gstat&quot;, &quot;plotKML&quot;, &quot;dismo&quot;, &quot;rJava&quot;))<br />
&gt; install.packages(&quot;GSIF&quot;, repos=c(&quot;http://R-Forge.R-project.org&quot;), type = &quot;source&quot;)</p>
<p>GSIF package extensively uses a number of external software, hence it is highly recommended that, prior to starting GSIF, you first download and install:</p>
<ol>
  <li><a href="http://fwtools.maptools.org" target="_blank">FWTools</a> (<em>optional</em>) &#8212; this software is called by several functions (e.g. <a href="make.3Dgrid.html">make.3Dgrid</a>);</li>
  <li><a href="http://www.saga-gis.org" target="_blank">SAGA GIS</a> (<em>optional</em>) &#8212; this software is highly recommend but not required;</li>
</ol>
<p>Read more:  installation and first steps with <strong><a href="http://plotkml.r-forge.r-project.org/">plotKML</a></strong>.</p>
<p><strong>News:</strong></p>
<ul>
  <li>June 2014: added <a href="SoilGrids.html">SoilGrids class</a> (3D SpatialPixels with uncertainty); </li>
  <li>May 2014: added procedure to derive <a href="OCSKGM.html">soil organic carbon stock</a>; </li>
  <li>Apr 2014: new version of <a href="http://www.soilgrids.org/">SoilGrids1km</a> released (several new functions in the GSIF package included); </li>
  <li>Jun 2013: added a function to map uncertainty using the <a href="http://cran.r-project.org/web/packages/quantregForest/">quantregForest</a> package; </li>
  <li>Apr 2013: added functionality for map tiling (updated to sp package 1.0-8); </li>
  <li>Mar 2013: finished producing soil property maps of Africa at 1 km resolution; </li>
  <li>Feb 2013: GSIF package tutorial moved to the <a href="http://gsif.isric.org/doku.php?id=wiki:tutorial_eberg">DokuWiki page</a>; </li>
  <li>Dec 2012: added functionality for tiling large spatial objects; </li>
  <li>Nov 2012: revised gstatModel fitting functionality and reduced package dependencies; </li>
  <li>Aug 2012: added examples of how to <a href="http://plotkml.r-forge.r-project.org/tutorial.php">visualize</a> various (spatial) soil data; </li>
  <li>July 2012: added functionality for <a href="tutorial_eberg.php">3D regression-kriging and soil-class mapping</a>; </li>
  <li>Apr 2012: first meeting of the package development team at the <a href="http://www.pedometrics.org/dsm_oz/" target="_blank">DSM conference</a> in Sydney; </li>
  <li>Mar 2012: first version of the package on R-forge; </li>
</ul>
<hr />
<p><a href="http://www.isric.org" target="_blank"><img src="ISRIC_logo.jpg" alt="ISRIC - World Soil Information" width="363" height="98" border="0" longdesc="http://www.isric.org" /></a></p>
</body>
</html>
