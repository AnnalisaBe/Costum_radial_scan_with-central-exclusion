//Author: Annalisa Bellandi
//Faulkner lab, John Innes Centre

///This macro allows for:
// Radial scan in selected angle in files with two colour channels
//PART1: registration - channel 1 is adjusted to account for movement of the sample, the adjustment matrix is then used on channel 2
//PART2: process each channel of the registered images by using metadata of original image


/// Main steps:
/// The user selects three points. Using those 3 points a circle is created with the points falling in on the circumference.
/// one line is created starting from point 1 to the centre of the circumference, and a second line is created from point 2 to the centre
/// of the circumference. 
/// Rays are created in the angle between line 1 and line 2, at a distance called step that can be set up
/// fluorescence intensity along the rays is saved for each time step


/// Information on the circle drawing here: http://imagej.1557.x6.nabble.com/macro-to-draw-a-circle-from-three-points-td5015618.html


//=======================================================================================================


///--------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------- PART 1: Run Hyperstackreg and register the series

folder = getDirectory("Choose folder with images to be registered")

//this is automated batch processing
setBatchMode(true);

list = getFileList(folder);
		for (i=0; i<list.length; i++) {
            print(list[i]);
            path = folder+list[i];
            run("Bio-Formats Importer", "open=[" + path + "] autoscale color_mode=Default open_files rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			imageTitle=getTitle;
			title_without_extension = substring(imageTitle, 0, lengthOf(imageTitle)-4);

run("HyperStackReg ", "transformation=[Rigid Body] channel1 show");

//make sure image selected is the output file of the plugin
selectWindow(title_without_extension + "-registered.czi");

//this converts the composite merged image into a normal hyperstack where the channels are visible separately
run("Channels Tool...");
Stack.setDisplayMode("grayscale");

saveAs("Tiff",  folder + title_without_extension + "_reg.tif");
close();

}



///----------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------     PART 2
									
//open not registered image
//note slap frame
//note stable frame
											
dir = getInfo("image.directory");// get the image directory
print(dir)
imageTitle=getTitle; //get the image name
											
//get pixel size
getPixelSize(unit, pixelWidth, pixelHeight);
print("Current image pixel width = " + pixelWidth + " " + unit +".");
pxw = pixelWidth;
print(pxw);
											
//check planeTimings.txt https://docs.openmicroscopy.org/bio-formats/5.7.3/users/imagej/
run("Bio-Formats Macro Extensions"); //to get metadata info 
selectWindow(imageTitle); //have to return to the original because I didn't save duplicate in folder
id = getInfo("image.directory") + getInfo("image.filename");
Ext.setId(id);
Ext.getImageCount(imageCount);
deltaT = newArray(nSlices); //create array with vector of time points
//memorize the title of the image
title_without_extension = substring(imageTitle, 0, lengthOf(imageTitle)-4);
											
close(); //close orignal
									
//import the registered version of the image
run("Bio-Formats Importer", "open=[" + dir + title_without_extension + "_reg.tif" + "] autoscale color_mode=Default open_files rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
											
									
									
//select image and duplicate to make measures
dir = getInfo("image.directory");// get the image directory
print(dir)
imageTitle=getTitle; //get the image name
run("Duplicate...", "duplicate"); //duplicate the entire stack for measurements
measuredImageTitle="dup"+imageTitle; //add the starting "dup" to the original image title
rename(measuredImageTitle);//rename it
selectWindow(measuredImageTitle); //select the duplicate image to measure
						
waitForUser("make 3 points on the edge of the trichome base");
									
//draw the circle
if ( selectionType() == 10 ){ // is MultiPoint selection
	getSelectionCoordinates(x, y);
									
	if (lengthOf(x) != 3){
	print("Must be 3 points !");
	return;
	}
									
	b1 = x[0]* x[0] + y[0]*y[0];
	b2 = x[1]* x[1] + y[1]*y[1];
	b3 = x[2]* x[2] + y[2]*y[2];
									
	xL1 = x[0];
	yL1 = y[0];
									
	xL2 = x[1];
	yL2 = y[1];
									
	detA = 4*( x[0]*(y[1]-y[2]) - y[0]*(x[1]-x[2]) + (x[1]*y[2]-x[2]*y[1]) );
									
	xc = 2*( b1*(y[1]-y[2]) - y[0]*(b2-b3) + (b2*y[2]-b3*y[1]) )/detA;
	yc = 2*( x[0]*(b2-b3) - b1*(x[1]-x[2]) + (x[1]*b3-x[2]*b2) )/detA;
	K = 4*( x[0]*(y[1]*b3-y[2]*b2) - y[0]*(x[1]*b3-x[2]*b2) + b1*(x[1]*y[2]-x[2]*y[1]) )/detA;
	r = sqrt(abs(K + xc*xc + yc*yc));
									
	print("xc : " + xc);
	print("yc : " + yc);
	print("radius of base [pixels number]: " + r);
	r_um = r*pxw;
	print("radius of base [um]: " + r_um);
									
	run("Split Channels");
	selectWindow("C1-"+measuredImageTitle);
									
	makeOval(xc-r, yc-r, 2*r, 2*r);
	Overlay.addSelection("magenta", 3);
									
	//line 1
	makeLine(xc, yc, xL1, yL1, 3);
									 
	HL1 = yL1-yc;
	LL1 = xL1-xc;
									
	alfaL1 = atan2(HL1,LL1);// the formulas take radiants so no need to convert in degrees
	Overlay.addSelection("yellow", 3);
										
	//line 2
	makeLine(xc, yc, xL2, yL2, 3);
									 
	HL2 = yL2-yc;
	LL2 = xL2-xc;
									
	alfaL2 = atan2(HL2,LL2);// the formulas take radiants so no need to convert in degrees
	Overlay.addSelection("yellow", 3);
									
	scan = 700;
	print("green lines lenght [um]:" + (scan+r_um));
	print("green lines - radius [um]:" + scan);
	len=scan/pxw;
	linelenght=len+r; //in pixels because the x and y coordinatesd are in pixels
	imin = alfaL1;
	imax = alfaL2;
	stepsize=2*PI/200; //define the distance from one ray to the next
									
	run("Clear Results");
			
selectWindow("C1-"+measuredImageTitle);
geno="iglu";
									
//Loop0: repeats the loop1 and 2 for each time point (slice) - the loops starts like this: for(k=0; k<=nSlices; k++){ , insert slice in place of 0
for(k=1; k<=nSlices; k++){
	setSlice(k);

	no=k-1; // arrays start at zero, so I subtract 1.
	Ext.getPlaneTimingDeltaT(deltaT[no], no); //....
	t=deltaT[no];
	print(t);
									
	nstep= abs(imax-imin)/stepsize;
	print("nstep:" +nstep);
	print(abs(imax-imin));
									
	for(i=0; i<=nstep; i+=1){
		alfa = i*stepsize+imin;
		makeLine(xc,yc,xc+linelenght*cos(alfa),yc+linelenght*sin(alfa));
		Overlay.addSelection("green", 2);
		profile = getProfile();
									
		//Loop2: for each  profile, saves the values of the profile and the position of the values along the lenght in a table
		for(j=0; j<profile.length; j++){
			setResult(i, j+linelenght*(k-1), profile[j]); //for each angle saves the profile along the lenght of the line
			setResult("frame", j+linelenght*(k-1), k);
			setResult("t", j+linelenght*(k-1), t); 
			setResult("d", j+linelenght*(k-1), j*pxw); // this should give the space in um in my d column instead fo the pxl n
			setResult("r drop", j+linelenght*(k-1), r_um); // saves the lenght of radius of the droplet in micron
			setResult("n radii", j+linelenght*(k-1), nstep); //saves the number of radii that fit in the scanned angle between line 1 and line 2
									        }
									    }
									}

title_without_extension = substring(imageTitle, 7, lengthOf(imageTitle)-8);
saveAs("Results", dir + geno + title_without_extension + ".csv");
run("Clear Results");

selectWindow("C2-"+measuredImageTitle);
geno="rgeco";


makeOval(xc-r, yc-r, 2*r, 2*r);
Overlay.addSelection("magenta", 3);
//line 1
makeLine(xc, yc, xL1, yL1, 3);
									 
HL1 = yL1-yc;
LL1 = xL1-xc;
									
alfaL1 = atan2(HL1,LL1);// the formulas take radiants so no need to convert in degrees
//print("angle line 1:" + alfaL1);
Overlay.addSelection("yellow", 3);
										
//line 2
makeLine(xc, yc, xL2, yL2, 3);
									 
HL2 = yL2-yc;
LL2 = xL2-xc;
									
alfaL2 = atan2(HL2,LL2);// the formulas take radiants so no need to convert in degrees
Overlay.addSelection("yellow", 3);
									
scan = 700;
print("green lines lenght [um]:" + (scan+r_um));
print("green lines - radius [um]:" + scan);
len=scan/pxw;
linelenght=len+r; //in pixels because the x and y coordinatesd are in pixels
imin = alfaL1;
imax = alfaL2;
stepsize=2*PI/200; //define the distance from one ray to the next
			
//Loop0: repeats the loop1 and 2 for each time point (slice) - the loops starts like this: for(k=0; k<=nSlices; k++){ , insert slice in place of 0
	for(k=1; k<=nSlices; k++){
		setSlice(k);

		no=k-1; // arrays start at zero, so I subtract 1.
		Ext.getPlaneTimingDeltaT(deltaT[no], no); //....
		t=deltaT[no];
		print(t);
									
		nstep= abs(imax-imin)/stepsize;
		print("nstep:" +nstep);
		print(abs(imax-imin));
									
		for(i=0; i<=nstep; i+=1){
			alfa = i*stepsize+imin;
			makeLine(xc,yc,xc+linelenght*cos(alfa),yc+linelenght*sin(alfa));
			Overlay.addSelection("green", 2);
			profile = getProfile();
									
//Loop2: for each  profile, saves the values of the profile and the position of the values along the lenght in a table
			for(j=0; j<profile.length; j++){
				setResult(i, j+linelenght*(k-1), profile[j]); //for each angle saves the profile along the lenght of the line
				setResult("frame", j+linelenght*(k-1), k);
				setResult("t", j+linelenght*(k-1), t); 
				setResult("d", j+linelenght*(k-1), j*pxw); // this should give the space in um in my d column instead fo the pxl n
				setResult("r drop", j+linelenght*(k-1), r_um); // saves the lenght of radius of the droplet in micron
				setResult("n radii", j+linelenght*(k-1), nstep); //saves the number of radii that fit in the scanned angle between line 1 and line 2
									        }
									    }
									}
									 
title_without_extension = substring(imageTitle, 7, lengthOf(imageTitle)-8);
saveAs("Results", dir + geno + title_without_extension + ".csv");		
run("Clear Results");

waitForUser("Happy with the scan?");
close(); 
close(); 