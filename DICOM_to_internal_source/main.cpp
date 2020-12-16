#include "DICOM.h"

double interp(double x, double x1, double x2, double y1, double y2) {
	return (y2*(x-x1)+y1*(x2-x))/(x2-x1);
}

void submerge(QVector <DICOM *> &data, int i, int c, int f) {
	// We have three indices, l for one subsection, r the other, and j for the new sorted array
	int l = i, r = c+1, j = 0;
	QVector <DICOM *> temp(f-i+1); // Set aside memory for sorted array
	
	// While we have yet to iterate through either subsection
	while (l <= c && r <= f) {
		// If value at r index is smaller then add it to temp and move to next r
		if (data[l]->z > data[r]->z)
			temp[j++] = data[r++];
		// If value at l index is smaller then add it to temp and move to next l
		else
			temp[j++] = data[l++];
	}
	
	// Add all the remaining elements
	while (r <= f)
		temp[j++] = data[r++];
	while (l <= c)
		temp[j++] = data[l++];

	// Reassign all the data values to the temp values
	for (int k = 0; k < j; k++) {
		data[i+k] = temp[k]; 
	}
}

void mergeSort(QVector <DICOM *> &data, int n) {
	// If our array is size 1 or less quit
	if (n <= 1)
		return;
	
	int subn = 1; // subn the size of subsections that are being submerged
	int i = 0; // i is the current index of the array at which we are at
	
	// While we are still submerging sections that are smaller than the size of the array
	while (subn < n) {
		// Reset the index to 0
		i = 0;
		
		// Iterate through n/(2*subn) sections, truncated
		while (i < n - subn) {
			
			// submerge two subn sized portions of data of the array
			if (i + (2 * subn) < n)
				submerge (data, i, i + subn - 1, i + 2 * subn - 1);
			
			// Or submerge a subn sized section and whatever is left of the array
			else 
				submerge (data, i, i + subn - 1, n - 1);
			
			// Move the index to submerge the next 2 subsections
			i += 2 * subn; 
		}
		
		// Double the size of subsection to be merged
		subn *= 2; 
	}
}

int main(int argc, char **argv) {
	// Start clock for timing
    std::clock_t start;
    double duration;
    start = std::clock();
	
	// ---------------------------------------------------------- //
	// PARSE INPUT                                                //
	// ---------------------------------------------------------- //
	/*
	In this section, all the inputs that the programmed is invoked
	with are parsed.  This section reads in any file name given,
	assuming its a DICOM file, unless the input has the format
	"-outputImages" or ends with ".egsphant"/".begsphant". Then
	those, appropriate options are enabled instead or the egsphant
	(likely generated in DICOM_to_egsphant) file is loaded into
	phant.  If a file fails to be read in as DICOM, the code
	terminates.
	
	filterLowDensity=X sets the activity to any voxel with density
	less than X equal to zero.
	
	filterLowActivity=Y sets any activity below Y*maxActivity equal
	to zero.
	
	outputImages outputs PNGs of each slice of the output phantoms
	for both media and density, as well as a blue->red wash of 
	activity assigned to phant in Activity.txt.  It is very time
	intensive, so best to be used to check activity and
	registration.
	*/
	bool outputImages = false;
	double filterLowDensity = 0;
	double filterLowActivity = 0.01;
	
	if (argc == 1) {
        std::cout << "Please call this program with one or more .dcm files and a .egsphant or .begsphant file.\n";
        return 0;
    }

    database dat;
    QVector <DICOM *> dicom;
    QVector <DICOM *> dicomExtra;
	EGSPhant phant;

    for (int i = 0; i < argc-1; i++) {
        QString path(argv[i+1]);
        if (!path.compare("-outputImages"))
			outputImages = true;
        else if (!path.left(17).compare("filterLowDensity="))
			filterLowDensity = path.right(path.size()-17).toDouble();
        else if (!path.left(18).compare("filterLowActivity="))
			filterLowActivity = path.right(path.size()-18).toDouble();
        else if (path.endsWith(".egsphant"))
			phant.loadEGSPhantFilePlus(path);
		else if (path.endsWith(".begsphant"))
			phant.loadbEGSPhantFilePlus(path);			
		else {			
			DICOM *d = new DICOM(&dat);
			if (d->parse(path)) {
				dicom.append(d);
			}
			else {
				std::cout << "Unsuccessfully parsed " << path.toStdString() << ", quitting...\n";
				for (int j = 0; j < dicom.size(); j++) {
					delete dicom[j];
				}
				dicom.clear();
				delete d;
				return -1;
			}
		}
    }
	if (!phant.nx && !phant.ny && !phant.nz) {
		std::cout << "Did not receive .egsphant or .begsphant input, quitting...\n";
		return -1;
	}
	
	duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
	std::cout << "Parsed the " << dicom.size() << " DICOM files.  Time elapsed is " << duration << " s.\n";
	
	// ---------------------------------------------------------- //
	// INITIAL DICOM INFORMATION FETCHING                         //
	// ---------------------------------------------------------- //
	/*
	In this section, we define all the variables that are used to
	store the DICOM data.  Then, we go through every DICOM file to
	check their tag[0008,0008] field to ensure that it contains 
	ORIGINAL, as a quick and easy way to see if we have an image
	file or not.  Then we seperate the DICOM files that don't have
	the field into dicomExtra, and now the QVector dicom, holds all
	the slice info.
	
	As it is not safe to assume that the data is given in proper
	ascending order, all the dicom data is sorted by the z value
	extracted as the code was parsed.  The z value is given by the
	tag[0020,1041] which provides slice height.  A classic, non-
	recursive merge sort is used to sort them.
	*/

	// Sort out all the DICOM data into the following
	double rescaleM = 1, rescaleB = 0, rescaleFlag;
	int bitsStored = 16, bytesStored = 2;
    QVector <QVector <QVector <unsigned short int> > > HU;
    QVector <unsigned short int> xPix;
    QVector <unsigned short int> yPix;
    QVector <QVector <double> > imagePos;
    QVector <QVector <double> > xySpacing;
    QVector <double> zSpacing;
	
	bool ctFlag = false;
	
	// Put not-CT dicom in seperate array
    for (int i = 0; i < dicom.size(); i++) {
		ctFlag = false;
        for (int j = 0; j < dicom[i]->data.size(); j++)
            if (dicom[i]->data[j]->tag[0] == 0x0008 && dicom[i]->data[j]->tag[1] == 0x0008) {
                QString temp = "";
                for (unsigned int s = 0; s < dicom[i]->data[j]->vl; s++) {
                    temp.append(dicom[i]->data[j]->vf[s]);
                }
				
                if (temp.contains("ORIGINAL")) {
					ctFlag = true;
				}
			}
		if (!ctFlag) {
			dicomExtra.append(dicom[i]);
			dicom.remove(i--);
		}
	}
		
	// Sort all CT slices by z height
	mergeSort(dicom,dicom.size());
	
	duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
	std::cout << "Sorted the " << dicom.size() << " DICOM CT files along z.  Time elapsed is " << duration << " s.\n";


	// ---------------------------------------------------------- //
	// FETCH HU DATA FROM DICOM FILES                             //
	// ---------------------------------------------------------- //
	/*
	In this section, we extract all the data we require from the
	dicom data.  Here are the tags and the relevant data we fetch
	from each of them:
	
	Tag			Representation
	0028,0030	Pixel thickness for x and y
	0018,0050	Pixel thickness for z of the current slice
	0028,0101	Bits stored for each HU
	0020,0032	Position of the center of the most -x and +y voxel
	0028,0010	Number of voxels in x
	0028,0011	Number of voxels in y
	0028,1053	The m in the equation "HU_f = m*HU_i+b"
	0028,1052	The b in the equation "HU_f = m*HU_i+b"
	7fe0,0010	The unscaled HU data (HU_i) above for all voxels
	
	Careful in that most DICOM data is in mm, and we are converting
	to cm for the egsphant.  The HU data starts at the most -x and
	+y voxel and iterates first by x, and then by y.
	
	The number of voxels in z is then therefore simply the number of
	files we have iterated through.  It is important to note that
	the z index	is now the most top-level iterator for all the
	above arrays.
	*/
	
    for (int i = 0; i < dicom.size(); i++) {
		rescaleFlag = 0;
        for (int j = 0; j < dicom[i]->data.size(); j++) {
            // Pixel Spacing (Decimal String), row spacing and then column spacing (in mm)
            if (dicom[i]->data[j]->tag[0] == 0x0028 && dicom[i]->data[j]->tag[1] == 0x0030) {
                xySpacing.resize(xySpacing.size()+1);
                xySpacing.last().resize(2);

                QString temp = "";
                for (unsigned int s = 0; s < dicom[i]->data[j]->vl; s++) {
                    temp.append(dicom[i]->data[j]->vf[s]);
                }

                xySpacing.last()[0] = (temp.split('\\',QString::SkipEmptyParts)[0]).toDouble();
                xySpacing.last()[1] = (temp.split('\\',QString::SkipEmptyParts)[1]).toDouble();
            } // Slice Thickness (Decimal String, in mm)
            else if (dicom[i]->data[j]->tag[0] == 0x0018 && dicom[i]->data[j]->tag[1] == 0x0050) {
                QString temp = "";
                for (unsigned int s = 0; s < dicom[i]->data[j]->vl; s++) {
                    temp.append(dicom[i]->data[j]->vf[s]);
                }

                zSpacing.append(temp.toDouble());
            } // Bits Stored
            else if (dicom[i]->data[j]->tag[0] == 0x0028 && dicom[i]->data[j]->tag[1] == 0x0101) {
				if (dicom[i]->isBigEndian)
                    bitsStored = ((unsigned short int)(((short int)(dicom[i]->data[j]->vf[0]) << 8) +
					(short int)(dicom[i]->data[j]->vf[1])));
                else
                    bitsStored = ((unsigned short int)(((short int)(dicom[i]->data[j]->vf[1]) << 8) +
					(short int)(dicom[i]->data[j]->vf[0])));
				
				bytesStored = bitsStored/8;
            } // Image Position [x,y,z] (Decimal String, in mm)
            else if (dicom[i]->data[j]->tag[0] == 0x0020 && dicom[i]->data[j]->tag[1] == 0x0032) {
                imagePos.resize(imagePos.size()+1);
                imagePos.last().resize(3);

                QString temp = "";
                for (unsigned int s = 0; s < dicom[i]->data[j]->vl; s++) {
                    temp.append(dicom[i]->data[j]->vf[s]);
                }
				
                imagePos.last()[0] = (temp.split('\\',QString::SkipEmptyParts)[0]).toDouble();
                imagePos.last()[1] = (temp.split('\\',QString::SkipEmptyParts)[1]).toDouble();
                imagePos.last()[2] = (temp.split('\\',QString::SkipEmptyParts)[2]).toDouble();
            } // Rows
            else if (dicom[i]->data[j]->tag[0] == 0x0028 && dicom[i]->data[j]->tag[1] == 0x0010) {
				if (dicom[i]->isBigEndian)
                    xPix.append((unsigned short int)(((short int)(dicom[i]->data[j]->vf[0]) << 8) +
					(short int)(dicom[i]->data[j]->vf[1])));
                else
                    xPix.append((unsigned short int)(((short int)(dicom[i]->data[j]->vf[1]) << 8) +
					(short int)(dicom[i]->data[j]->vf[0])));	
            } // Columns
            else if (dicom[i]->data[j]->tag[0] == 0x0028 && dicom[i]->data[j]->tag[1] == 0x0011) {
				if (dicom[i]->isBigEndian)
                    yPix.append((unsigned short int)(((short int)(dicom[i]->data[j]->vf[0]) << 8) +
					(short int)(dicom[i]->data[j]->vf[1])));
                else
                    yPix.append((unsigned short int)(((short int)(dicom[i]->data[j]->vf[1]) << 8) +
					(short int)(dicom[i]->data[j]->vf[0])));		
            } // Rescale HU slope (assuming type is HU)
            else if (dicom[i]->data[j]->tag[0] == 0x0028 && dicom[i]->data[j]->tag[1] == 0x1053) {
                QString temp = "";
                for (unsigned int s = 0; s < dicom[i]->data[j]->vl; s++) {
                    temp.append(dicom[i]->data[j]->vf[s]);
                }
				
				rescaleM = temp.toDouble();
				rescaleFlag++;
            } // Rescale HU intercept (assuming type is HU)
            else if (dicom[i]->data[j]->tag[0] == 0x0028 && dicom[i]->data[j]->tag[1] == 0x1052) {
                QString temp = "";
                for (unsigned int s = 0; s < dicom[i]->data[j]->vl; s++) {
                    temp.append(dicom[i]->data[j]->vf[s]);
                }
				
				rescaleB = temp.toDouble();
				rescaleFlag++;
            } // HU values (assuming 2-bytes as I have yet to encounter anything different, ie, assumes TAG (0028,0100) = 16)
            else if (dicom[i]->data[j]->tag[0] == 0x7fe0 && dicom[i]->data[j]->tag[1] == 0x0010) {
                HU.resize(HU.size()+1);
				if (HU.size() == xPix.size() && HU.size() == yPix.size()) {
                    HU.last().resize(yPix.last());
                    for (unsigned int k = 0; k < yPix.last(); k++) {
                        HU.last()[k].resize(xPix.last());
                    }
					
					unsigned short int temp;
                    if (dicom[i]->isBigEndian)
                        for (unsigned int s = 0; s < dicom[i]->data[j]->vl; s+=bytesStored) {
							temp = 0;
							for (int ss = 0; ss < bytesStored; ss++)
								temp += (dicom[i]->data[j]->vf[s+ss]) << (bitsStored-((ss+1)*8));
							temp = rescaleFlag == 2 ? rescaleM*temp+rescaleB : temp;
							
							HU.last()[int(int(s/bytesStored)/xPix.last())][int(s/bytesStored)%xPix.last()] = temp;
						}
                    else
                        for (unsigned int s = 0; s < dicom[i]->data[j]->vl; s+=bytesStored) {
							temp = 0;
							for (int ss = 0; ss < bytesStored; ss++)
								temp += (dicom[i]->data[j]->vf[s+ss]) << (ss*8);
							temp = rescaleFlag == 2 ? rescaleM*temp+rescaleB : temp;
														
							HU.last()[int(int(s/bytesStored)/xPix.last())][int(s/bytesStored)%xPix.last()] = temp;
						}
                }
            }
        }
    }
	
	if (HU.size() > 0) {
		duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
		std::cout << "Extracted all HU data for the " << xPix[0] << "x" << yPix[0] << " slices.  Time elapsed is " << duration << " s.\n";
	}
	else
		std::cout << "Did not find HU data.\n";
	
	// ---------------------------------------------------------- //
	// CONVERTING HU TO ACTIVITY LIST                             //
	// ---------------------------------------------------------- //
	/*
	The boundaries of the egsphant activity are set based on DICOM
	data.  The HU data for each voxel [i,j,k] is then used to
	assign activity to the [i,ny-1-j,k] voxel of the egsphant.
	This inverts the y-axis, so that it counts up from -y rather
	than from +y.
	*/
	
	EGSPhant activity;
	activity.nx = xPix[0];
	activity.ny = yPix[0];
	activity.nz = dicom.size();
    activity.x.fill(0,activity.nx+1);
    activity.y.fill(0,activity.ny+1);
    activity.z.fill(0,activity.nz+1);
	{
		QVector <char> mz(activity.nz, 1);
		QVector <QVector <char> > my(activity.ny, mz);
		QVector <QVector <QVector <char> > > mx(activity.nx, my);
		activity.m = mx;
		QVector <double> dz(activity.nz, 0);
		QVector <QVector <double> > dy(activity.ny, dz);
		QVector <QVector <QVector <double> > > dx(activity.nx, dy);
		activity.d = dx;
	}
	activity.media.append("DUMMY");
	
	// Define xy bound values, still assuming first slice matches the rest
    for (int i = 0; i <= activity.nx; i++)
		activity.x[i] = (imagePos[0][0]+(i-0.5)*xySpacing[0][0])/10.0;
    for (int i = 0; i <= activity.ny; i++)
		activity.y[i] = (imagePos[0][1]+(i-0.5)*xySpacing[0][1])/10.0;
	
	// Define z bound values
	double prevZ, nextZ;
	nextZ = imagePos[0][2]-zSpacing[0]/2.0;
	for (int i = 0; i < activity.nz; i++) {
		prevZ = nextZ/2.0 + (imagePos[i][2]-zSpacing[i]/2.0)/2.0;
		nextZ = imagePos[i][2]+zSpacing[i]/2.0;
		activity.z[i] = prevZ/10.0;
	}
	activity.z.last() = nextZ/10.0;
	
	// Write HU into density array (which is actually activity)
	double maxAct = 0;
	int nj = 0;
	for (int k = 0; k < activity.nz; k++) // Z //
		for (int j = 0; j < activity.ny; j++) { // Y //
			nj = activity.ny-1-j; // Reversed j index for density and media assignment
			for (int i = 0; i < activity.nx; i++) { // X //
				activity.d[i][nj][k] = HU[k][j][i];
				maxAct = (maxAct<HU[k][j][i])?HU[k][j][i]:maxAct;
			}
		}
			
			
	duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
    std::cout << "Succesfully generated activity matrix (dimensions x: ["
			  << activity.x[0] << "," << activity.x[activity.nx] << "], y:["
			  << activity.y[0] << "," << activity.y[activity.ny] << "], z:["
			  << activity.z[0] << "," << activity.z[activity.nz] << "]).  Time elapsed is " << duration << " s.\n";
	
	// ---------------------------------------------------------- //
	// OUTPUT ACTIVITY FILE                                       //
	// ---------------------------------------------------------- //
	/*
	In this section, we iterate through the imported egsphant
	phant and check the activity from the DICOM data (which is
	stored as density in the egsphant activity).
	
	We iterate through z, y, and x of phant.  At the center of
	each voxel, we invoke interpDen (interpolate density) in
	egsphant activity to find the activity at the center of the
	phant voxel.  Then, we output the region of the voxel and the
	activity we just computed to Activity.txt.
	*/
		
	// Output activity to file
	double minAct = maxAct*filterLowActivity, tempAct, zMid, yMid, xMid; // Threshold cutoff defined here
	QFile outputF("Activity.txt");
	
    if (outputF.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream output(&outputF);
		for (int k = 0; k < phant.nz; k++) { // Z //
			zMid = (phant.z[k]+phant.z[k+1])/2.0;
			for (int j = 0; j < phant.ny; j++) { // Y //
				yMid = (phant.y[j]+phant.y[j+1])/2.0;
				for (int i = 0; i < phant.nx; i++) { // X //
					if (phant.d[i][j][k] >= filterLowDensity) {
						xMid = (phant.x[i]+phant.x[i+1])/2.0;
						tempAct = activity.interpDen(xMid, yMid, zMid);
						
						if (tempAct > minAct)
							output << i+j*phant.nx+k*phant.nx*phant.ny << " " << tempAct << "\n";
					}
				}
			}
		}
	}
	else {
		std::cout << "Failed to open Activity.txt, quitting...\n";
		return 0;
	}
	outputF.close();
	
	duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
    std::cout << "Activity file successfully output.  Time elapsed is " << duration << " s.\n";

	// ---------------------------------------------------------- //
	// OUTPUT IMAGES                                              //
	// ---------------------------------------------------------- //
	/*
	In this section, we output phant images using the prebuilt
	image functions in the egsphant class.  There is also an
	additional step in drawing a blue to red color wash of
	activity as assigned in the Activity.txt file.
	*/
	if (outputImages) {
		double xi = (phant.x[0]+phant.x[1])/2.0;
		double xf = (phant.x[phant.nx-1]+phant.x[phant.nx])/2.0;
		double yi = (phant.y[0]+phant.y[1])/2.0;
		double yf = (phant.y[phant.ny-1]+phant.y[phant.ny])/2.0;
		double res = 1/(phant.x[1]-phant.x[0])*2; // This sets resolution to be 2 pixels for each voxel in x
		QImage tempM, tempD;
		QPen pen;
		pen.setWidth(1);
		for (int k = 0; k < phant.nz; k++) {
			zMid = (phant.z[k]+phant.z[k+1])/2.0;
			tempM = phant.getEGSPhantPicMed("z axis", yi, yf, xi, xf, zMid, res);
			tempD = phant.getEGSPhantPicDen("z axis", yi, yf, xi, xf, zMid, res);
			
			QPainter paintM (&tempM), paintD (&tempD);
			for (int j = 0; j < phant.ny; j++) { // Y //
				yMid = (phant.y[j]+phant.y[j+1])/2.0;
				for (int i = 0; i < phant.nx; i++) { // X //
					xMid = (phant.x[i]+phant.x[i+1])/2.0;
					
					tempAct = activity.interpDen(xMid, yMid, zMid);
					if (tempAct > minAct && phant.d[i][j][k] >= filterLowDensity) {
						tempAct = (tempAct-minAct)/(maxAct-minAct);
						tempAct = tempAct > 1.0 ? 1.0 : tempAct;
						pen.setColor(QColor(tempAct*255.0,0,(1.0-tempAct)*255.0));
						
						paintM.setPen(pen); paintM.drawPoint(i*2, (phant.ny-1-j)*2); paintM.drawPoint(i*2+1, (phant.ny-1-j)*2+1);
						paintD.setPen(pen); paintD.drawPoint(i*2, (phant.ny-1-j)*2); paintD.drawPoint(i*2+1, (phant.ny-1-j)*2+1);
					}
				}
			}
			tempM.save(QString("Image/MedPic")+QString::number(k+1)+".png");
			tempD.save(QString("Image/DenPic")+QString::number(k+1)+".png");
		}
		
		duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
		std::cout << "Image data successfully output.  Time elapsed is " << duration << " s.\n";
	}
	
	// ---------------------------------------------------------- //
	// CLEAR MEMORY ALLOCATIONS                                   //
	// ---------------------------------------------------------- //
	/*
	This section just deallocates the assigned memory.  Hope I did
	not miss anything.
	*/
	
    for (int j = 0; j < dicom.size(); j++) {
        delete dicom[j];
    }
    for (int j = 0; j < dicomExtra.size(); j++) {
        delete dicomExtra[j];
    }
    dicom.clear();
	dicomExtra.clear();
    return 1;
}