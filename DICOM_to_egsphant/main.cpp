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
	"-makeMasks", "-outputImages", "-nominalDensity", or
	"tag=string", then those appropriate options are enabled on
	instead.  If a file fails to be read in as DICOM, the code terminates.
	
	tag=string changes the lookup of the priority and TAS files 
	from Default to string, ie,
		Default_priority.txt -> string_priority.txt 
		Default_all_TAS.txt -> string_all_TAS.txt
		
	makeMasks outputs begsphant files for each structure with
	uniform density and two possible media, NONTARGET and TARGET,
	used to identify structures later for post-processing apps
	such as 3ddose_tools.
	
	outputImages outputs PNGs of each slice of the output phantoms
	for both media and density, as well as outlines of structures
	over the media images.  It is very time intensive, so best to
	be used to check TAS and registration.
	
	nominalDensity uses the file Default_mediaDensity.txt (where
	the file lookup would change for the appropriate tag) to
	assign densities to all media.
	*/
	
	bool makeMasks = false;
	bool outputImages = false;
	bool nominalDensity = false;
	QString TAS_tag("Default");
	
	if (argc == 1) {
        std::cout << "Please call this program with one or more .dcm files.\n";
        return 0;
    }

    database dat;
    QVector <DICOM *> dicom;
    QVector <DICOM *> dicomExtra;

    for (int i = 0; i < argc-1; i++) {
        QString path(argv[i+1]);
        DICOM *d = new DICOM(&dat);
        if (!path.compare("-outputImages"))
			outputImages = true;
		else if (!path.compare("-makeMasks"))
			makeMasks = true;
		else if (!path.compare("-nominalDensity"))
			nominalDensity = true;
		else if (!path.left(4).compare("tag="))
			TAS_tag = path.right(path.size()-4);
		else if (d->parse(path)) {
            //debug//std::cout << std::dec << "Successfully parsed " << path.toStdString() << ".\n";
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
	
	// Options
	
	duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
	std::cout << "Parsed the " << dicom.size() << " DICOM files.  Time elapsed is " << duration << " s.\n";
	
	// ---------------------------------------------------------- //
	// INITIAL DICOM INFORMATION FETCHING                         //
	// ---------------------------------------------------------- //
	/*
	In this section, we define all the variables that are used to
	store the DICOM data.  Then, we go through every DICOM file to
	check their tag[0008,0008] field to ensure that it contains 
	AXIAL, as a quick and easy way to see if we have an image file
	or not.  Then we seperate the DICOM files that don't have the
	field into dicomExtra, and now the QVector dicom, holds all the
	ct slice info.
	
	As it is not safe to assume that the CT data is given in proper
	ascending order, all the dicom data is sorted by the z value
	extracted as the code was parsed.  The z value is given by the
	tag[0020,1041] which provides slice height.  A classic, non-
	recursive merge sort is used to sort them.
	*/

	// Sort out all the DICOM data into the following
	double rescaleM = 1, rescaleB = 0, rescaleFlag;
    QVector <QVector <QVector <short int> > > HU;
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
				
                if (temp.contains("AXIAL")) {
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
	// FETCH CT DATA FROM DICOM FILES                             //
	// ---------------------------------------------------------- //
	/*
	In this section, we extract all the data we require from the
	dicom data.  Here are the tags and the relevant data we fetch
	from each of them:
	
	Tag			Representation
	0028,0030	Pixel thickness for x and y
	0018,0050	Pixel thickness for z of the current slice
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
            } // HU values
            else if (dicom[i]->data[j]->tag[0] == 0x7fe0 && dicom[i]->data[j]->tag[1] == 0x0010) {
                HU.resize(HU.size()+1);
				if (HU.size() == xPix.size() && HU.size() == yPix.size()) {
                    HU.last().resize(yPix.last());
                    for (unsigned int k = 0; k < yPix.last(); k++) {
                        HU.last()[k].resize(xPix.last());
                    }
					
					short int temp;
                    if (dicom[i]->isBigEndian)
                        for (unsigned int s = 0; s < dicom[i]->data[j]->vl; s+=2) {
                            temp  = (dicom[i]->data[j]->vf[s+1]);
							temp += (short int)(dicom[i]->data[j]->vf[s]) << 8;
							
							HU.last()[int(int(s/2)/xPix.last())][int(s/2)%xPix.last()] =
								rescaleFlag == 2 ? rescaleM*temp+rescaleB : temp;
						}
                    else
                        for (unsigned int s = 0; s < dicom[i]->data[j]->vl; s+=2) {
                            temp  = (dicom[i]->data[j]->vf[s]);
							temp += (short int)(dicom[i]->data[j]->vf[s+1]) << 8;
							
							HU.last()[int(int(s/2)/xPix.last())][int(s/2)%xPix.last()] =
								rescaleFlag == 2 ? rescaleM*temp+rescaleB : temp;
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
	// FETCH RS FILE STRUCTURE DATA                               //
	// ---------------------------------------------------------- //
	/*
	In this section, we look for structure data in the dicomExtra
	array.  Like above, here are the tags and the relevant data we
	fetch from each of them:
	
	Tag			Representation
	3006,0020	Parent to the following:
	*3006,0026	*Structure name
	*3006,0022	*Structure number
	3006,0039	Parent to the following:
	*3006,0040	*Parent to the following:
	**3006,0050	**Full contour data
	*3006,0084	*Structure number
	
	* nested sequence, ** double-nested sequence
	
	As the structure data that we are interested in (name and
	contours) are in seperate tag groups, we need to extract the
	structure number twice, to match structure name to the
	contour data.
	
	This section involves reading SQ at one and two layers deep,
	so the parseSequence function is invoked several times to
	parse the appropriate data.
	
	The structures (which are many sets of [x,y,z] positions)
	into an array of polygons and an array of z positions, for
	each structure.
	*/
	
	// Create struct arrays ---------------------------------- //
	QVector <QVector <QPolygonF> > structPos;
	QVector <QVector <double> > structZ;
	QVector <QString> structName;
	QVector <int> structNum;
	QMap <int, int> structLookup;
	QVector <int> structReference;
	
	// Data for parsing singly and doubly nested SQ sets and point data strings
	QVector <Attribute*> *att, *att2;
	QByteArray tempData, tempData2;
	QStringList pointData;
	
    for (int i = 0; i < dicomExtra.size(); i++) {
        for (int j = 0; j < dicomExtra[i]->data.size(); j++) {
            // Structure info (looking for structure names and nums)
            if (dicomExtra[i]->data[j]->tag[0] == 0x3006 && dicomExtra[i]->data[j]->tag[1] == 0x0020) {
				for (int k = 0; k < dicomExtra[i]->data[j]->seq.items.size(); k++) {
					tempData = QByteArray((char*)dicomExtra[i]->data[j]->seq.items[k]->vf,dicomExtra[i]->data[j]->seq.items[k]->vl);
					QDataStream dataStream(tempData);
					
					att = new QVector <Attribute*>;
					if (!dicomExtra[i]->parseSequence(&dataStream, att)) {
						std::cout << "Failed to parse sequence data for tag (3006,0020), quitting...\n";
						return 0;
					}
					
					QString tempS = ""; // Get the name
					QString tempI = ""; // Get the number
					for (int l = 0; l < att->size(); l++) {
						if (att->at(l)->tag[0] == 0x3006 && att->at(l)->tag[1] == 0x0026)
							for (unsigned int s = 0; s < att->at(l)->vl; s++)
								tempS.append(att->at(l)->vf[s]);
						else if (att->at(l)->tag[0] == 0x3006 && att->at(l)->tag[1] == 0x0022)
							for (unsigned int s = 0; s < att->at(l)->vl; s++)
								tempI.append(att->at(l)->vf[s]);
					}
					tempS = tempS.trimmed().replace(' ','_');
					
					structName.append(tempS.trimmed());
					structNum.append(tempI.toInt());
					structLookup[tempI.toInt()] = structName.size()-1;
					
					for (int l = 0; l < att->size(); l++)
						delete att->at(l);
					delete att;
				}
			} // Structure data (looking for contour definitions)
			else if (dicomExtra[i]->data[j]->tag[0] == 0x3006 && dicomExtra[i]->data[j]->tag[1] == 0x0039) {
				for (int k = 0; k < dicomExtra[i]->data[j]->seq.items.size(); k++) {
					tempData = QByteArray((char*)dicomExtra[i]->data[j]->seq.items[k]->vf,dicomExtra[i]->data[j]->seq.items[k]->vl);
					QDataStream dataStream(tempData);
					
					att = new QVector <Attribute*>;
					if (!dicomExtra[i]->parseSequence(&dataStream, att)) {
						std::cout << "Failed to parse sequence data for tag (3006,0039), quitting...\n";
						return 0;
					}
					
					QString tempS = ""; // Get the contour, it's another nested sequence, so we must go deeper with parseSequence
					for (int l = 0; l < att->size(); l++) {
						if (att->at(l)->tag[0] == 0x3006 && att->at(l)->tag[1] == 0x0040)
						{
							structZ.resize(structZ.size()+1);
							structPos.resize(structPos.size()+1);
							for (int k = 0; k < att->at(l)->seq.items.size(); k++) {
								structPos.last().resize(structPos.last().size()+1);
								tempData2 = QByteArray((char*)att->at(l)->seq.items[k]->vf,att->at(l)->seq.items[k]->vl);
								QDataStream dataStream2(tempData2);
								
								att2 = new QVector <Attribute*>;
								if (!dicomExtra[i]->parseSequence(&dataStream2, att2)) {
									std::cout << "Failed to parse sequence data for tag (3006,0040), quitting...\n";
									return 0;
								}
								
								QString tempS = ""; // Get the points
								for (int m = 0; m < att2->size(); m++)
									if (att2->at(m)->tag[0] == 0x3006 && att2->at(m)->tag[1] == 0x0050)
										for (unsigned int s = 0; s < att2->at(m)->vl; s++)
											tempS.append(att2->at(m)->vf[s]);
								
								pointData = tempS.split('\\');
								structZ.last().append(pointData[2].toDouble()/10.0);
								for (int m = 0; m < pointData.size(); m+=3)
									structPos.last().last() << QPointF(pointData[m].toDouble()/10.0, pointData[m+1].toDouble()/10.0);
								
								for (int m = 0; m < att2->size(); m++)
									delete att2->at(m);
								delete att2;
							}
						}
						else if (att->at(l)->tag[0] == 0x3006 && att->at(l)->tag[1] == 0x0084) {
							QString tempI = ""; // Get the number
							for (unsigned int s = 0; s < att->at(l)->vl; s++)
								tempI.append(att->at(l)->vf[s]);
							structReference.append(tempI.toInt());
						}
					}
										
					for (int l = 0; l < att->size(); l++)
						delete att->at(l);
					delete att;
				}
			}
		}
	}
	duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
	std::cout << "Extracted data for all " << structName.size() << " structures.  Time elapsed is " << duration << " s.\n";

	// ---------------------------------------------------------- //
	// GENERATE EGSPHANT AND HU CONVERSION MAPS/THRESHOLDS        //
	// ---------------------------------------------------------- //
	/*
	In this section, we fetch the density curve, the priorities,
	and all of the TAS files associated with the struct names.
	Though there is a lot of code, it is all very simple and
	similar.
	
	The latter half of this section defines an empty egsphant of 
	the correct dimensions based on the previously retrieved CT
	data.  It also adds all the different media in the given TAS
	files to the egsphant.
	
	The z dimensions are taken as the average of the z bounds of
	two neighbouring slices (as z values don't always line up).
	*/

	// Generate an egsphant with the given data -------------- //
	EGSPhant phant;
	
	// Assignment schemes
	QString density("Assignment_Schemes/"+TAS_tag+"_densityCurve.txt");
	QString medium("Assignment_Schemes/"+TAS_tag+"_all_TAS.txt");
	QString priority("Assignment_Schemes/"+TAS_tag+"_priority.txt");
	QString nominal("Assignment_Schemes/"+TAS_tag+"_mediaDensity.txt");
	
	QVector <double> HUMap, denMap; // HU to density lookups
	
	QVector <QString> medThreshold; // Default media names for TAS threshold lookup
	QVector <QVector <QString> > medThresholds; // Alternate media names for TAS threshold lookup
	QVector <double> denThreshold; // Default densities for TAS threshold lookup
	QVector <QVector <double> > denThresholds; // Alternate densities for TAS threshold lookup
	QVector <int> structPrio; // This array holds the overlap priority of different structures
	QVector <int> structUnique; // This array dictates whether or not to do alternate TAS
	
	QVector <QString> medNom; // List of media to use with denNom to assign density
	QVector <double> denNom; // Nominal densities to be assigned to voxels
	
	// Read in the HU to density conversion data
	QFile* file = new QFile (density);
	QString tempS = "";
	if (file->open(QIODevice::ReadOnly | QIODevice::Text)) {
		QTextStream input(file);
		
		while (!input.atEnd()) {
			tempS = input.readLine();
			
			if (tempS.contains(" ")) {
				denMap << tempS.split(' ', QString::SkipEmptyParts)[0].toDouble();
				HUMap << tempS.split(' ', QString::SkipEmptyParts)[1].toDouble();
			}
			else if (tempS.contains("\t")) {
				denMap << tempS.split('\t', QString::SkipEmptyParts)[0].toDouble();
				HUMap << tempS.split('\t', QString::SkipEmptyParts)[1].toDouble();				
			}
			else {
				std::cout << "Did not find a space or tab delimiter in file " << density.toStdString() << ", quitting...\n";
				return -1;
			}
		}
    }
	else {
		std::cout << "Did not find general media TAS file " << density.toStdString() << ", quitting...\n";	
		return -1;	
	}
	delete file;
	
	// Read in the default TAS
	file = new QFile (medium);
	if (file->open(QIODevice::ReadOnly | QIODevice::Text)) {
		QTextStream input(file);
		
		while (!input.atEnd()) {
			tempS = input.readLine();
			
			if (tempS.contains(" ")) {
				medThreshold << tempS.split(' ', QString::SkipEmptyParts)[0];
				denThreshold << tempS.split(' ', QString::SkipEmptyParts)[1].toDouble();
			}
			else if (tempS.contains("\t")) {
				medThreshold << tempS.split('\t', QString::SkipEmptyParts)[0];
				denThreshold << tempS.split('\t', QString::SkipEmptyParts)[1].toDouble();				
			}
			else {
				std::cout << "Did not find a space or tab delimiter in file " << medium.toStdString() << ", quitting...\n";
				return -1;
			}
		}
    }
	else {
		std::cout << "Did not find general media TAS file " << medium.toStdString() << ", quitting...\n";	
		return -1;	
	}
	delete file;
	
	// Assign default priorities
	for (int i = 0; i < structName.size(); i++)
		structPrio.append(structName.size()-i+1);
	
	// Read in any alternate priorities
	file = new QFile(priority);
	if (file->open(QIODevice::ReadOnly | QIODevice::Text) && structName.size() > 0) {
		QTextStream input(file);
		
		while (!input.atEnd()) {
			tempS = input.readLine();
			QString tempS2 = "";
			int i = 0;
			
			if (tempS.contains(" ")) {
				tempS2 = tempS.split(' ', QString::SkipEmptyParts)[0];
				for (i = 0; i < structName.size(); i++)
					if (!structName[i].compare(tempS2))
						break;
					
				if (i < structName.size() - 1 || !structName.last().compare(tempS2))
					structPrio[i] = tempS.split(' ', QString::SkipEmptyParts)[1].toInt();
			}
			else if (tempS.contains("\t")) {
				tempS2 = tempS.split('\t', QString::SkipEmptyParts)[0];
				for (i = 0; i < structName.size(); i++)
					if (!structName[i].compare(tempS2))
						break;
					
				if (i < structName.size() - 1 || !structName.last().compare(tempS2))
					structPrio[i] = tempS.split('\t', QString::SkipEmptyParts)[1].toInt();				
			}
			else {
				std::cout << "Did not find a space or tab delimiter in file " << priority.toStdString() << ", quitting...\n";
				return -1;
			}
		}
    }
	else {
		std::cout << "Did not find general priority file " << priority.toStdString() << ", assuming ordered weighting...\n";
	}
	delete file;
	
	std::cout << "Structure priorities (higher score overrules):\n";
	for (int i = 0; i < structName.size(); i++)
		std::cout << "\t" << structName[i].toStdString() << " - " << structPrio[i] << "\n";
	
	// Read in alternate TAS
	for (int i = 0; i < structName.size(); i++) {
		file = new QFile("Assignment_Schemes/"+TAS_tag+"_"+structName[i]+"_TAS.txt");
		medThresholds.resize(medThresholds.size()+1);
		denThresholds.resize(denThresholds.size()+1);
		
		if (file->open(QIODevice::ReadOnly | QIODevice::Text) && structName.size() > 0) {
			QTextStream input(file);
			while (!input.atEnd()) {
				tempS = input.readLine();
				
				if (tempS.contains(" ")) {
					medThresholds.last() << tempS.split(' ', QString::SkipEmptyParts)[0];
					denThresholds.last() << tempS.split(' ', QString::SkipEmptyParts)[1].toDouble();
				}
				else if (tempS.contains("\t")) {
					medThresholds.last() << tempS.split('\t', QString::SkipEmptyParts)[0];
					denThresholds.last() << tempS.split('\t', QString::SkipEmptyParts)[1].toDouble();				
				}
				else {
					std::cout << "Did not find a space or tab delimiter in file " << ("Assignment_Schemes/"+TAS_tag+"_"+structName[i]+"_TAS.txt").toStdString() << ", quitting...\n";
					return -1;
				}
			}
			//std::cout << "Found " << ("Assignment_Schemes/"+TAS_tag+"_"+structName[i]+"_TAS.txt").toStdString() << ".\n";
			structUnique << true;
		}
		else {
			std::cout << "Did not find " << ("Assignment_Schemes/"+TAS_tag+"_"+structName[i]+"_TAS.txt").toStdString() << ", assuming default media assignment.\n";
			medThresholds.last() = medThreshold;
			denThresholds.last() = denThreshold;
			structUnique << false;
		}
		delete file;		
	}
	
	// Read in nominal densities (if requested)
	if (nominalDensity) {
		file = new QFile(nominal);
		if (file->open(QIODevice::ReadOnly | QIODevice::Text)) {
			QTextStream input(file);
			
			while (!input.atEnd()) {
				tempS = input.readLine();
				
				if (tempS.contains(" ")) {
					medNom.append(tempS.split(' ', QString::SkipEmptyParts)[0]);
					denNom.append(tempS.split(' ', QString::SkipEmptyParts)[1].toDouble());
				}
				else if (tempS.contains("\t")) {
					medNom.append(tempS.split('\t', QString::SkipEmptyParts)[0]);
					denNom.append(tempS.split('\t', QString::SkipEmptyParts)[1].toDouble());			
				}
				else {
					std::cout << "Did not find a space or tab delimiter in file " << nominal.toStdString() << ", quitting...\n";
					return -1;
				}
			}
		}
		else {
			std::cout << "Did not find nominal density file " << nominal.toStdString() << ", assigning density based on HU...\n";
		}
		delete file;
	}
	
	// Assume first slice matches the rest and set x, y, and z boundaries
	phant.nx = xPix[0];
	phant.ny = yPix[0];
	phant.nz = dicom.size();
    phant.x.fill(0,phant.nx+1);
    phant.y.fill(0,phant.ny+1);
    phant.z.fill(0,phant.nz+1);
	
	// Track the number of media, and create a lookup for media to ASCII character representation
	QMap <QString,unsigned char> mediaMap;
	int medNum = 0;
	for (int i = 0; i < medThreshold.size(); i++) {
		phant.media.append(medThreshold[i]);
		mediaMap.insert(medThreshold[i], 49 + medNum + (medNum>8?7:0) + (medNum>34?6:0));
		medNum++; // 49 is the ASCII value of '1', 7 jumps from ';' to 'A', 6 jumps from '[' to 'a'
	}
	
	for (int i = 0; i < structName.size(); i++) {
		for (int j = 0; j < medThresholds[i].size(); j++)
			if (!mediaMap.contains(medThresholds[i][j])) {
				phant.media.append(medThresholds[i][j]);
				mediaMap.insert(medThresholds[i][j], 49 + medNum + (medNum>8?7:0) + (medNum>34?6:0));
				medNum++; // 49 is the ASCII value of '1', 7 jumps from ';' to 'A', 6 jumps from '[' to 'a'
			}
	}
	
	{
		QVector <char> mz(phant.nz, 0);
		QVector <QVector <char> > my(phant.ny, mz);
		QVector <QVector <QVector <char> > > mx(phant.nx, my);
		phant.m = mx;
		QVector <double> dz(phant.nz, 0);
		QVector <QVector <double> > dy(phant.ny, dz);
		QVector <QVector <QVector <double> > > dx(phant.nx, dy);
		phant.d = dx;
	}
	
	// Define xy bound values, still assuming first slice matches the rest
    for (int i = 0; i <= phant.nx; i++)
		phant.x[i] = (imagePos[0][0]+(i-0.5)*xySpacing[0][0])/10.0;
    for (int i = 0; i <= phant.ny; i++)
		phant.y[i] = (imagePos[0][1]+(i-0.5)*xySpacing[0][1])/10.0;

	// Define z bound values
	double prevZ, nextZ;
	nextZ = imagePos[0][2]-zSpacing[0]/2.0;
    for (int i = 0; i < phant.nz; i++) {
		prevZ = nextZ/2.0 + (imagePos[i][2]-zSpacing[i]/2.0)/2.0;
		nextZ = imagePos[i][2]+zSpacing[i]/2.0;
		phant.z[i] = prevZ/10.0;
	}
	phant.z.last() = nextZ/10.0;
		
	// ---------------------------------------------------------- //
	// CONVERTING HU TO APPROPRIATE DENSITY AND MEDIUM            //
	// ---------------------------------------------------------- //
	/*
	In this section, we define a list of bounding rectangles for
	each structure, which will be used to determine whether or not
	we need to check against structures when going through
	different slices and columns for efficiency.
	
	If masks are requested, they are generated based on the
	egsphant dimensions.
	
	Here is the main loop, that iterates through z, y, and x HU
	data.  For every new z value, a list of structures that could
	occur on this slice is made.  Then, for each new y, the list
	is checked again to see if a struct occurs in this column.
	
	Then, the HU data for the [i,j,k] voxel is used to assign
	media and density to the [i,ny-1-j,k] voxel of the egsphant.
	This inverts the y-axis, so that it counts up from -y rather
	than from +y.
	
	Finally, the egsphant file and masks if requested are output.
	*/
	
	// Get bounding rectangles over each struct
	QVector <QVector <QRectF> > structRect;
	
	for (int i = 0; i < structPos.size(); i++) {
		structRect.resize(i+1);
		for (int j = 0; j < structPos[i].size(); j++) {
			structRect[i].resize(j+1);
			structRect[i][j] = structPos[i][j].boundingRect();
		}
	}
	
	// Setup masks if makeMasks is set
	QVector <EGSPhant*> masks;
	if (makeMasks && !structName.isEmpty()) {
		for (int i = 0; i < structName.size(); i++) {
			EGSPhant* temp = new EGSPhant;
			temp->makeMask(&phant);
			masks << temp;
		}
	}
	
	// Arrays that hold the struct numbers and center voxel values to be used
	QList<QPoint> zIndex, yIndex;
	QList<QPoint>::iterator p;
	double zMid, yMid, xMid;
	int tempHU = 0, n = 0, q = 0, inStruct = 0, prio = 0, nj = 0;
	
	// Convert HU to density and media without masks
	for (int k = 0; k < phant.nz; k++) { // Z //
		if (structZ.size() > 0) {
			zIndex.clear(); // Reset lookup
			zMid = (phant.z[k]+phant.z[k+1])/2.0;
			for (int l = 0; l < structZ.size(); l++)
				if (structUnique[l])
					for (int m = 0; m < structZ[l].size(); m++) {
						// If slice j of struct i on the same plane as slice k of the phantom
						if (abs(structZ[l][m] - zMid) < (phant.z[k+1]-phant.z[k])/2.0) {
							zIndex << QPoint(l,m); // Add it to lookup
						}
					}
		}
		
		for (int j = 0; j < phant.ny; j++) { // Y //
			if (zIndex.size() > 0) {
				yIndex.clear(); // Reset lookup
				yMid = (phant.y[j]+phant.y[j+1])/2.0;
				for (p = zIndex.begin(); p != zIndex.end(); p++) {
					// If column p->y() of struct p->x() on the same column as slice k,j of the phantom
					if (structRect[p->x()][p->y()].top() <= yMid && yMid <= structRect[p->x()][p->y()].bottom()) {
						yIndex << *p;
					}
				}
			}
			
			for (int i = 0; i < phant.nx; i++) { // X //
				tempHU = HU[k][j][i];
				nj = phant.ny-1-j; // Reversed j index for density and media assignment
				xMid = (phant.x[i]+phant.x[i+1])/2.0;
				
				// Linear search because I don't think these arrays every get big
				// get the right density
				for (n = 0; n < HUMap.size()-1; n++)
					if (HUMap[n] <= tempHU && tempHU < HUMap[n+1])
						break;
				if (tempHU < HUMap[0])
					n = 0;
					
				double temp = interp(tempHU,HUMap[n],HUMap[n+1],denMap[n],denMap[n+1]);
				temp = temp<=0?0.000001:temp; // Set min density to 0.000001
				
				if (temp > phant.maxDensity) // Track max density for images
					phant.maxDensity = temp;
				
				phant.d[i][nj][k] = temp;
				
				// get the right media
				if (yIndex.size() > 0) {
					inStruct = 0;
					prio = 0;
					for (p = yIndex.begin(); p != yIndex.end(); p++) { // Check through each
					// If row p->y() of struct p->x() on the same row as slice k,j,i of the phantom
						if (structRect[p->x()][p->y()].left() <= xMid && xMid <= structRect[p->x()][p->y()].right())
							if (structPos[p->x()][p->y()].containsPoint(QPointF(xMid,yMid), Qt::OddEvenFill))
								if (structPrio[p->x()] > prio) {
									inStruct = p->x()+1; // Inflate index for the next check
									prio = structPrio[p->x()]; // Set priority at this struct's prio
								}
					}
				}
				
				if (inStruct) { 
					inStruct--; // Deflate index again
					q = structLookup[structReference[inStruct]]; // get the structName index which matches denThresholds
					for (n = 0; n < denThresholds[q].size()-1; n++)
						if (temp < denThresholds[q][n])
							break;
					
					phant.m[i][phant.ny-1-j][k] = mediaMap[medThresholds[q][n]];
					if (makeMasks)
						masks[inStruct]->m[i][nj][k] = 1;
				}
				else {
					for (n = 0; n < denThreshold.size()-1; n++)
						if (temp < denThreshold[n])
							break;
						
					phant.m[i][nj][k] = 49 + n + (n>8?7:0) + (n>34?6:0);
				}
				
				
				if (nominalDensity) {
					n = phant.m[i][nj][k] - 49;
					n = n>8?n-7:n;
					n = n>26?n-6:n;
					for (int l = 0; l < medNom.size(); l++)
						if (!medNom[l].compare(phant.media[n]))
							phant.d[i][nj][k] = denNom[l];
				}
			}
		}
	}
	
	duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
    std::cout << "Succesfully generated egsphant (dimensions x: [" << phant.x[0] << "," << phant.x[phant.nx] << "], y:["
			  << phant.y[0] << "," << phant.y[phant.ny] << "], z:["
			  << phant.z[0] << "," << phant.z[phant.nz] << "]).  Time elapsed is " << duration << " s.\n";
	
	// Save file
	phant.saveEGSPhantFile("PrimaryOutput.egsphant");
	//phant.savebEGSPhantFile("PrimaryOutput.begsphant");
	duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
    std::cout << "File successfully output.  Time elapsed is " << duration << " s.\n";
	
	// Output and delete masks
	if (makeMasks && !structName.isEmpty()) {
		for (int i = masks.size()-1; i >= 0; i--) {
			masks[i]->saveEGSPhantFile(structName[i]+"_mask.egsphant");
			delete masks[i];
		}
		masks.clear();
		duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
		std::cout << "Masks successfully output.  Time elapsed is " << duration << " s.\n";
	}

	// ---------------------------------------------------------- //
	// OUTPUT IMAGES                                              //
	// ---------------------------------------------------------- //
	/*
	In this section, we simply output images using the prebuilt
	image functions in the egsphant class.  There is also an
	additional step in drawing the structure points overtop of the
	media images.
	*/
	
	double xi = (phant.x[0]+phant.x[1])/2.0;
	double xf = (phant.x[phant.nx-1]+phant.x[phant.nx])/2.0;
	double yi = (phant.y[0]+phant.y[1])/2.0;
	double yf = (phant.y[phant.ny-1]+phant.y[phant.ny])/2.0;
	double z;
	double res = 10.0/xySpacing[0][0]*2; // This sets resolution to be 2 pixels for each voxel in x
	QImage temp;
	QList <QPointF> tempF;
	QList<QPointF>::iterator tempIt;
	QPen pen;
	pen.setWidth(2);
	if (outputImages) {
		for (int i = 0; i < phant.nz; i++) {
			z = imagePos[i][2]/10.0;
			phant.getEGSPhantPicDen("z axis", yi, yf, xi, xf, z, res).save(QString("Image/DenPic")+QString::number(i+1)+".png");
			
			temp = phant.getEGSPhantPicMed("z axis", yi, yf, xi, xf, z, res);
			QPainter paint (&temp);
			for (int j = 0; j < structZ.size(); j++)
				for (int k = 0; k < structZ[j].size(); k++) {
					zMid = (phant.z[i]+phant.z[i+1])/2.0;
					pen.setColor(QColor(double(j)/double(structZ.size())*255.0,0,255.0-double(j)/double(structZ.size())*255.0));
					paint.setPen(pen);
					if (abs(structZ[j][k] - zMid) < (phant.z[i+1]-phant.z[i])/2.0) {
						tempF = structPos[j][k].toList();
						for (tempIt = tempF.begin(); tempIt != tempF.end(); tempIt++)
							paint.drawPoint(phant.getIndex("x axis", tempIt->x())*2, phant.getIndex("y axis", tempIt->y())*2);
					}
				}
			temp.save(QString("Image/MedPic")+QString::number(i+1)+".png");
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