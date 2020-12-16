#include "dose.h"

int main(int argc, char **argv) {
	// Start clock for timing
    std::clock_t start;
    double duration;
    start = std::clock();
	
	if (argc == 1) {
        std::cout << "Please call this program with a .3ddose or .b3ddose file.\n";
        return 0;
    }
	
	Dose output;
	double scalingFactor = 0;
	QString name = "egs_mird";
	
    for (int i = 0; i < argc-1; i++) {
        QString path(argv[i+1]);
		
        if (!path.left(8).compare("scaling=")) {
			scalingFactor = (path.right(path.size()-8)).toDouble();
		}
        else if (!path.left(5).compare("name=")) {
			name = (path.right(path.size()-5));
		}
		else if (!path.right(7).compare(".3ddose")) {
			output.readIn(path,1);
		}
		else if (!path.right(8).compare(".b3ddose")) {
			output.readBIn(path,1);
		}
		else {
			std::cout << "DICOM_from_3ddose invoked with arguments which are not \"scaling=X\", \"name=Y\",  and a 3ddose file, exiting.\n";
			return 0;
		}
    }
	
	if (!output.x || !output.y || !output.z) {
		std::cout << "Did not succesfully read in 3ddose file, exiting.\n";
		return 0;
	}
	if (scalingFactor <= 0) {
		std::cout << "Did not read in positive, non-zero value for scaling factor, exiting.\n";
		return 0;
	}
	output.scale(scalingFactor);
	
	duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
	std::cout << "Read the dose file and scaled it by " << scalingFactor << ".  Time elapsed is " << duration << " s.\n";
	
	// Open the DICOM file
    QFile file(name+".RTDOSE.dcm");
    if (file.open(QIODevice::WriteOnly)) {
        unsigned char *dat;
        QDataStream out(&file);
        out.setByteOrder(QDataStream::LittleEndian);
		
		// Output the header, can be anything, needs to be 128 characters long
		dat = new unsigned char[128];
		strncpy((char*)dat,"     This file is written by the DICOM_from_3ddose application intended for use with egs_mird, written by Martin Martinov.      ",128);
		out.writeRawData((char*)dat, 128);
		delete[] dat;
		
		
		// Output the DICM header
		dat = new unsigned char[4];
		strncpy((char*)dat,"DICM",4);
		out.writeRawData((char*)dat, 4);
		delete[] dat;
		
		duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
		std::cout << "Finished outputting header.  Time elapsed is " << duration << " s.\n";
		
		int size;
		// Below are data elements taken from an example, some of which are used for output
		
		//(Group,Element)  	TAG Description                    	VR	VM	Length    	Value
		//(0002,0000)      	FileMetaInformationGroupLength     	UL	1	4         	192
		size = 8 + 4;
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x2; dat[3] = 0x0; dat[2] = 0x0;
		strncpy((char*)(&dat[4]),"UL",2); dat[6] = size-8; dat[7] = 0;
		//dat[8] = 192; dat[9] = 0; dat[10] = 0; dat[11] = 0;
		dat[8] = 72; dat[9] = 0; dat[10] = 0; dat[11] = 0;
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0002,0001)      	FileMetaInformationVersion         	OB	1	2 		
		// Weird exception //////////////////////////////////////////////////
		size = 8 + 6;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x2; dat[3] = 0x0; dat[2] = 0x1;
		strncpy((char*)(&dat[4]),"OB",2); dat[6] = 0; dat[7] = 0;
		dat[8] = 2; dat[9] = 0; dat[10] = 0; dat[11] = 0;
		dat[12] = 0; dat[13] = 1;
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0002,0002)      	MediaStorageSOPClassUID            	UI	1	30        	1.2.840.10008.5.1.4.1.1.481.2
		size = 8 + 30;
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x2; dat[3] = 0x0; dat[2] = 0x2;
		strncpy((char*)(&dat[4]),"UI",2); dat[6] = size-8; dat[7] = 0;
		strncpy((char*)(&dat[8]),"1.2.840.10008.5.1.4.1.1.481.2",size-8); // RT Dose Storage
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0002,0003)      	MediaStorageSOPInstanceUID         	UI	1	58        	2.16.840.1.114362.1.11753351.23780827009.551244094.509.207
		//size = 8 + 58;  
		//dat = new unsigned char[size];
		//
		//dat[1] = 0x0; dat[0] = 0x2; dat[3] = 0x0; dat[2] = 0x3;
		//strncpy((char*)(&dat[4]),"UI",2); dat[6] = size-8; dat[7] = 0;
		//strncpy((char*)(&dat[8]),"2.16.840.1.114362.1.11753351.23780827009.551244094.509.207",size-8);
		//
		//out.writeRawData((char*)dat, size);
		//delete[] dat;
		
		//(0002,0010)      	TransferSyntaxUID                  	UI	1	18        	1.2.840.10008.1.2
		size = 8 + 18;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x2; dat[3] = 0x0; dat[2] = 0x10;
		strncpy((char*)(&dat[4]),"UI",2); dat[6] = size-8; dat[7] = 0;
		strncpy((char*)(&dat[8]),"1.2.840.10008.1.2",size-8); // Implicit VR Endian: Default Transfer Syntax for DICOM
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0002,0012)      	ImplementationClassUID             	UI	1	20        	2.16.840.1.114362.1
		//size = 8 + 20;  
		//dat = new unsigned char[size];
		//
		//dat[1] = 0x0; dat[0] = 0x2; dat[3] = 0x0; dat[2] = 0x12;
		//strncpy((char*)(&dat[4]),"UI",2); dat[6] = size-8; dat[7] = 0;
		//strncpy((char*)(&dat[8]),"2.16.840.1.114362.1",size-8);
		//
		//out.writeRawData((char*)dat, size);
		//delete[] dat;
		
		//(0002,0013)      	ImplementationVersionName          	SH	1	12        	MIM683I53007
		//size = 8 + 12;  
		//dat = new unsigned char[size];
		//
		//dat[1] = 0x0; dat[0] = 0x2; dat[3] = 0x0; dat[2] = 0x13;
		//strncpy((char*)(&dat[4]),"SH",2); dat[6] = size-8; dat[7] = 0;
		//strncpy((char*)(&dat[8]),"MIM683I53007",size-8);
		//
		//out.writeRawData((char*)dat, size);
		//delete[] dat;
		
		
		duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
		std::cout << "Output group 0002 elements, swapping to transfer syntax.  Time elapsed is " << duration << " s.\n";
		// Explicit Value Representation stops here
		
		//(0008,0005)      	SpecificCharacterSet               	CS	1	10        	ISO_IR 100
		size = 8 + 10;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x8; dat[3] = 0x0; dat[2] = 0x5;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"ISO_IR 100",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0008,0012)      	InstanceCreationDate               	DA	1	8         	20200629
		auto now = std::chrono::system_clock::now();
		std::time_t now_c = std::chrono::system_clock::to_time_t(now);
		struct tm *parts = std::localtime(&now_c);
		
		std::string date = std::to_string(1900+parts->tm_year);
		std::string temp = std::to_string(1+parts->tm_mon);
		if (temp.size() == 1) temp = "0"+temp;
		date = date + temp;
		temp=std::to_string(parts->tm_mday);
		if (temp.size() == 1) temp = "0"+temp;
		date = date + temp;
		
		size = 8 + date.size();  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x8; dat[3] = 0x0; dat[2] = 0x12;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)&(dat[8]),date.data(),size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0008,0013)      	InstanceCreationTime               	TM	1	6         	161448
		std::string time = std::to_string(parts->tm_hour);
		if (time.size() == 1) time = "0"+time;
		temp = std::to_string(parts->tm_min);
		if (temp.size() == 1) temp = "0"+temp;
		time = time + temp;
		temp = std::to_string(parts->tm_sec);
		if (temp.size() == 1) temp = "0"+temp;
		time = time + temp;
		
		size = 8 + time.size();  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x8; dat[3] = 0x0; dat[2] = 0x13;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)&(dat[8]),time.data(),size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0008,0016)      	SOPClassUID                        	UI	1	30        	1.2.840.10008.5.1.4.1.1.481.2
		size = 8 + 30;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x8; dat[3] = 0x0; dat[2] = 0x16;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"1.2.840.10008.5.1.4.1.1.481.2",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0008,0018)      	SOPInstanceUID                     	UI	1	58        	2.16.840.1.114362.1.11753351.23780827009.551244094.509.207
		//size = 8 + 58;  
		//dat = new unsigned char[size];
		//
		//dat[1] = 0x0; dat[0] = 0x8; dat[3] = 0x0; dat[2] = 0x18;
		//dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		//strncpy((char*)(&dat[8]),"2.16.840.1.114362.1.11753351.23780827009.551244094.509.20",size-8);
		//
		//out.writeRawData((char*)dat, size);
		//delete[] dat;
		
		//(0008,0020)      	StudyDate                          	DA	1	8         	20150728
		//(0008,0023)      	ContentDate                        	DA	1	8         	20200629
		//(0008,0030)      	StudyTime                          	TM	1	6         	161008
		//(0008,0033)      	ContentTime                        	TM	1	6         	161448
		//(0008,0050)      	AccessionNumber                    	SH	0	0         	
		size = 8;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x8; dat[3] = 0x0; dat[2] = 0x50;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0008,0060)      	Modality                           	CS	1	6         	RTDOSE
		size = 8+6;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x8; dat[3] = 0x0; dat[2] = 0x60;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"RTDOSE",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0008,0070)      	Manufacturer                       	LO	1	4         	ADAC
		//(0008,0080)      	InstitutionName                    	LO	0	0         	
		//(0008,0090)      	ReferringPhysicianName             	PN	0	0         	
		//(0008,1030)      	StudyDescription                   	LO	1	6         	LUNR1
		size = 8+20;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x8; dat[3] = 0x10; dat[2] = 0x30;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"EGS MIRD CALCULATION",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0008,103E)      	SeriesDescription                  	LO	1	10        	Loc: LUNR1
		//(0008,1090)      	ManufacturerModelName              	LO	1	10        	Pinnacle3
		//(0008,1110)      	ReferencedStudySequence            	SQ	1	-1        	
		//(0008,1150)   	ReferencedSOPClassUID              	UI	1	24        	1.2.840.10008.3.1.2.3.2
		size = 8 + 24;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x8; dat[3] = 0x11; dat[2] = 0x50;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"1.2.840.10008.3.1.2.3.2",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0008,1155)   	ReferencedSOPInstanceUID           	UI	1	58        	2.16.840.1.114362.1.11753351.23780827009.551244090.852.92
		//size = 8 + 58;  
		//dat = new unsigned char[size];
		//
		//dat[1] = 0x0; dat[0] = 0x8; dat[3] = 0x11; dat[2] = 0x55;
		//dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		//strncpy((char*)(&dat[8]),"2.16.840.1.114362.1.11753351.23780827009.551244090.852.92",size-8);
		//
		//out.writeRawData((char*)dat, size);
		//delete[] dat;
		
		//(0010,0010)      	PatientName                        	PN	1	10        	MISSILE_8
		//(0010,0020)      	PatientID                          	LO	1	10        	MISSILE_8
		//(0010,0030)      	PatientBirthDate                   	DA	0	0         	
		//(0010,0040)      	PatientSex                         	CS	1	2         	F
		//(0012,0062)      	PatientIdentityRemoved             	CS	1	4         	YES
		//(0012,0063)      	DeidentificationMethod             	LO	1	36        	Limited Data Set: MIM.6.8.3.I530-07
		//(0018,0050)      	SliceThickness                     	DS	1	2         	3
		std::string thick = std::to_string((output.cz[1]-output.cz[0])*10);
		if (!(thick.size()%2)) thick = thick+" ";
		size = 8 + thick.size();  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x18; dat[3] = 0x0; dat[2] = 0x50;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)&(dat[8]),thick.data(),size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0018,1020)      	SoftwareVersions                   	LO	2	8         	9.8\9.8
		//(0020,000D)      	StudyInstanceUID                   	UI	1	58        	2.16.840.1.114362.1.11753351.23780827009.551244090.852.92
		//size = 8 + 58;  
		//dat = new unsigned char[size];
		//
		//dat[1] = 0x0; dat[0] = 0x20; dat[3] = 0x0; dat[2] = 0xD;
		//dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		//strncpy((char*)(&dat[8]),"2.16.840.1.114362.1.11753351.23780827009.551244090.852.92",size-8);
		//
		//out.writeRawData((char*)dat, size);
		//delete[] dat;
		
		//(0020,000E)      	SeriesInstanceUID                  	UI	1	58        	2.16.840.1.114362.1.11753351.23780827009.551244094.509.208
		//size = 8 + 58;  
		//dat = new unsigned char[size];
		//
		//dat[1] = 0x0; dat[0] = 0x20; dat[3] = 0x0; dat[2] = 0xE;
		//dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		//strncpy((char*)(&dat[8]),"2.16.840.1.114362.1.11753351.23780827009.551244094.509.208",size-8);
		//
		//out.writeRawData((char*)dat, size);
		//delete[] dat;
		
		//(0020,0010)      	StudyID                            	SH	1	10        	MISSILE_8
		size = 8+20;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x20; dat[3] = 0x0; dat[2] = 0x10;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"EGS MIRD CALCULATION",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0020,0011)      	SeriesNumber                       	IS	1	2         	1
		size = 8+2;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x20; dat[3] = 0x0; dat[2] = 0x11;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"1 ",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0020,0013)      	InstanceNumber                     	IS	1	2         	1
		size = 8+2;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x20; dat[3] = 0x0; dat[2] = 0x13;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"1 ",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0020,0032)      	ImagePositionPatient               	DS	3	26        	-254.665\-110.579\-705.34
		std::string position = std::to_string((output.cx[1]+output.cx[0])/2*10)+"\\"+
							   std::to_string((output.cy[1]+output.cy[0])/2*10)+"\\"+
							   std::to_string((output.cz[1]+output.cz[0])/2*10);
		if (!(position.size()%2)) position = position+" ";
		size = 8 + position.size();  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x20; dat[3] = 0x0; dat[2] = 0x32;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)&(dat[8]),position.data(),size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0020,0037)      	ImageOrientationPatient            	DS	6	12        	1\0\0\0\1\0
		size = 8 + 12;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x20; dat[3] = 0x0; dat[2] = 0x37;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"1\\0\\0\\0\\1\\0 ",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0020,0052)      	FrameOfReferenceUID                	UI	1	58        	2.16.840.1.114362.1.11753351.23780827009.551244090.868.204
		//size = 8 + 58;  
		//dat = new unsigned char[size];
		//
		//dat[1] = 0x0; dat[0] = 0x20; dat[3] = 0x0; dat[2] = 0x52;
		//dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		//strncpy((char*)(&dat[8]),"2.16.840.1.114362.1.11753351.23780827009.551244090.868.204",size-8);
		//
		//out.writeRawData((char*)dat, size);
		//delete[] dat;
		
		//(0020,1040)      	PositionReferenceIndicator         	LO	0	0         	
		size = 8 + 0;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x20; dat[3] = 0x10; dat[2] = 0x40;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0028,0002)      	SamplesPerPixel                    	US	1	2         	1
		size = 8 + 2;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x28; dat[3] = 0x0; dat[2] = 0x2;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		dat[8] = 1; dat[9] = 0;
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0028,0004)      	PhotometricInterpretation          	CS	1	12        	MONOCHROME2
		size = 8+12;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x28; dat[3] = 0x0; dat[2] = 0x4;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"MONOCHROME2 ",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0028,0008)      	NumberOfFrames                     	IS	1	2         	93
		std::string zCount = std::to_string(output.z);
		if (!(zCount.size()%2)) zCount = zCount+" ";
		size = 8+zCount.size();  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x28; dat[3] = 0x0; dat[2] = 0x08;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)&(dat[8]),zCount.data(),size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0028,0009)      	FrameIncrementPointer              	AT	1	4         	(3004,000c)
		size = 8+4;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x28; dat[3] = 0x0; dat[2] = 0x9;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		dat[9] = 0x30; dat[8] = 0x4; dat[11] = 0x0; dat[10] = 0xc;
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0028,0010)      	Rows                               	US	1	2         	115
		size = 8+2;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x28; dat[3] = 0x0; dat[2] = 0x10;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		dat[8] = output.x%(1<<8); dat[9] = int(output.x/(1<<8));
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0028,0011)      	Columns                            	US	1	2         	165
		size = 8+2;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x28; dat[3] = 0x0; dat[2] = 0x11;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		dat[8] = output.y%(1<<8); dat[9] = int(output.y/(1<<8));
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
			
		//(0028,0030)      	PixelSpacing                       	DS	2	4         	3\3
		std::string xyThick = std::to_string((output.cx[1]-output.cx[0])*10)+"\\"+
							  std::to_string((output.cy[1]-output.cy[0])*10)+" ";
		if (!(xyThick.size()%2)) xyThick = xyThick+" ";
		size = 8+xyThick.size();
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x28; dat[3] = 0x0; dat[2] = 0x30;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)&(dat[8]),xyThick.data(),size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0028,0100)      	BitsAllocated                      	US	1	2         	16
		size = 8+2;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x28; dat[3] = 0x1; dat[2] = 0x0;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		dat[8] = 16; dat[9] = 0;
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0028,0101)      	BitsStored                         	US	1	2         	16
		size = 8+2;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x28; dat[3] = 0x1; dat[2] = 0x1;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		dat[8] = 16; dat[9] = 0;
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0028,0102)      	HighBit                            	US	1	2         	15
		size = 8+2;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x28; dat[3] = 0x1; dat[2] = 0x2;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		dat[8] = 15; dat[9] = 0;
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(0028,0103)      	PixelRepresentation                	US	1	2         	0
		size = 8+2;  
		dat = new unsigned char[size];
		
		dat[1] = 0x0; dat[0] = 0x28; dat[3] = 0x1; dat[2] = 0x3;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		dat[8] = 0; dat[9] = 0;
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(3004,0002)      	DoseUnits                          	CS	1	2         	GY
		size = 8+2;  
		dat = new unsigned char[size];
		
		dat[1] = 0x30; dat[0] = 0x4; dat[3] = 0x0; dat[2] = 0x2;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"GY",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(3004,0004)      	DoseType                           	CS	1	8         	PHYSICAL
		size = 8+8;  
		dat = new unsigned char[size];
		
		dat[1] = 0x30; dat[0] = 0x4; dat[3] = 0x0; dat[2] = 0x4;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"PHYSICAL",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(3004,000A)      	DoseSummationType                  	CS	1	8         	FRACTION
		size = 8+6;  
		dat = new unsigned char[size];
		
		dat[1] = 0x30; dat[0] = 0x4; dat[3] = 0x0; dat[2] = 0xA;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"RECORD",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(3004,000C)      	GridFrameOffsetVector              	DS	93	334       	0\3\6\9\12\15\18\21\24\27\30\33\36\39\42\45\48\51\54\57\60\63\66\69\72\75\78\81\84\87\90\93\96\99\102\105\108\111\114\117\120\123\126\129\132\135\138\141\144\147\150\153\156\159\162\165\168\171\174\177\180\183\186\189\192\195\198\201\204\207\210\213\216\219\222\225\228\231\234\237\240\243\246\249\252\255\258\261\264\267\270\273\276
		std::string zPlanes = "";
		for (int i = 0; i < output.z-1; i++)
			zPlanes += std::to_string((output.cz[i+1]+output.cz[i])/2*10)+"\\";
		zPlanes += std::to_string((output.cz[output.z+1]+output.cz[output.z])/2*10);
		if (!(zPlanes.size()%2)) zPlanes = zPlanes+" ";
		size = 8 + zPlanes.size();
		dat = new unsigned char[size];
		
		dat[1] = 0x30; dat[0] = 0x4; dat[3] = 0x0; dat[2] = 0xC;
		dat[4] = (size-8)%(1<<8); dat[5] = int((size-8)/(1<<8)); dat[6] = int((size-8)/(1<<16))%(1<<8); dat[7] = int((size-8)/(1<<24))%(1<<8);
		strncpy((char*)&(dat[8]),zPlanes.data(),size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(3004,000E)      	DoseGridScaling                    	DS	1	10        	0.00510632
		double scaling = (0xEFFF)/output.getMax();
		output.scale(scaling);
		
		std::string scalingString = std::to_string(1/scaling);
		if (!(scalingString.size()%2)) scalingString = scalingString+" ";
		size = 8 + scalingString.size();
		dat = new unsigned char[size];
		
		dat[1] = 0x30; dat[0] = 0x4; dat[3] = 0x0; dat[2] = 0xE;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)&(dat[8]),scalingString.data(),size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		//(3004,0014)      	TissueHeterogeneityCorrection      	CS	1	6         	IMAGE
		size = 8+6;  
		dat = new unsigned char[size];
		
		dat[1] = 0x30; dat[0] = 0x4; dat[3] = 0x0; dat[2] = 0x14;
		dat[4] = size-8; dat[5] = 0; dat[6] = 0; dat[7] = 0;
		strncpy((char*)(&dat[8]),"IMAGE ",size-8);
		
		out.writeRawData((char*)dat, size);
		delete[] dat;
		
		duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
		std::cout << "Output all non-pixel data elements.  Time elapsed is " << duration << " s.\n";
		
		//(300C,0002)      	ReferencedRTPlanSequence           	SQ	1	-1        	
		//(0008,1150)   	ReferencedSOPClassUID              	UI	1	30        	1.2.840.10008.5.1.4.1.1.481.5
		//(0008,1155)   	ReferencedSOPInstanceUID           	UI	1	58        	2.16.840.1.114362.1.11753351.23780827009.551244092.853.205
		//(300C,0020)   	ReferencedFractionGroupSequence    	SQ	1	-1
		//(300C,0022)   	ReferencedFractionGroupNumber      	IS	1	2         	1
		//(7FE0,0010)      	PixelData                          	OW	1	37950
		dat = new unsigned char[8];
		unsigned long int bSize = output.x*output.y*output.z*2;
		
		dat[1] = 0x7F; dat[0] = 0xE0; dat[3] = 0x0; dat[2] = 0x10;
		dat[4] = bSize%(1<<8); dat[5] = int(bSize/(1<<8))%(1<<8); dat[6] = int(bSize/(1<<16))%(1<<8); dat[7] = int(bSize/(1<<24));
		
		out.writeRawData((char*)dat, 8);
		delete[] dat;
		
		unsigned short int bDat;
		for (int k = 0; k < output.z; k++)
			for (int j = output.y-1; j >= 0; j--)
				for (int i = 0; i < output.x; i++) {
					bDat = output.val[i][j][k];
					out << bDat;
				}
	}
	
	duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
	std::cout << "Successfully output DICOM file.  Time elapsed is " << duration << " s.\n";
    return 1;
}