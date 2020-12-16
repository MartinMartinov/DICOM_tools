#include "DICOM.h"

//#define OUTPUT_ALL
//#define OUTPUT_TAG
//#define OUTPUT_SQ
#define MAX_DATA_PRINT 0 // 0 means any size

Attribute::Attribute() {
    vf = NULL; // This stops seg faults when calling the destructor below
}

Attribute::~Attribute() {
    if (vf != NULL) {
        delete[] vf;
    }
}

SequenceItem::SequenceItem(unsigned long int size, unsigned char *data) {
    vl = size;
    vf = data;
}

SequenceItem::~SequenceItem() {
    if (vf != NULL) {
		delete[] vf;
    }
}

Sequence::~Sequence() {
    for (int i = 0; i < items.size(); i++) {
        delete items[i];
    }
    items.clear();
}

DICOM::DICOM(database *l) {
    lib = l;
    isImplicit = isBigEndian = false;
}

DICOM::~DICOM() {
    for (int i = 0; i < data.size(); i++) {
        delete data[i];
    }
    data.clear();
}

int DICOM::readSequence(QDataStream *in, Attribute *att) {
    bool flag = true;
	int depth = 0;
    unsigned int *tag, size;
    unsigned char *dat;
    while (flag) {
        tag = new unsigned int;
        if (in->readRawData((char*)tag, 4) != 4) {
            // Not a DICOM file
            delete tag;
            return 0;
        }
        if (in->readRawData((char*)(&size), 4) != 4) {
            // Not a DICOM file
            delete tag;
            return 0;
        }

        if (*tag == (unsigned int)0xE0DDFFFE) { // sequence delimiter
            flag = false;
        }
        else if (size != (unsigned int)0xFFFFFFFF) {
            // sequence item with defined size
            dat = new unsigned char[size];
            if ((unsigned int)in->readRawData((char*)dat, size) != size) {
                // Not a DICOM file
                delete tag;
                delete[] dat;
                return 0;
            }
            att->seq.items.append(new SequenceItem(size, dat));
        }
        else if (size == (unsigned int)0xFFFFFFFF) {
            // sequence item with undefined size
            QByteArray buffer;
			depth = 0;
            dat = new unsigned char[1];
            while (true) { // Keeping appending data one char at a time...
                if (in->readRawData((char*)dat, 1) != 1) {
                    // Not a DICOM file
                    delete tag;
                    delete[] dat;
                    return 0;
                }
                buffer.append(*dat);
				
				// Is there a subsequence we are parsing
				if (buffer.endsWith("����")) {
					if (!isImplicit) {
						if (buffer.endsWith("SQ\0\0����")) {
							depth++; // Increase depth to skip delimiters until we exit subsequence	
						}
					}
					else {
						unsigned short int tag[2];
						tag[0] = *(unsigned short int*)(buffer.right(8).left(2).data());
						tag[1] = *(unsigned short int*)(buffer.right(6).left(2).data());
						Reference nearest = lib->binSearch(tag[0], tag[1], 0, lib->lib.size()-1);
						if (nearest.tag[0] == tag[0] && nearest.tag[1] == tag[1] && !nearest.vr.compare("SQ")) {
							depth++; // Increase depth to skip delimiters until we exit subsequence								
						}
					}
				}
				
                // ...until we reach the sequence item delimiter
                if (*((unsigned int*)buffer.right(4).data()) == (unsigned int)0xE00DFFFE && !depth) {
					buffer.chop(4);
                    delete[] dat;
                    dat = new unsigned char[4];
                    if (in->readRawData((char*)dat, 4) != 4) {
                        // Not a DICOM file
                        delete tag;
                        delete[] dat;
                        return 0;
                    }
                    delete[] dat;
                    dat = new unsigned char[buffer.size()];
                    for (int i = 0; i < buffer.size(); i++) {
                        dat[i] = buffer[i];
                    }
                    att->seq.items.append(new SequenceItem(buffer.size(),dat));
                    break;
                }
				else if (*((unsigned int*)buffer.right(4).data()) == (unsigned int)0xE0DDFFFE) {
					depth--;
				}
            }
        }

        delete tag;
    }
    return 1;
}

int DICOM::readDefinedSequence(QDataStream *in, Attribute *att, unsigned long int n) {
	int depth = 0;
    unsigned int *tag, size;
    unsigned char *dat;
    while (n > 0) {
        tag = new unsigned int;
        if (in->readRawData((char*)tag, 4) != 4) {
            // Not a DICOM file
            delete tag;
            return 0;
        }
        if (in->readRawData((char*)(&size), 4) != 4) {
            // Not a DICOM file
            delete tag;
            return 0;
        }
		n-=8;

        if (size != (unsigned int)0xFFFFFFFF) {
            // sequence item with defined size
            dat = new unsigned char[size];
            if ((unsigned int)in->readRawData((char*)dat, size) != size) {
                // Not a DICOM file
                delete tag;
                delete[] dat;
                return 0;
            }
            att->seq.items.append(new SequenceItem(size, dat));
			n-=size;
        }
        else if (size == (unsigned int)0xFFFFFFFF) {
            // sequence item with undefined size
            QByteArray buffer;
            dat = new unsigned char[1];
            while (true) { // Keeping appending data one char at a time...
                if (in->readRawData((char*)dat, 1) != 1) {
                    // Not a DICOM file
                    delete tag;
                    delete[] dat;
                    return 0;
                }
                buffer.append(*dat);
				
				// Is there a subsequence we are parsing
				if (buffer.endsWith("����")) {
					if (!isImplicit) {
						if (buffer.endsWith("SQ\0\0����")) {
							depth++; // Increase depth to skip delimiters until we exit subsequence	
						}
					}
					else {
						unsigned short int tag[2];
						tag[0] = *(unsigned short int*)(buffer.right(8).left(2).data());
						tag[1] = *(unsigned short int*)(buffer.right(6).left(2).data());
						Reference nearest = lib->binSearch(tag[0], tag[1], 0, lib->lib.size()-1);
						if (nearest.tag[0] == tag[0] && nearest.tag[1] == tag[1] && !nearest.vr.compare("SQ")) {
							depth++; // Increase depth to skip delimiters until we exit subsequence								
						}
					}
				}
				
                // ...until we reach the sequence item delimiter
                if (*((unsigned int*)buffer.right(4).data()) == (unsigned int)0xE00DFFFE && !depth) {
                    delete[] dat;
                    dat = new unsigned char[4];
                    if (in->readRawData((char*)dat, 4) != 4) {
                        // Not a DICOM file
                        delete tag;
                        delete[] dat;
                        return 0;
                    }
                    delete[] dat;
                    dat = new unsigned char[buffer.size()];
                    for (int i = 0; i < buffer.size(); i++) {
                        dat[i] = buffer[i];
                    }
                    att->seq.items.append(new SequenceItem(buffer.size(),dat));
					n-=buffer.size();
                    break;
                }
				else if (*((unsigned int*)buffer.right(4).data()) == (unsigned int)0xE0DDFFFE)
					depth--;
            }
        }

        delete tag;
    }
    return 1;
}
	
int DICOM::parse(QString p) {
	path = p;
    QFile file(path);
    int k = 0, l = 0;
    if (file.open(QIODevice::ReadOnly)) {
        unsigned char *dat;
        QDataStream in(&file);
        in.setByteOrder(QDataStream::LittleEndian);

        /*============================================================================*/
        /*DICOM HEADER READER=========================================================*/
        // Skip the first bit of white space in DICOM
        dat = new unsigned char[128];
        if (in.readRawData((char*)dat, 128) != 128) {
            // Not a DICOM file
            delete[] dat;
            file.close();
            return 0;
        }
        delete[] dat;

        // Read in DICM characters at start of file
        dat = new unsigned char[4];
        if (in.readRawData((char*)dat, 4) != 4) {
            // Not a DICOM file
            delete[] dat;
            file.close();
            return 0;
        }
        else if ((QString(dat[0])+dat[1]+dat[2]+dat[3]) != "DICM") {
            // Not a DICOM file
            delete[] dat;
            file.close();
            return 0;
        }
        delete[] dat;

        /*============================================================================*/
        /*BEGINNING OF DATA ELEMENT READING LOOP======================================*/
        Attribute *temp;
        unsigned int size;
        QString VR;
        bool nested = false;
        while (!in.atEnd()) {
            temp = new Attribute();

            /*============================================================================*/
            /*RETRIEVE ELEMENT TAG========================================================*/
            // Get the tag
			k++; // iterate
			#if defined(OUTPUT_ALL) || defined(OUTPUT_TAG)
                std::cout << std::dec << k << ") " << "Tag "; 
            #endif
            dat = new unsigned char[4];
            if (in.readRawData((char*)dat,4) != 4) {
                // Not a DICOM file
                delete[] dat;
                file.close();
                return 0;
            }
            temp->tag[0]= ((unsigned short int)(dat[1]) << 8) +
                          (unsigned short int)dat[0];
            temp->tag[1]= ((unsigned short int)(dat[3]) << 8) +
                          (unsigned short int)dat[2];
			#if defined(OUTPUT_ALL) || defined(OUTPUT_TAG)
                std::cout << std::hex << temp->tag[0] << ","
                          <<  temp->tag[1] << " | Representation ";
			#endif
            delete[] dat;
			
			if (temp->tag[0] == 0xFFFE && (temp->tag[1] == 0xE0DD || temp->tag[1] == 0xE00D)) {
                // Not a DICOM file
				std::cout << "Misreading sequence delimiters as top level data elements, something has gone wrong, quitting...\n";
                delete[] dat;
                file.close();
                return 0;
			}

            /*============================================================================*/
            /*NON-NESTED PROCEDURE: GET VR, SIZE AND DATA=================================*/
            // Normal data elements
            // Get the VR
            dat = new unsigned char[4];
			
            if (!isImplicit || temp->tag[0] == 0x0002) {
                if (in.readRawData((char*)dat,4) != 4) {
                    // Not a DICOM file
                    delete[] dat;
                    file.close();
                    return 0;
                }

                VR = QString(dat[0])+dat[1];
					
				#if defined(OUTPUT_ALL) || defined(OUTPUT_TAG)
					std::cout << ((unsigned short int)(dat[0]) << 8) +
							  (unsigned short int)dat[1]
							  << " -> " << VR.toStdString() << " | Size ";
				#endif
            }
            else {
                VR = lib->binSearch(temp->tag[0], temp->tag[1], 0, lib->lib.size()-1).vr;
				#if defined(OUTPUT_ALL) || defined(OUTPUT_TAG)
					std::cout << VR.toStdString() << " (implicit) | Size ";
				#endif
            }
			
            // Get size
			if ((temp->tag[0] != 0x0002 && isImplicit) || (lib->implicitVR.contains(VR))) {
                if (in.readRawData((char*)dat,4) != 4) { //Reread for size
                    // Not a DICOM file
                    delete[] dat;
                    file.close();
                    return 0;
                }
                temp->vl = ((unsigned int)(dat[3]) << 24) +
                           ((unsigned int)(dat[2]) << 16) +
                           ((unsigned int)(dat[1]) << 8) +
                           (unsigned int)dat[0];

                // We have a sequence
                if (!VR.compare("SQ") && temp->vl == (unsigned int)0xFFFFFFFF) {
                    nested = true;
                    if (!readSequence(&in, temp)) {
                        return 0;
                    }
                }
				else if (!VR.compare("SQ")) {
                    nested = true;
                    if (!readDefinedSequence(&in, temp, temp->vl)) {
                        return 0;
                    }
                }					
			}
            else {
				if (isImplicit && temp->tag[0] != 0x0002)
					if (in.readRawData((char*)dat,4) != 4) {
						// Not a DICOM file
						delete[] dat;
						file.close();
						return 0;
					}

				if (lib->validVR.contains(VR))
					temp->vl = ((unsigned short int)(dat[3]) << 8) +
							   (unsigned short int)dat[2];
				else
					temp->vl = ((unsigned int)(dat[3]) << 24) +
							   ((unsigned int)(dat[2]) << 16) +
							   ((unsigned int)(dat[1]) << 8) +
							   (unsigned int)dat[0];
							   
				// We have a sequence
                if (!VR.compare("SQ") && temp->vl == (unsigned int)0xFFFFFFFF) {
                    nested = true;
                    if (!readSequence(&in, temp)) {
                        return 0;
                    }
                }
				else if (!VR.compare("SQ")) {
                    nested = true;
                    if (!readDefinedSequence(&in, temp, temp->vl)) {
                        return 0;
                    }
                }	
			}

            if (temp->vl == (unsigned int)0xFFFFFFFF) {
                temp->vl = 0;
            }

            size = temp->vl;

			#if defined(OUTPUT_ALL) || defined(OUTPUT_TAG)
                std::cout << temp->vl << " -> " << std::dec << size << "\n";
			#endif
			
            Reference closest = lib->binSearch(temp->tag[0], temp->tag[1], 0, lib->lib.size()-1);
            if (closest.tag[0] == temp->tag[0] && closest.tag[1] == temp->tag[1]) {
                temp->desc = closest.title;
                l++;
            }
            else {
                temp->desc = "Unknown Tag";
            }

			#ifdef OUTPUT_ALL
                std::cout << temp->desc.toStdString() << ": ";
			#endif
            delete[] dat;

            // Get data
            if (!nested) {
                temp->vf = new unsigned char[size];
                if (size > 0 && size < (unsigned long int)INT_MAX) {
                    if (in.readRawData((char*)temp->vf,size) != (long int)size) {
                        // Not a DICOM file
                        file.close();
                        return 0;
                    }
                }
                else if (size > 0) { // In case we need to read in more data
                    // then the buffer can handle
                    unsigned char *pt = temp->vf;
                    for (int i = 0; (unsigned long int)i <
                            size/((unsigned long int)INT_MAX); i++) {
                        if (in.readRawData((char*)pt,INT_MAX) != INT_MAX) {
                            // Not a DICOM file
                            file.close();
                            return 0;
                        }
                        pt += sizeof(char)*INT_MAX;
                    }
                }

				#ifdef OUTPUT_ALL
                    unsigned long int avoidWarning =
                        (unsigned long int)MAX_DATA_PRINT;
                    if (avoidWarning == 0 || size < avoidWarning)
						// It's a string
						if (!VR.compare("UI") || !VR.compare("SH") || !VR.compare("AE") || !VR.compare("DA") ||
							!VR.compare("TM") || !VR.compare("LO") || !VR.compare("ST") || !VR.compare("PN") ||
							!VR.compare("DT") || !VR.compare("LT") || !VR.compare("UT") || !VR.compare("IS") ||
							!VR.compare("OW") || !VR.compare("DS") || !VR.compare("CS") || !VR.compare("AS"))
							for (unsigned long int i = 0; i < size; i++)
								std::cout << dat[i];
						// It's a tag
						else if (!VR.compare("AT"))
							std::cout << ((unsigned int)(dat[3]) << 24) +
										 ((unsigned int)(dat[2]) << 16) +
										 ((unsigned int)(dat[1]) << 8) +
										  (unsigned int)(dat[0]);
						else if (!VR.compare("FL"))
							if (isBigEndian)
								std::cout << std::dec << float(((int)(dat[0]) << 24) +
										 ((int)(dat[1]) << 16) +
										 ((int)(dat[2]) << 8) +
										  (int)(dat[3])) << std::hex;
							else
								std::cout << std::dec << float(((int)(dat[3]) << 24) +
										 ((int)(dat[2]) << 16) +
										 ((int)(dat[1]) << 8) +
										  (int)(dat[0])) << std::hex;
						else if (!VR.compare("FD"))
							if (isBigEndian)
								std::cout << std::dec << double(((long int)(dat[0]) << 56) +
										 ((long int)(dat[1]) << 48) +
										 ((long int)(dat[2]) << 40) +
										 ((long int)(dat[3]) << 32) +
										 ((long int)(dat[4]) << 24) +
										 ((long int)(dat[5]) << 16) +
										 ((long int)(dat[6]) << 8) +
										  (long int)(dat[7])) << std::hex;
							else
								std::cout << std::dec << double(((long int)(dat[7]) << 56) +
										 ((long int)(dat[6]) << 48) +
										 ((long int)(dat[5]) << 40) +
										 ((long int)(dat[4]) << 32) +
										 ((long int)(dat[3]) << 24) +
										 ((long int)(dat[2]) << 16) +
										 ((long int)(dat[1]) << 8) +
										  (long int)(dat[0])) << std::hex;
						else if (!VR.compare("SL"))
							if (isBigEndian)
								std::cout << std::dec << (((int)(dat[0]) << 24) +
										 ((int)(dat[1]) << 16) +
										 ((int)(dat[2]) << 8) +
										  (int)(dat[3])) << std::hex;
							else
								std::cout << std::dec << (((int)(dat[3]) << 24) +
										 ((int)(dat[2]) << 16) +
										 ((int)(dat[1]) << 8) +
										  (int)(dat[0])) << std::hex;
						else if (!VR.compare("SS"))
							if (isBigEndian)
								std::cout << std::dec << (((short int)(dat[0]) << 8) +
										 (short int)(dat[1])) << std::hex;
							else
								std::cout << std::dec << (((short int)(dat[1]) << 8) +
										 (short int)(dat[0])) << std::hex;
						else if (!VR.compare("UL"))
							if (isBigEndian)
								std::cout << std::dec << (unsigned int)(((int)(dat[0]) << 24) +
										 ((int)(dat[1]) << 16) +
										 ((int)(dat[2]) << 8) +
										  (int)(dat[3])) << std::hex;
							else
								std::cout << std::dec << (unsigned int)(((int)(dat[3]) << 24) +
										 ((int)(dat[2]) << 16) +
										 ((int)(dat[1]) << 8) +
										  (int)(dat[0])) << std::hex;
						else if (!VR.compare("US"))
							if (isBigEndian)
								std::cout << std::dec << (unsigned short int)(((short int)(dat[0]) << 8) +
										 (short int)(dat[1])) << std::hex;
							else
								std::cout << std::dec << (unsigned short int)(((short int)(dat[1]) << 8) +
										 (short int)(dat[0])) << std::hex;
						else if (!VR.compare("SQ"))
							std::cout << "Sequence printed as strings below";
						else 
							std::cout << "Unsupported format";
                    else
                        std::cout << "Data larger than " << std::dec
                                  << avoidWarning << std::hex;
                    std::cout << "\n";
				#endif
				#if defined(OUTPUT_ALL) || defined(OUTPUT_TAG)
                    std::cout << std::dec << "\n";
				#endif

                // Save proper transfer syntax for farther parsing
                if (temp->tag[0] == 0x0002 && temp->tag[1] == 0x0010) {
                    QString TransSyntax((char*)temp->vf);
                    if (!TransSyntax.compare("1.2.840.10008.1.2.1")) {
                        isImplicit = false;
                        isBigEndian = false;
                    }
                    else if (!TransSyntax.compare("1.2.840.10008.1.2.2")) {
                        isImplicit = false;
                        isBigEndian = true;
                    }
                    else if (!TransSyntax.compare("1.2.840.10008.1.2")) {
                        isImplicit = true;
                        isBigEndian = false;
                    }
                    else {
                        std::cout << "Unknown transfer syntax, assuming explicit and little endian\n";
                        isImplicit = false;
                        isBigEndian = false;
                    }
                }

                // Save slice height for later sorting
                if (temp->tag[0] == 0x0020 && temp->tag[1] == 0x1041) {
					QString tempS = "";
					for (unsigned int s = 0; s < temp->vl; s++) {
						tempS.append(temp->vf[s]);
					}

					z = tempS.toDouble();
                }
            }
            else if (nested) {
                #if defined(OUTPUT_ALL) || defined(OUTPUT_TAG)
                    std::cout << "Nested data\n";
                    for (int i = 0; i < temp->seq.items.size(); i++) {
                        std::cout << "\t" << std::dec << i+1 << ") ";
						unsigned long int avoidWarning = (unsigned long int)MAX_DATA_PRINT;
						if (avoidWarning == 0 || temp->seq.items[i]->vl < avoidWarning)
							for (unsigned int j = 0; j < temp->seq.items[i]->vl; j++)
								std::cout << std::hex << (*(temp->seq.items[i])).vf[j];
						else
							std::cout << "Data larger than " << std::dec << avoidWarning << std::hex;
							
                        std::cout << std::hex << "\n";
					}
                    std::cout << "\n";
                #endif
                nested = false;
            }
            data.append(temp);
            /*============================================================================*/
            /*REPEAT UNTIL EOF============================================================*/
        }
        file.close();
        return l;
    }
    return 0;
}

int DICOM::parseSequence(QDataStream *in, QVector <Attribute*> *att) {
	unsigned char *dat;
	in->setByteOrder(QDataStream::LittleEndian);
	Attribute *temp;
	unsigned int size;
	QString VR;
	bool nested = false;
	#if defined(OUTPUT_SQ)
		std::cout << "\nEntering the parsing loop\n"; std::cout.flush();
	#endif
	while (!in->atEnd()) {
		nested = false;
		temp = new Attribute();

		// Get the tag
		dat = new unsigned char[4];
		if (in->readRawData((char*)dat,4) != 4) {
			// Not a DICOM file
			delete[] dat;
			delete temp;
			return 1;
		}
		#if defined(OUTPUT_SQ)
		    std::cout << "Tag "; std::cout.flush();
		#endif
		
		temp->tag[0]= ((unsigned short int)(dat[1]) << 8) +
					  (unsigned short int)dat[0];
		temp->tag[1]= ((unsigned short int)(dat[3]) << 8) +
					  (unsigned short int)dat[2];
		delete[] dat;
		#if defined(OUTPUT_SQ)
		    std::cout << std::hex << temp->tag[0] << "," <<  temp->tag[1] << " | Representation " << std::dec; std::cout.flush();
		#endif
		
		// Get the VR
		dat = new unsigned char[4];		
		if (!isImplicit || temp->tag[0] == 0x0002) {
			if (in->readRawData((char*)dat,4) != 4) {
				// Not a DICOM file
				delete[] dat;
				delete temp;
				return 0;
			}

			VR = QString(dat[0])+dat[1];
			#if defined(OUTPUT_SQ)
			    std::cout << ((unsigned short int)(dat[0]) << 8) + (unsigned short int)dat[1] << " -> " << VR.toStdString() << " | Size "; std::cout.flush(); 
			#endif
		}
		else {
			VR = lib->binSearch(temp->tag[0], temp->tag[1], 0, lib->lib.size()-1).vr;
			#if defined(OUTPUT_SQ)
			    std::cout << VR.toStdString() << " (implicit) | Size ";  std::cout.flush();
			#endif
		}
		
		// Get size
		if ((temp->tag[0] != 0x0002 && isImplicit) || (lib->implicitVR.contains(VR))) {
			if (in->readRawData((char*)dat,4) != 4) { //Reread for size
				// Not a DICOM file
				delete[] dat;
				delete temp;
				return 0;
			}
			temp->vl = ((unsigned int)(dat[3]) << 24) +
					   ((unsigned int)(dat[2]) << 16) +
					   ((unsigned int)(dat[1]) << 8) +
					   (unsigned int)dat[0];

			// We have a sequence
			if (!VR.compare("SQ") && temp->vl == (unsigned int)0xFFFFFFFF) {
				nested = true;
				if (!readSequence(in, temp)) {
					delete[] dat;
					delete temp;
					return 0;
				}
			}
			else if (!VR.compare("SQ")) {
				nested = true;
				if (!readDefinedSequence(in, temp, temp->vl)) {
					delete[] dat;
					delete temp;
					return 0;
				}
			}
		}
		else {
			if (isImplicit && temp->tag[0] != 0x0002)
				if (in->readRawData((char*)dat,4) != 4) {
					// Not a DICOM file
					delete[] dat;
					delete temp;
					return 0;
				}

			if (lib->validVR.contains(VR))
				temp->vl = ((unsigned short int)(dat[3]) << 8) +
						   (unsigned short int)dat[2];
			else
				temp->vl = ((unsigned int)(dat[3]) << 24) +
						   ((unsigned int)(dat[2]) << 16) +
						   ((unsigned int)(dat[1]) << 8) +
						   (unsigned int)dat[0];
						   
			// We have a sequence
			if (!VR.compare("SQ") && temp->vl == (unsigned int)0xFFFFFFFF) {
				nested = true;
				if (!readSequence(in, temp)) {
					delete[] dat;
					delete temp;
					return 0;
				}
			}
			else if (!VR.compare("SQ")) {
				nested = true;
				if (!readDefinedSequence(in, temp, temp->vl)) {
					delete[] dat;
					delete temp;
					return 0;
				}
			}	
		}

		if (temp->vl == (unsigned int)0xFFFFFFFF) {
			temp->vl = 0;
		}

		size = temp->vl;
		#if defined(OUTPUT_SQ)
		    std::cout << temp->vl << " -> " << std::dec << size << "\n";
		#endif
		
		Reference closest = lib->binSearch(temp->tag[0], temp->tag[1], 0, lib->lib.size()-1);
		if (closest.tag[0] == temp->tag[0] && closest.tag[1] == temp->tag[1])
			temp->desc = closest.title;
		else
			temp->desc = "Unknown Tag";

		delete[] dat;

		// Get data
		if (!nested) {
			temp->vf = new unsigned char[size];
			if (size > 0 && size < (unsigned long int)INT_MAX) {
				if (in->readRawData((char*)temp->vf,size) != (long int)size) {
					// Not a DICOM file
					delete temp;
					return 0;
				}
			}
			else if (size > 0) { // In case we need to read in more data then the buffer can handle
				unsigned char *pt = temp->vf;
				for (int i = 0; (unsigned long int)i < size/((unsigned long int)INT_MAX); i++) {
					if (in->readRawData((char*)pt,INT_MAX) != INT_MAX) {
						// Not a DICOM file
						delete temp;
						return 0;
					}
					pt += sizeof(char)*INT_MAX;
				}
			}
		}
		
		#if defined(OUTPUT_SQ)
		    if (!nested) {
		    	QString tempS = "";
		    	for (unsigned int s = 0; s < temp->vl; s++)
		    		tempS.append(temp->vf[s]);
		    	std::cout << tempS.toStdString() << "\n";
		    }
		    else {
		    	for (int i = 0; i < temp->seq.items.size(); i++) {
		    		std::cout << i << ") ";
		    		for (unsigned int j = 0; j < temp->seq.items[i]->vl; j++)
		    			std::cout << (*(temp->seq.items[i])).vf[j];
		    		std::cout << "\n";
		    	}
		    }
		#endif
		
		att->append(temp);
	}
	return att->size();
}