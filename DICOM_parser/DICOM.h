#ifndef DICOM_H
#define DICOM_H

#include <QtGui>
#include <iostream>
#include <math.h>

// These need to be declared ahead of time, they are needed for nested sequences
class Sequence;
class SequenceItem;
class Attribute;

// The following two classes are used to hold a sequence of items (and yes, you
// can have nested sequences, cause, you know, why not?)
class Sequence {
public:
    QVector <SequenceItem *> items;
    ~Sequence();
};

class SequenceItem {
public:
    unsigned long int vl; // Value Length
    unsigned char *vf; // Value Field
    Sequence seq; // Contains potential sequences

    SequenceItem(unsigned long int size, unsigned char *data);
    SequenceItem(unsigned long int size, Attribute *data);
    ~SequenceItem();
};

// Might as well be a struct, but I might want some methods in the future
class Attribute {
public:
    unsigned short int tag[2]; // Element Identifier
    QString desc; // Desciption
    unsigned short int vr; // Value Representation
    unsigned long int vl; // Value Length
    unsigned char *vf; // Value Field
    Sequence seq; // Contains potential sequences

    Attribute();
    ~Attribute();
};

// These are all defined in database.cpp so as to save alot of recompiling
// hassle
struct Reference {
    unsigned short int tag[2]; // Element Identifier
    QString vr; // Element Identifier
    QString title; // Title of element
};

class database : public QObject {
    Q_OBJECT

public:
    // Contains a list of known attribute entries
    QVector <Reference *> lib;
    // Contains all the accepteable value representations
    QStringList validVR;
    QStringList implicitVR;

    database();
    ~database();

    Reference binSearch(unsigned short int one, unsigned short int two, int min,
                        int max);
};

class DICOM : public QObject {
    Q_OBJECT

public:
    // Contains all the data read in from a dicom file sorted into attributes
    QVector <Attribute *> data;
	
    // Pointer to precompiled DICOM library
    database *lib;
	
	// Transfer syntax
    bool isImplicit, isBigEndian;
	
	// z height (default to NaN, only change if slice height tag is found)
	double z = std::nan("1");

	// file location for later lookup
	QString path;

    DICOM(database *);
    ~DICOM();

    int parse(QString p);
    int readSequence(QDataStream *in, Attribute *att);
    int readDefinedSequence(QDataStream *in, Attribute *att, unsigned long int n = 0);
	
	int parseSequence(QDataStream *in, QVector <Attribute*> *att);
};

#endif
