/***********************************************************************
************************************************************************
*    This code is part of a PRE-RELEASE VERSION of 3ddose_tools,       *
*    a *.3ddose file analysis code.                                    *
*    Copyright 2014 by Carleton University, Ottawa, Canada             *
*    GNU GENERAL PUBLIC LICENSE                                        *
*    Version 3, 29 June 2007                                           *
*                                                                      *
*    Please report all problems to:                                    *
*    Martin Martinov martinov@physics.carleton.ca                      *
*    Rowan Thomson rthomson@physics.carleton.ca                        *
************************************************************************
***********************************************************************/

#include "egsphant.h"

EGSPhant::EGSPhant() {
    nx = ny = nz = 0;
}

// Make a mask template from another EGSPhant
void EGSPhant::makeMask(EGSPhant* mask) {
    nx = mask->nx;
	ny = mask->ny;
	nz = mask->nz;
    x = mask->x;
	y = mask->y;
	z = mask->z;
    maxDensity = mask->maxDensity;
	{
		QVector <char> mz(nz, 49);
		QVector <QVector <char> > my(ny, mz);
		QVector <QVector <QVector <char> > > mx(nx, my);
		m = mx;
		QVector <double> dz(nz, 0);
		QVector <QVector <double> > dy(ny, dz);
		QVector <QVector <QVector <double> > > dx(nx, dy);
		d = dx;
	}
    media << "OTHER" << "TARGET";
}

void EGSPhant::loadEGSPhantFile(QString path) {
    QFile file(path);

    // Increment size of the status bar
    double increment;

    // Open up the file specified at path
    if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream input(&file);
        QString line = input.readLine();

        // read in the number of media
        int num = line.trimmed().toInt();
        media.resize(num);

        // read the media into an array
        for (int i = 0; i < num; i++) {
            media[i] = input.readLine().trimmed();
        }

        // skim over the the ESTEP info
        line = input.readLine().trimmed();

        // read in the dimensions of the egsphant file and
        // store the size and resize the matrices holding the boundaries
        input.skipWhiteSpace();
        input >> nx;
        x.fill(0,nx+1);
        input.skipWhiteSpace();
        input >> ny;
        y.fill(0,ny+1);
        input.skipWhiteSpace();
        input >> nz;
        z.fill(0,nz+1);

        // resize the 3D matrix to hold all densities
        {
            QVector <char> mz(nz, 0);
            QVector <QVector <char> > my(ny, mz);
            QVector <QVector <QVector <char> > > mx(nx, my);
            m = mx;
            QVector <double> dz(nz, 0);
            QVector <QVector <double> > dy(ny, dz);
            QVector <QVector <QVector <double> > > dx(nx, dy);
            d = dx;
        }

        // read in all the boundaries of the phantom
        input.skipWhiteSpace();
        for (int i = 0; i <= nx; i++) {
            input >> x[i];
        }
        for (int i = 0; i <= ny; i++) {
            input >> y[i];
        }
        for (int i = 0; i <= nz; i++) {
            input >> z[i];
        }

        // Determine the increment this egsphant file gets
        increment = MAX_PROGRESS/double(nz-1);

        // Read in all the media
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++) {
                    input.skipWhiteSpace();
                    input >> m[i][j][k];
                }
            emit progressMade(increment); // Update progress bar
        }

        file.close();
    }
}

void EGSPhant::loadEGSPhantFilePlus(QString path) {
    QFile file(path);

    // Increment size of the status bar
    double increment;

    // Open up the file specified at path
    if (file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QTextStream input(&file);
        QString line = input.readLine();

        // read in the number of media
        int num = line.trimmed().toInt();
        media.resize(num);

        // read the media into an array
        for (int i = 0; i < num; i++) {
            media[i] = input.readLine().trimmed();
        }

        // skim over the the ESTEP info
        line = input.readLine().trimmed();

        // read in the dimensions of the egsphant file
        input.skipWhiteSpace();
        input >> nx;
        x.fill(0,nx+1);
        input.skipWhiteSpace();
        input >> ny;
        y.fill(0,ny+1);
        input.skipWhiteSpace();
        input >> nz;
        z.fill(0,nz+1);

        // resize the 3D matrix to hold all densities
        {
            QVector <char> mz(nz, 0);
            QVector <QVector <char> > my(ny, mz);
            QVector <QVector <QVector <char> > > mx(nx, my);
            m = mx;
            QVector <double> dz(nz, 0);
            QVector <QVector <double> > dy(ny, dz);
            QVector <QVector <QVector <double> > > dx(nx, dy);
            d = dx;
        }

        // read in all the boundaries of the phantom
        input.skipWhiteSpace();
        for (int i = 0; i <= nx; i++) {
            input >> x[i];
        }
        for (int i = 0; i <= ny; i++) {
            input >> y[i];
        }
        for (int i = 0; i <= nz; i++) {
            input >> z[i];
        }

        // Determine the increment this egsphant file gets
        increment = MAX_PROGRESS/double(nz-1);

        // Read in all the media
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++) {
                    input.skipWhiteSpace();
                    input >> m[i][j][k];
                }
            emit progressMade(increment/100.0*10.0); // Update progress bar
        }

        // Read in all the densities
        maxDensity = 0;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++) {
                    input.skipWhiteSpace();
                    input >> d[i][j][k];
                    if (d[i][j][k] > maxDensity) {
                        maxDensity = d[i][j][k];
                    }

                }
            emit progressMade(increment/100.0*90.0); // Update progress bar
        }

        file.close();
    }
}

void EGSPhant::saveEGSPhantFile(QString path) {
    QFile file(path);

    // Increment size of the status bar
    double increment;

    // Open up the file specified at path
    if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream output(&file);
		
        // read out the number of media
        output << media.size() << "\n";

        // read out the media into an array
        for (int i = 0; i < media.size(); i++) {
            output << media[i] << "\n";
        }

        // print generic ESTEP info
        for (int i = 0; i < media.size(); i++) {
            output << "0.50 ";
        }
		output << "\n";

        // read out the dimensions of the egsphant file
        output << nx << " ";
        output << ny << " ";
        output << nz << "\n";

        // read out all the boundaries of the phantom
        for (int i = 0; i <= nx; i++) {
            output << x[i] << " ";
        }
		output << "\n";
        for (int i = 0; i <= ny; i++) {
            output << y[i] << " ";
        }
		output << "\n";
        for (int i = 0; i <= nz; i++) {
            output << z[i] << " ";
        }
		output << "\n";

        // Determine the increment this egsphant file gets
        increment = MAX_PROGRESS/double(nz-1);

        // Read out all the media
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    output << m[i][j][k];
                }
				output << "\n";
			}
            emit progressMade(increment/100.0*10.0); // Update progress bar
            output << "\n";
        }

        // Read out all the densities
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    output << d[i][j][k] << " ";
                }
				output << "\n";
            }
            emit progressMade(increment/100.0*90.0); // Update progress bar
            output << "\n";
        }

        file.close();
    }
}

void EGSPhant::loadbEGSPhantFile(QString path) {
    QFile file(path);

    // Increment size of the status bar
    double increment;

    // Open up the file specified at path
    if (file.open(QIODevice::ReadOnly)) {
        QDataStream input(&file);
		input.setByteOrder(QDataStream::LittleEndian);
		
        // read in the number of media
        unsigned char num;
        input >> num;
        media.resize(num);

        // read the media into an array
        char *temp2;
        for (int i = 0; i < num; i++) {
            input >> temp2;
            media[i] = QString(temp2);
        }

        // skim over the the ESTEP info
        double temp;
        for (int i = 0; i < num; i++) {
            input >> temp;
        }

        // read in the dimensions of the egsphant file
        // store the size and resize the matrices holding the boundaries
        input >> nx;
        x.fill(0,nx+1);
        input >> ny;
        y.fill(0,ny+1);
        input >> nz;
        z.fill(0,nz+1);

        // resize the 3D matrix to hold all densities
        {
            QVector <char> mz(nz, 0);
            QVector <QVector <char> > my(ny, mz);
            QVector <QVector <QVector <char> > > mx(nx, my);
            m = mx;
            QVector <double> dz(nz, 0);
            QVector <QVector <double> > dy(ny, dz);
            QVector <QVector <QVector <double> > > dx(nx, dy);
            d = dx;
        }

        // read in all the boundaries of the phantom
        for (int i = 0; i <= nx; i++) {
            input >> x[i];
        }
        for (int i = 0; i <= ny; i++) {
            input >> y[i];
        }
        for (int i = 0; i <= nz; i++) {
            input >> z[i];
        }

        // Determine the increment this egsphant file gets
        increment = MAX_PROGRESS/double(nz-1);

        // Read in all the media
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++) {
                    input >> num;
                    m[i][j][k] = num;
                }
            emit progressMade(increment); // Update progress bar
        }

        file.close();
    }
}

void EGSPhant::loadbEGSPhantFilePlus(QString path) {
    QFile file(path);

    // Increment size of the status bar
    double increment;

    // Open up the file specified at path
    if (file.open(QIODevice::ReadOnly)) {
        QDataStream input(&file);
		input.setByteOrder(QDataStream::LittleEndian);
		
        // read in the number of media
        unsigned char num;
        input >> num;
        media.resize(num);

        // read the media into an array
        char *temp2;
        for (int i = 0; i < num; i++) {
            input >> temp2;
            media[i] = QString(temp2);
        }

        // skim over the the ESTEP info
        double temp;
        for (int i = 0; i < num; i++) {
            input >> temp;
        }

        // read in the dimensions of the egsphant file
        // store the size and resize the matrices holding the boundaries
        input >> nx;
        x.fill(0,nx+1);
        input >> ny;
        y.fill(0,ny+1);
        input >> nz;
        z.fill(0,nz+1);

        // resize the 3D matrix to hold all densities
        {
            QVector <char> mz(nz, 0);
            QVector <QVector <char> > my(ny, mz);
            QVector <QVector <QVector <char> > > mx(nx, my);
            m = mx;
            QVector <double> dz(nz, 0);
            QVector <QVector <double> > dy(ny, dz);
            QVector <QVector <QVector <double> > > dx(nx, dy);
            d = dx;
        }

        // read in all the boundaries of the phantom
        for (int i = 0; i <= nx; i++) {
            input >> x[i];
        }
        for (int i = 0; i <= ny; i++) {
            input >> y[i];
        }
        for (int i = 0; i <= nz; i++) {
            input >> z[i];
        }

        // Determine the increment this egsphant file gets
        increment = MAX_PROGRESS/double(nz-1);

        // Read in all the media
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++) {
                    input >> num;
                    m[i][j][k] = num;
                }
            emit progressMade(increment/100.0*50.0); // Update progress bar
        }

        // Read in all the densities
        maxDensity = 0;
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++) {
                    input >> d[i][j][k];
                    if (d[i][j][k] > maxDensity) {
                        maxDensity = d[i][j][k];
                    }
                }
            emit progressMade(increment/100.0*50.0); // Update progress bar
        }

        file.close();
    }
}

void EGSPhant::savebEGSPhantFile(QString path) {
    QFile file(path);

    // Increment size of the status bar
    double increment;

    // Open up the file specified at path
    if (file.open(QIODevice::WriteOnly)) {
        QDataStream output(&file);
		output.setByteOrder(QDataStream::LittleEndian);
		
        // read out the number of media
        output << media.size();

        // read out the media into an array
        for (int i = 0; i < media.size(); i++) {
            output << media[i];
        }

        // print generic ESTEP info
        double temp = 0.5;
        for (int i = 0; i < media.size(); i++) {
            output << temp;
        }

        // read out the dimensions of the egsphant file
        // store the size and resize the matrices holding the boundaries
        output << nx;
        output << ny;
        output << nz;

        // read out all the boundaries of the phantom
        for (int i = 0; i <= nx; i++)
            output << x[i];
        for (int i = 0; i <= ny; i++)
            output << y[i];
        for (int i = 0; i <= nz; i++)
            output << z[i];

        // Determine the increment this egsphant file gets
        increment = MAX_PROGRESS/double(nz-1);

        // Read out all the media
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                    output << m[i][j][k];
            emit progressMade(increment/100.0*50.0); // Update progress bar
        }

        // Read out all the densities
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                    output << d[i][j][k];
            emit progressMade(increment/100.0*50.0); // Update progress bar
        }

        file.close();
    }
}

char EGSPhant::getMedia(double px, double py, double pz) {
    int ix, iy, iz;
    ix = iy = iz = -1;

    // Find the index of the boundary that is less than px, py and pz
    for (int i = 0; i < nx; i++)
        if (px <= x[i+1]) {
            ix = i;
            break;
        }
    if (px < x[0]) {
        ix = -1;
    }

    for (int i = 0; i < ny; i++)
        if (py <= y[i+1]) {
            iy = i;
            break;
        }
    if (py < y[0]) {
        iy = -1;
    }

    for (int i = 0; i < nz; i++)
        if (pz <= z[i+1]) {
            iz = i;
            break;
        }
    if (pz < z[0]) {
        iz = -1;
    }

    // This is to insure that no area outside the vectors is accessed
    if (ix < nx && ix >= 0 && iy < ny && iy >= 0 && iz < nz && iz >= 0) {
        return m[ix][iy][iz];
    }

    return 0; // We are not within our bounds
}

double EGSPhant::getDensity(double px, double py, double pz) {
    int ix, iy, iz;
    ix = iy = iz = -1;

    // Find the index of the boundary that is less than px, py and pz
    for (int i = 0; i < nx; i++)
        if (px <= x[i+1]) {
            ix = i;
            break;
        }
    if (px < x[0]) {
        ix = -1;
    }

    for (int i = 0; i < ny; i++)
        if (py <= y[i+1]) {
            iy = i;
            break;
        }
    if (py < y[0]) {
        iy = -1;
    }

    for (int i = 0; i < nz; i++)
        if (pz <= z[i+1]) {
            iz = i;
            break;
        }
    if (pz < z[0]) {
        iz = -1;
    }

    // This is to insure that no area outside the vectors is accessed
    if (ix < nx && ix >= 0 && iy < ny && iy >= 0 && iz < nz && iz >= 0) {
        return d[ix][iy][iz];
    }

    return 0; // We are not within our bounds
}

double EGSPhant::interpDen(double xp, double yp, double zp) {
	if (xp <= x[1] || x[nx-1] <= xp || yp <= y[1] || y[ny-1] <= yp || zp <= z[1] || z[nz-1] <= zp)
		return 0;
	
    // All the positions needed (in addition to the ones passed in)
    double x0, y0, z0, x1, y1, z1;

    // Indices
    int xi, yi, zi;

    // Define X
    xi = getIndex("x axis", xp);
    if (xp < (x[xi] + x[xi-1])/2.0) {
        x0 = (x[xi] + x[xi-1])/2.0;
        x1 = (x[xi] + x[xi+1])/2.0;
    }
    else {
        x0 = (x[xi] + x[xi+1])/2.0;
        x1 = (x[xi+1] + x[xi+2])/2.0;
    }

    // Define Y
    yi = getIndex("y axis", yp);
    if (yp < (y[yi] + y[yi-1])/2.0) {
        y0 = (y[yi] + y[yi-1])/2.0;
        y1 = (y[yi] + y[yi+1])/2.0;
    }
    else {
        y0 = (y[yi] + y[yi+1])/2.0;
        y1 = (y[yi+1] + y[yi+2])/2.0;
    }

    // Define Z
    zi = getIndex("z axis", zp);
    if (zp < (z[zi] + z[zi-1])/2.0) {
        z0 = (z[zi] + z[zi-1])/2.0;
        z1 = (z[zi] + z[zi+1])/2.0;
    }
    else {
        z0 = (z[zi] + z[zi+1])/2.0;
        z1 = (z[zi+1] + z[zi+2])/2.0;
    }

    // Precompute lengths
    double ix = xp-x0, fx = x1-xp, iy = yp-y0, fy = y1-yp, iz = zp-z0, fz = z1-zp;
    double val = 0, vol = 0, tVol = 0;

    // Add all 8 values together, weighted by the volume of the rectangular
    // prism formed by the point that was passed in and the one in the opposite
    // corner of the total rectangular prism volume
    vol = fx*fy*fz;
    val += getDensity(x0,y0,z0)*vol;
	tVol += vol;

    vol = ix*fy*fz;
    val += getDensity(x1,y0,z0)*vol;
	tVol += vol;

    vol = fx*iy*fz;
    val += getDensity(x0,y1,z0)*vol;
	tVol += vol;

    vol = fx*fy*iz;
    val += getDensity(x0,y0,z1)*vol;
	tVol += vol;

    vol = fx*iy*iz;
    val += getDensity(x0,y1,z1)*vol;
	tVol += vol;

    vol = ix*fy*iz;
    val += getDensity(x1,y0,z1)*vol;
	tVol += vol;

    vol = ix*iy*fz;
    val += getDensity(x1,y1,z0)*vol;
	tVol += vol;

    vol = ix*iy*iz;
    val += getDensity(x1,y1,z1)*vol;
	tVol += vol;

    val /= tVol;
	//if  (val > 10000)
	//	std::cout << ix << " - " << fx << ", " << iy << " - " << fy << ", " << iz << " - " << fz << "\n";

    return val;
}

int EGSPhant::getIndex(QString axis, double p) {
    int index = -1;

    if (!axis.compare("x axis")) {
        for (int i = 0; i < nx; i++)
            if (p <= x[i+1]) {
                index = i;
                break;
            }
        if (p < x[0]) {
            index = -1;
        }
    }
    else if (!axis.compare("y axis")) {
        for (int i = 0; i < ny; i++)
            if (p <= y[i+1]) {
                index = i;
                break;
            }
        if (p < y[0]) {
            index = -1;
        }
    }
    else if (!axis.compare("z axis")) {
        for (int i = 0; i < nz; i++)
            if (p <= z[i+1]) {
                index = i;
                break;
            }
        if (p < z[0]) {
            index = -1;
        }
    }

    return index; // -1 if we are out of bounds
}

QImage EGSPhant::getEGSPhantPicMed(QString axis, double ai, double af,
                                   double bi, double bf, double d, double res) {
    // Create a temporary image
    int width  = (af-ai)*res;
    int height = (bf-bi)*res;
    QImage image(height, width, QImage::Format_RGB32);
    double hInc, wInc, cInc;
    double h, w, c = 0;
	
	// Get the axis
	int ax;
	if (!axis.compare("x axis"))
		ax = 1;
	else if (!axis.compare("y axis"))
		ax = 2;
	else if (!axis.compare("z axis"))
		ax = 3;
	else
		return image;

    // Calculate the size (in cm) of pixels, and then the range for grayscaling
    wInc = 1/res;
    hInc = 1/res;
    cInc = 255.0/(media.size()-1);

    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++) {
            // determine the location of the current pixel in the phantom
            h = (double(bi)) + hInc * double(i);
            w = (double(ai)) + wInc * double(j);

            // get the media, which differs based on axis through which image is
            // sliced
            if (ax == 1) {
                c = getMedia(d, h, w) - 49;
				c -= (c>9?17:0);
            }
            else if (ax == 2) {
                c = getMedia(h, d, w) - 49;
				c -= (c>9?17:0);
            }
            else if (ax == 3) {
                c = getMedia(h, w, d) - 49;
				c -= (c>9?17:0);
            }
            // finally, paint the pixel
            image.setPixel(i, width-1-j, qRgb(int(cInc*c), int(cInc*c), int(cInc*c)));
        }

    return image; // return the image created
}

QImage EGSPhant::getEGSPhantPicDen(QString axis, double ai, double af,
                                   double bi, double bf, double d, double res) {
    // Create a temporary image
    int width  = (af-ai)*res;
    int height = (bf-bi)*res;
    QImage image(height, width, QImage::Format_RGB32);
    double hInc, wInc, cInc;
    double h, w, c = 0;

	// Get the axis
	int ax;
	if (!axis.compare("x axis"))
		ax = 1;
	else if (!axis.compare("y axis"))
		ax = 2;
	else if (!axis.compare("z axis"))
		ax = 3;
	else
		return image;
	
    // Calculate the size (in cm) of pixels, and then the range for grayscaling
    wInc = 1/res;
    hInc = 1/res;
    cInc = 255.0/maxDensity;

    for (int i = 0; i < height; i++)
        for (int j = 0; j < width; j++) {
            // determine the location of the current pixel in the phantom
            h = (double(bi)) + hInc * double(i);
            w = (double(ai)) + wInc * double(j);

            // get the density, which differs based on axis through which image
            // os sliced
            if (ax == 1) {
                c = getDensity(d, h, w);
            }
            else if (ax == 2) {
                c = getDensity(h, d, w);
            }
            else if (ax == 3) {
                c = getDensity(h, w, d);
            }

            // finally, paint the pixel
            image.setPixel(i, width-1-j, qRgb(int(cInc*c), int(cInc*c), int(cInc*c)));
        }

    return image; // return the image created
}
