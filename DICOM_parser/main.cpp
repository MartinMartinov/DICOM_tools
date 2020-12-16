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
	
	if (argc == 1) {
        std::cout << "Please call this program with one or more .dcm files.\n";
        return 0;
    }

    database dat;
    QVector <DICOM *> dicom;

    for (int i = 0; i < argc-1; i++) {
        QString path(argv[i+1]);
        DICOM *d = new DICOM(&dat);
        if (!d->parse(path)) {
            std::cout << "Unsuccessfully parsed " << path.toStdString() << ", quitting...\n";
            for (int j = 0; j < dicom.size(); j++) {
                delete dicom[j];
            }
            dicom.clear();
            delete d;
            return -1;
        }
		dicom.append(d);
    }
	
	duration = (std::clock()-start)/(double)CLOCKS_PER_SEC;
	std::cout << "Parsed the " << dicom.size() << " DICOM files.  Time elapsed is " << duration << " s.\n";
	
    for (int j = 0; j < dicom.size(); j++) {
        delete dicom[j];
    }
    dicom.clear();
    return 1;
}