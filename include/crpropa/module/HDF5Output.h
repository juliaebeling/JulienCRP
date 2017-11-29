#ifdef CRPROPA_HAVE_HDF5

#ifndef CRPROPA_HDF5OUTPUT_H
#define CRPROPA_HDF5OUTPUT_H

#include "crpropa/module/Output.h"
#include "stdint.h"
#include <ctime>

#include <H5Ipublic.h>

namespace crpropa {

const size_t propertyBufferSize = 1024;

class HDF5Output: public Output {

	typedef struct OutputRow {
		double D;
		double z;
		uint64_t SN;
		int32_t ID;
		double E;
		double X;
		double Y;
		double Z;
		double Px;
		double Py;
		double Pz;
		uint64_t SN0;
		int32_t ID0;
		double E0;
		double X0;
		double Y0;
		double Z0;
		double P0x;
		double P0y;
		double P0z;
		uint64_t SN1;
		int32_t ID1;
		double E1;
		double X1;
		double Y1;
		double Z1;
		double P1x;
		double P1y;
		double P1z;
		double weight;
		unsigned char propertyBuffer[propertyBufferSize];
	} OutputRow;

	std::string filename;

	hid_t file, sid;
	hid_t dset, dataspace;
	mutable std::vector<OutputRow> buffer;

	time_t lastFlush;
public:
	HDF5Output(const std::string &filename);
	HDF5Output(const std::string &filename, OutputType outputtype);
	~HDF5Output();

	void process(Candidate *candidate) const;
	std::string getDescription() const;

	void open(const std::string &filename);
	void close();
	void flush() const;

};

} // namespace crpropa

#endif // CRPROPA_HDF5OUTPUT_H

#endif // CRPROPA_HAVE_HDF5
