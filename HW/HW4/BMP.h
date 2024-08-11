// include header files
#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <algorithm> 
#include <vector>
#include <climits>
#include "fftw/fftw3.h"

using namespace std;

class BMP
{
	public:
		
		// variable
		char* 	_header;
		uint16_t	_identifier;			// 0x0000
		uint32_t	_file_size;				// 0x0002
		uint32_t	_reserved;				// 0x0006
		uint32_t	_bitmap_data_offset;	// 0x000A
		uint32_t	_header_size;			// 0x000E
		uint32_t	_width;					// 0x0012
	    uint32_t	_height;				// 0x0016
	    uint16_t	_planes;				// 0x001A
	    uint16_t	_bits_per_pixel;		// 0x001C
	    uint32_t	_compression;			// 0x001E
	    uint32_t	_bitmap_data_size;		// 0x0022
	    uint32_t	_h_resolution;			// 0x0026
	    uint32_t	_v_resolution;			// 0x002A
	    uint32_t	_used_colors;			// 0x002E
	    uint32_t	_important_colors;		// 0x0032
	    uint32_t	_channel;
	    
	    char* 	    _image;
	    
	    // function
		bool    		read(const string filename);
		void    		show();
		bool    		write(const string filename);
		double** 		gaussianFilter(int kernelSize);
		double**		motionFilter(int kernelSize);
		double*			pad(int kernelSize, double** f);
		void			wiener(fftw_complex* f, fftw_complex* finv, double k);
		void			fMultiply(fftw_complex* fB, fftw_complex* fG, fftw_complex* fR, fftw_complex* f);
		void			restore(int gSize, int mSize, double k);
};

bool BMP::read(const string filename)
{
	// file open
	fstream in_bmp(filename, ios::in | ios::binary);
	if(!in_bmp.is_open())
		cerr << "BMP Open Fail\n";
	
	// read header
	_header = new char[54];
	in_bmp.read(_header, sizeof(char) * 54);		
	_identifier 		= *(uint16_t*)  &_header[0];
	_file_size 			= *(uint32_t*)  &_header[2];
	_reserved 			= *(uint32_t*)  &_header[6];
	_bitmap_data_offset = *(uint32_t*)  &_header[10];
	_header_size 		= *(uint32_t*)  &_header[14];
	_width 				= *(uint32_t*)  &_header[18];
    _height 			= *(uint32_t*)  &_header[22];
    _planes 			= *(uint32_t*)  &_header[26];
    _bits_per_pixel 	= *(uint32_t*)  &_header[28];
    _compression 		= *(uint32_t*)  &_header[30];
    _bitmap_data_size 	= *(uint32_t*)  &_header[34];
    _h_resolution 		= *(uint32_t*)  &_header[38];
    _v_resolution 		= *(uint32_t*)  &_header[42];
    _used_colors 		= *(uint32_t*)  &_header[46];
    _important_colors 	= *(uint32_t*)  &_header[50];
    
    // calculate the number of channels
    _channel			= _bits_per_pixel / 8;
    
    // read image
    _image = new char[_height * _width * _channel];
	in_bmp.read(_image, sizeof(char) * _height * _width * _channel);
	
	// file close
	in_bmp.close();
	return true;
}

void BMP::show()
{
	// show the information of the input BMP
	cout << "ID: " 				 << _identifier 		<< '\n';
	cout << "File Size: " 		 << _file_size 			<< '\n';
	cout << "Reserved: " 		 << _reserved 			<< '\n';
	cout << "BMP Data Offset: "  << _bitmap_data_offset << '\n';
	cout << "Header Size: " 	 << _header_size 		<< '\n';
	cout << "Width: " 			 << _width 				<< '\n';
	cout << "Height: " 			 << _height 			<< '\n';
	cout << "Planes: " 			 << _planes 			<< '\n';
	cout << "BPP: " 			 << _bits_per_pixel 	<< '\n';
	cout << "Compression: " 	 << _compression 		<< '\n';
	cout << "BMP Data Size: " 	 << _bitmap_data_size 	<< '\n';
	cout << "H Resolution: " 	 << _h_resolution 		<< '\n';
	cout << "V Resolution: " 	 << _v_resolution 		<< '\n';
	cout << "Used Colors: " 	 << _used_colors 		<< '\n';
	cout << "Important Colors: " << _important_colors 	<< '\n';
	cout << "---------------------------------\n";
}

bool BMP::write(const string filename)
{
	// file open
	fstream out_bmp(filename, ios::out | ios::binary);
	if(!out_bmp)
		cerr << "BMP Write Fail\n";
	
	// file write
	out_bmp.write(_header, sizeof(char) * 54);
	out_bmp.write(_image , sizeof(char) * _height * _width * _channel);

	// file close
	out_bmp.close();
	return true;
}

// compute gaussian filter
double** BMP::gaussianFilter(int kernelSize)
{
	// compute radius
	int radius = kernelSize / 2;

	// declare dynamic array to store gaussian filter
	double** g = new double*[kernelSize];
	for(int i=0;i<kernelSize;i++)
		g[i] = new double[kernelSize];

	// declare sum to normalize later
	double sum = 0;

	// compute value of every point in gaussian filter
	// g(x,y) = e^(-(x^2 + y^2))
	for(int i=-radius;i<=radius;i++)
	{
		for(int j=-radius;j<=radius;j++)
		{
			g[radius-i][radius-j] = exp(-1 * (pow(i, 2) + pow(j, 2))); 
			sum += g[radius-i][radius-j];
		}
	}

	// normalize
	for(int i=0;i<kernelSize;i++)
		for(int j=0;j<kernelSize;j++)
			g[i][j] /= sum;

	return g;
}

// transfer filter size from kernelSize * kernelSize to height * width of the image
double* BMP::pad(int kernelSize, double** f)
{
	// compute radius
	int radius = kernelSize / 2;

	// declare dynamic array to store the result after pad and shift
	double* p = new double[_height * _width];

	// initialize
	for(int i=0;i<_height * _width;i++)
	{
		p[i] = 0;
	}

	// pad and shift
	for(int i=-radius;i<=radius;i++)
	{
		for(int j=-radius;j<=radius;j++)
		{
			if(i >=0 && j >= 0)
				p[i*_width+j] = f[radius-i][radius-j];
			else if(i >= 0 && j < 0)
				p[i*_width+(_width+j)] = f[radius-i][radius-j];
			else if(i < 0 && j >= 0)
				p[(_height+i)*_width+j] = f[radius-i][radius-j];
			else				
				p[(_height+i)*_width+(_width+j)] = f[radius-i][radius-j];
		}
	}

	return p;
}

// compute motion filter
double** BMP::motionFilter(int kernelSize)
{
	// declare dynamic array to store motion filter
	double** m = new double*[kernelSize];
	for(int i=0;i<kernelSize;i++)
		m[i] = new double[kernelSize];

	// get motion filter
	for(int i=0;i<kernelSize;i++)
	{
		for(int j=0;j<kernelSize;j++)
		{
            if(i - j == 0)
				// normalize
                m[i][j] = 1.0 / kernelSize;
            else 
				m[i][j] = 0;
        }
	}
	return m;
}

// compute wiener filter
// F(u,v) = [(1/H(u,v)) * (|H(u,v)|^2 / (H(u,v)|^2 + k)]
void BMP::wiener(fftw_complex* f, fftw_complex* finv, double k)
{
	double real, imag;
	for(int i=0;i<_height * _width;i++)
	{
		real = f[i][0];
		imag = f[i][1];
		finv[i][0] =  real / (real * real + imag * imag + k);
		finv[i][1] = -imag / (real * real + imag * imag + k);
	}
}

// compute multiply in frequency domain which is convolution in spatial domain
void BMP::fMultiply(fftw_complex* fB, fftw_complex* fG, fftw_complex* fR, fftw_complex* f)
{
	double real, imag, x, y;
	
	for(int i=0;i<_height * _width;i++)
	{
		x = f[i][0];
		y = f[i][1];

		real = fB[i][0];
		imag = fB[i][1];
		fB[i][0] = real * x - imag * y;
		fB[i][1] = real * y + imag * x;

		real = fG[i][0];
		imag = fG[i][1];
		fG[i][0] = real * x - imag * y;
		fG[i][1] = real * y + imag * x;

		real = fR[i][0];
		imag = fR[i][1];
		fR[i][0] = real * x - imag * y;
		fR[i][1] = real * y + imag * x;
	}
}


void BMP::restore(int gSize, int mSize, double k)
{
	// declare gaussian filter size = gSize
	int kernelSize = gSize;

	// gaussian filter
	double** gaussianF = new double*[kernelSize];
	for(int i=0;i<kernelSize;i++)
		gaussianF[i] = new double[kernelSize];
	gaussianF = gaussianFilter(kernelSize);

	// pad and shift
	double* gaussianPad = new double[_height * _width];
	gaussianPad = pad(kernelSize, gaussianF);

	// declare motion filter size = mSize
	kernelSize = mSize;
	
	// motion filter
	double** motionF = new double*[kernelSize];
	for(int i=0;i<kernelSize;i++)
		motionF[i] = new double[kernelSize];
	motionF = motionFilter(kernelSize);

	// pad and shift
	double* motionPad = new double[_height * _width];
	motionPad = pad(kernelSize, motionF);

	double* B = new double[_height * _width];
	double* G = new double[_height * _width];
	double* R = new double[_height * _width];
	for(int i=0;i<_height * _width;i++)
	{
		B[i] = *(uint8_t*) &_image[i*_channel  ];
		G[i] = *(uint8_t*) &_image[i*_channel+1];
		R[i] = *(uint8_t*) &_image[i*_channel+2];
	}

	fftw_complex* fB = new fftw_complex[_height * _width];
	fftw_complex* fG = new fftw_complex[_height * _width];
	fftw_complex* fR = new fftw_complex[_height * _width];
	fftw_complex* fmf = new fftw_complex[_height * _width];
	fftw_complex* fgf = new fftw_complex[_height * _width];

	fftw_plan dft;
	dft = fftw_plan_dft_r2c_2d(_height, _width, gaussianPad, fgf, FFTW_ESTIMATE);

	fftw_execute_dft_r2c(dft, B, fB);
	fftw_execute_dft_r2c(dft, G, fG);
	fftw_execute_dft_r2c(dft, R, fR);
	fftw_execute_dft_r2c(dft, gaussianPad, fgf);
	fftw_execute_dft_r2c(dft, motionPad, fmf);

	fftw_complex* fmfInv = new fftw_complex[_height * _width];
	fftw_complex* fgfInv = new fftw_complex[_height * _width];

	wiener(fmf, fmfInv, k);
	fMultiply(fB, fG, fR, fmfInv);
	wiener(fgf, fgfInv, k);
	fMultiply(fB, fG, fR, fgfInv);
	
	fftw_plan idft;
	idft = fftw_plan_dft_c2r_2d(_height, _width, fB, B, FFTW_ESTIMATE);

	fftw_execute_dft_c2r(idft, fB, B);
	fftw_execute_dft_c2r(idft, fG, G);
	fftw_execute_dft_c2r(idft, fR, R);

	for(int i=0;i<_height * _width;i++)
	{
        B[i] /= (_height * _width);
        G[i] /= (_height * _width);
        R[i] /= (_height * _width);
    }

	for(int i=0;i<_height * _width;i++)
	{
		if(B[i] >= 255)
        	*(uint8_t*) &_image[i*_channel  ] = 255;
		else if(B[i] <= 0)
			*(uint8_t*) &_image[i*_channel  ] = 0;
		else
			*(uint8_t*) &_image[i*_channel  ] = uint8_t(B[i]);
		
		if(G[i] >= 255)
        	*(uint8_t*) &_image[i*_channel+1] = 255;
		else if(G[i] <= 0)
			*(uint8_t*) &_image[i*_channel+1] = 0;
		else
			*(uint8_t*) &_image[i*_channel+1] = uint8_t(G[i]);

		if(R[i] >= 255)
        	*(uint8_t*) &_image[i*_channel+2] = 255;
		else if(R[i] <= 0)
			*(uint8_t*) &_image[i*_channel+2] = 0;
		else
			*(uint8_t*) &_image[i*_channel+2] = uint8_t(R[i]);
    }

	fftw_destroy_plan(dft);
    fftw_destroy_plan(idft);
}
