#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <algorithm> 
#include <vector>

// define up scaling is true, down scaling is false
#define UP_SCAL   true
#define DOWN_SCAL false

using namespace std;

typedef unsigned char  BYTE;
typedef unsigned short WORD;
typedef unsigned int   DWORD;

class BMP
{
	public:
		
		// variable
		char* 	_header;
		WORD	_identifier;			// 0x0000
		DWORD	_file_size;				// 0x0002
		DWORD	_reserved;				// 0x0006
		DWORD	_bitmap_data_offset;	// 0x000A
		DWORD	_header_size;			// 0x000E
		int		_width;					// 0x0012
	    int		_height;				// 0x0016
	    WORD	_planes;				// 0x001A
	    WORD	_bits_per_pixel;		// 0x001C
	    DWORD	_compression;			// 0x001E
	    DWORD	_bitmap_data_size;		// 0x0022
	    DWORD	_h_resolution;			// 0x0026
	    DWORD	_v_resolution;			// 0x002A
	    DWORD	_used_colors;			// 0x002E
	    DWORD	_important_colors;		// 0x0032
	    int		_channel;
	    
	    char* 	_image;
			    
	    char*	_out_header;
	    char*	_out_image;
	    
	    // function
		bool read(const string filename);
		void show();
		bool write(const string filename);
		void flip();
		void quan(int factor);
		BYTE bilinear(int i, int j, int k);
		void scal(bool up);
		void bruteForce(BYTE factor);
		void gammaCorrection(double gamma);
		void repliPadding(int i, int j, int k, BYTE PAD[3][3]);
		BYTE meanFilter(int i, int j, int k, BYTE PAD[3][3]);
		BYTE medianFilter(int i, int j, int k, BYTE PAD[3][3]);
		int  gaussianFilter(int i, int j, int k, BYTE PAD[3][3]);
		int  laplacianFilter(int i, int j, int k, BYTE PAD[3][3]);
		void sharp();
		void denoise();
};

bool BMP::read(const string filename)
{
	// file open
	fstream in_bmp(filename, ios::in | ios::binary);
	if(!in_bmp.is_open())
		cerr << "BMP Open Fail\n";
	
	// read header
	_header = new char[54];
	in_bmp.read(_header, sizeof(BYTE) * 54);		
	_identifier 		= *(WORD*)  &_header[0];
	_file_size 			= *(DWORD*) &_header[2];
	_reserved 			= *(DWORD*) &_header[6];
	_bitmap_data_offset = *(DWORD*) &_header[10];
	_header_size 		= *(DWORD*) &_header[14];
	_width 				= *(int*)   &_header[18];
    _height 			= *(int*)   &_header[22];
    _planes 			= *(WORD*)  &_header[26];
    _bits_per_pixel 	= *(WORD*)  &_header[28];
    _compression 		= *(DWORD*) &_header[30];
    _bitmap_data_size 	= *(DWORD*) &_header[34];
    _h_resolution 		= *(DWORD*) &_header[38];
    _v_resolution 		= *(DWORD*) &_header[42];
    _used_colors 		= *(DWORD*) &_header[46];
    _important_colors 	= *(DWORD*) &_header[50];
    
    // calculate the number of channels
    _channel			= _bits_per_pixel / 8;
    
    // read image
    _image = new char[_height * _width * _channel];
	in_bmp.read(_image, sizeof(BYTE) * _height * _width * _channel);
	
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
	out_bmp.write(_out_header, sizeof(BYTE) * 54);
	out_bmp.write(_out_image , sizeof(BYTE) * (*(int*) &_out_header[22]) * (*(int*) &_out_header[18]) * _channel);

	// file close
	out_bmp.close();
	return true;
}

void BMP::flip()
{
	_out_header = _header;
	_out_image = new char[_height * _width * _channel];
	
	// flip the image horizontally
	for(int i=0;i<_height;i++)
		for(int j=0;j<_width;j++)
			for(int k=0;k<_channel;k++)
				_out_image[i*_width*_channel+j*_channel+k] = _image[i*_width*_channel+(_width-1-j)*_channel+k];
}

void BMP::quan(int factor)
{
	_out_header = _header;
	_out_image = _image;
	
	// quantization resolution
	for(int i=0;i<_height * _width * _channel;i++)
		_out_image[i] = (_out_image[i] >> factor) << factor;
}

BYTE BMP::bilinear(int i, int j, int k)
{
	double x, y, xRatio, yRatio;
	
	// get inverse ratio, - 1.0 to avoid x0, x1, y0 and y1 from exceeding the boundaries
	xRatio = (_height - 1.0) / ((*(int*) &_out_header[22]) - 1.0);
	yRatio = (_width  - 1.0) / ((*(int*) &_out_header[18]) - 1.0);

	// get the original position where each new pixel point should be on the original image
	x = i * xRatio;
	y = j * yRatio;

	// get 4 nearest pixels of the original position
	int x0 = floor(x);
	int x1 = ceil(x);
	int y0 = floor(y);
	int y1 = ceil(y);
	
	// bilinear interpolation
	float dx0 = x - x0;
	float dx1 = 1 - dx0;
	float dy0 = y - y0;
	float dy1 = 1 - dy0;
	
	// get image infomation of the 4 nearest pixels
	BYTE a0 = *(BYTE*) &_image[x0*_width*_channel+y0*_channel+k];
    BYTE a1 = *(BYTE*) &_image[x0*_width*_channel+y1*_channel+k];
    BYTE a2 = *(BYTE*) &_image[x1*_width*_channel+y0*_channel+k];
    BYTE a3 = *(BYTE*) &_image[x1*_width*_channel+y1*_channel+k];
    
    // calculate by bilinear interpolation
    // first calculate in direction of y
    BYTE b0 = a0 * dy1 + a1 * dy0;       
    BYTE b1 = a2 * dy1 + a3 * dy0;
    
    // then calculate in direction of x
    return b0 * dx1 + b1 * dx0;
}

void BMP::scal(bool up)
{
	_out_header = _header;
	
	// scaling 
	
	// up scaling
	if(up)
	{
		*(int*) &_out_header[22] = _height * 1.5;
		*(int*) &_out_header[18] = _width  * 1.5;
	}
	
	// down scaling
	else
	{
		*(int*) &_out_header[22] = _height / 1.5;
		*(int*) &_out_header[18] = _width  / 1.5;
	} 
	
	// make the width a multiple of 4
	*(int*) &_out_header[18] = (*(int*) &_out_header[18]) / 4 * 4;
	
	// resize output image
	_out_image = new char[ (*(int*) &_out_header[22]) * (*(int*) &_out_header[18]) * _channel];

	// update the value of file size
	*(int*) &_out_header[2] = (*(int*) &_out_header[22]) * (*(int*) &_out_header[18]) * 3+ 54;
	
	// perform bilinear interpolation according to each new pixel point
	for(int i=0;i<*(int*) &_out_header[22];i++)
		for(int j=0;j<*(int*) &_out_header[18];j++)
			for(int k=0;k<_channel;k++)
				*(BYTE*) &_out_image[i*(*(int*) &_out_header[18])*_channel+j*_channel+k] = bilinear(i,j,k);
}

void BMP::bruteForce(BYTE factor)
{
	_out_header = _header;
	_out_image = new char[ (*(int*) &_out_header[22]) * (*(int*) &_out_header[18]) * _channel];
	
	// luminosity enhancement
	for(int i=0;i<_height * _width * _channel;i++)
	{
		// to avoid going beyond the range
		if((int (*(BYTE*) &_image[i]) + factor) >= 255)
			*(BYTE*) &_out_image[i] = 255;
		else
			*(BYTE*) &_out_image[i] = (*(BYTE*) &_image[i]) + factor;
	}
}

void BMP::gammaCorrection(double gamma)
{
	_out_header = _header;
	_out_image = new char[ (*(int*) &_out_header[22]) * (*(int*) &_out_header[18]) * _channel];
	
	// luminosity enhancement
	for(int i=0;i<_height * _width * _channel;i++)
		*(BYTE*) &_out_image[i] = pow(int(*(BYTE*) &_image[i])/255.,gamma) * 255;
}

void BMP::repliPadding(int i, int j, int k, BYTE PAD[3][3])
{
	// get nine pixels need to convolution and do replication padding
	if(i == 0 && j == 0)
	{
		PAD[0][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[0][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[0][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
		PAD[1][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
		PAD[2][0] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][1] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][2] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
	}
	else if ((i == *(int*) &_out_header[22] - 1) && j == 0)
	{
		PAD[0][0] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[0][1] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[0][2] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
		PAD[1][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
		PAD[2][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
	}
	else if(i == 0 && (j == *(int*) &_out_header[18] - 1))
	{
		PAD[0][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[0][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[0][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[1][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][0] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[2][1] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][2] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
	}
	else if((i == *(int*) &_out_header[22] - 1) && (j == *(int*) &_out_header[18] - 1))
	{
		PAD[0][0] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[0][1] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[0][2] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[1][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[2][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
	}
	else if (i == 0)
	{
		PAD[0][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[0][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[0][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
		PAD[1][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[1][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
		PAD[2][0] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[2][1] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][2] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
	}
	else if (i == *(int*) &_out_header[22] - 1)
	{
		PAD[0][0] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[0][1] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[0][2] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
		PAD[1][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[1][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
		PAD[2][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[2][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
	}
	else if(j == 0)
	{
		PAD[0][0] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[0][1] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[0][2] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
		PAD[1][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
		PAD[2][0] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][1] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][2] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
	}
	else if (j == *(int*) &_out_header[18] - 1)
	{
		PAD[0][0] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[0][1] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[0][2] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[1][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][0] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[2][1] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][2] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
	}
	else
	{
		PAD[0][0] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[0][1] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[0][2] = *(BYTE*) &_image[(i-1)*(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
		PAD[1][0] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[1][1] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[1][2] = *(BYTE*) &_image[ i   *(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
		PAD[2][0] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+(j-1)*_channel+k];
		PAD[2][1] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+ j   *_channel+k];
		PAD[2][2] = *(BYTE*) &_image[(i+1)*(*(int*) &_out_header[18])*_channel+(j+1)*_channel+k];
	}
}

BYTE BMP::meanFilter(int i, int j, int k, BYTE PAD[3][3])
{
	return (PAD[0][0] + PAD[0][1] + PAD[0][2] + 
			PAD[1][0] + PAD[1][1] + PAD[1][2] +
			PAD[2][0] + PAD[2][1] + PAD[2][2] ) / 9;
}

BYTE BMP::medianFilter(int i, int j, int k, BYTE PAD[3][3])
{
	vector <BYTE> tmp;
	for(int a=0;a<3;a++)
		for(int b=0;b<3;b++)
			tmp.push_back(PAD[a][b]);
	
	// sort nine pixels
	sort(tmp.begin(),tmp.end());
	
	// get median
	return tmp[4];
}

int BMP::gaussianFilter(int i, int j, int k, BYTE PAD[3][3])
{
	// Gaussian filter
	return (PAD[0][0] * 1 + PAD[0][1] * 2 + PAD[0][2] * 1 + 
			PAD[1][0] * 2 + PAD[1][1] * 4 + PAD[1][2] * 2 +
			PAD[2][0] * 1 + PAD[2][1] * 2 + PAD[2][2] * 1 ) / 16;
}

int BMP::laplacianFilter(int i, int j, int k, BYTE PAD[3][3])
{
	// Laplacian filter
	return 	PAD[0][0] * (-1) + PAD[0][1] * (-1) + PAD[0][2] * (-1) + 
			PAD[1][0] * (-1) + PAD[1][1] * ( 9) + PAD[1][2] * (-1) +
			PAD[2][0] * (-1) + PAD[2][1] * (-1) + PAD[2][2] * (-1) ;
}

void BMP::sharp()
{
	_out_header = _header;
	_out_image = new char[ (*(int*) &_out_header[22]) * (*(int*) &_out_header[18]) * _channel];
	
	// tmp is used to store convolution result
	int tmp;
	
	// PAD is used to store the nine pixels for the current convolution operation
	BYTE PAD[3][3];
	
	// convolution
	for(int i=0;i<*(int*) &_out_header[22];i++)
	{
		for(int j=0;j<*(int*) &_out_header[18];j++)
		{
			for(int k=0;k<_channel;k++)
			{
				// tmp initialization
				tmp = 0;
				
				// replication padding
				repliPadding(i,j,k,PAD);
				
				// convolution
				// Laplacian filter
				tmp = laplacianFilter(i,j,k,PAD);
						
				// if convolution result goes beyond the range, let it equal to 255 or 0
				if(tmp > 255)
					tmp = 255;
				else if(tmp < 0)
					tmp = 0;
				*(BYTE*) &_out_image[i*(*(int*) &_out_header[18])*_channel+j*_channel+k] = (BYTE) tmp;
			}
		}
	}
	
	// update
	_image = _out_image;
}

void BMP::denoise()
{
	_out_header = _header;
	_out_image = new char[ (*(int*) &_out_header[22]) * (*(int*) &_out_header[18]) * _channel];
	
	// PAD is used to store the nine pixels for the current convolution operation
	BYTE PAD[3][3]; 
	
	// tmp is used to store convolution result
	int tmp;
	
	// convolution
	for(int i=0;i<*(int*) &_out_header[22];i++)
	{
		for(int j=0;j<*(int*) &_out_header[18];j++)
		{
			for(int k=0;k<_channel;k++)
			{
				// tmp initialization
				tmp = 0;
				
				// replication padding
				repliPadding(i,j,k,PAD);
				
				// convolution mean filter
				*(BYTE*) &_out_image[i*(*(int*) &_out_header[18])*_channel+j*_channel+k] = meanFilter(i,j,k,PAD);
			
				// convolution median filter
				//*(BYTE*) &_out_image[i*(*(int*) &_out_header[18])*_channel+j*_channel+k] = medianFilter(i,j,k,PAD);
				
				// convolution Gaussian Filter
				//tmp = gaussianFilter(i,j,k,PAD);
				
				// if convolution result goes beyond the range, let it equal to 255 or 0
				/*if(tmp > 255)
					tmp = 255;
				else if(tmp < 0)
					tmp = 0;
					
				*(BYTE*) &_out_image[i*(*(int*) &_out_header[18])*_channel+j*_channel+k] = tmp; */ 
			} 
		}
	}
	
	// update
	_image = _out_image;
}

