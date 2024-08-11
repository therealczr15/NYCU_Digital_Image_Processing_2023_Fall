// define up scaling is true, down scaling is false
#define UP_SCAL   true
#define DOWN_SCAL false

// define _USE_MATH_DEFINES to use M_PI
#define _USE_MATH_DEFINES

// include header files
#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <algorithm> 
#include <vector>
#include <climits>

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
	    uint8_t*    _Y;
	    uint8_t*	_Cr;
	    uint8_t*	_Cb;
	    
	    // function
		bool    	read(const string filename);
		void    	show();
		bool    	write(const string filename);
		void    	flip();
		void    	quan(int factor);
		uint8_t 	bilinear(int i, int j, int k, int oH, int oW);
		void    	scal(bool up);
		void    	YCrCb2RGB();
		void 		RGB2YCrCb();
		void		HSI2RGB(double* H, double* S, double* I);
		void		RGB2HSI(double* H, double* S, double* I);
		void    	brightnessEnhance(double alpha, double beta);
		void    	gammaCorrection(double gamma);
		void    	histEqual();
		void    	convMask3      		(int i, int j, uint8_t mask[3][3]);
		uint8_t 	laplacianAFilter	(uint8_t mask[3][3]);
		uint8_t 	laplacianBFilter	(uint8_t mask[3][3]);
		uint8_t		prewittXFilter 		(uint8_t mask[3][3]);
		uint8_t		prewittYFilter 		(uint8_t mask[3][3]);
		uint8_t		sobelXFilter		(uint8_t mask[3][3]);
		uint8_t		sobelYFilter 		(uint8_t mask[3][3]);
		uint8_t 	gaussianFilter 		(uint8_t mask[3][3]);
		uint8_t 	medianFilter   		(uint8_t mask[3][3]);
		uint8_t 	meanFilter     		(uint8_t mask[3][3]);
		void    	sharp(int mode);
		void 		denoise(int mode);
		void 		maxRGB();
		void 		grayWorld();
		void		saturationEnhance(double gamma);
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
	
	// read Y, Cr, Cb
	RGB2YCrCb();
	
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

void BMP::flip()
{
	// declare tmpImg = original image
	char* tmpImg = new char[_height * _width * _channel];
	for(int i=0;i<_height*_width*_channel;i++)
		*(uint8_t*) &tmpImg[i] = *(uint8_t*) &_image[i];
		
	// flip the image horizontally
	for(int i=0;i<_height;i++)
		for(int j=0;j<_width;j++)
			for(int k=0;k<_channel;k++)
				*(uint8_t*) &_image[i*_width*_channel+j*_channel+k] = *(uint8_t*) &tmpImg[i*_width*_channel+(_width-1-j)*_channel+k];
	
	// update Y Cr Cb info
	RGB2YCrCb();
	
	// free memory space
	delete []tmpImg;
}

void BMP::quan(int factor)
{
	// quantization resolution
	for(int i=0;i<_height * _width * _channel;i++)
		*(uint8_t*) &_image[i] = (*(uint8_t*) &_image[i] >> factor) << factor;
		
	// update Y Cr Cb info
	RGB2YCrCb();
}

uint8_t BMP::bilinear(int i, int j, int k, int oH, int oW)
{
	double x, y, xRatio, yRatio;
	
	// get inverse ratio, - 1.0 to avoid x0, x1, y0 and y1 from exceeding the boundaries
	xRatio = (oH - 1.0) / (_height - 1.0);
	yRatio = (oW - 1.0) / (_width  - 1.0);

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
	uint8_t a0 = *(uint8_t*) &_image[x0*oW*_channel+y0*_channel+k];
    uint8_t a1 = *(uint8_t*) &_image[x0*oW*_channel+y1*_channel+k];
    uint8_t a2 = *(uint8_t*) &_image[x1*oW*_channel+y0*_channel+k];
    uint8_t a3 = *(uint8_t*) &_image[x1*oW*_channel+y1*_channel+k];
    
    // calculate by bilinear interpolation
    // first calculate in direction of y
    uint8_t b0 = a0 * dy1 + a1 * dy0;       
    uint8_t b1 = a2 * dy1 + a3 * dy0;
    
    // then calculate in direction of x
    return b0 * dx1 + b1 * dx0;
}

void BMP::scal(bool up)
{
	// store original height and width
	int originalH = _height;
	int originalW = _width;
	
	// scaling 
	
	// up scaling
	if(up)
	{
		*(uint32_t*) &_header[22] = _height * 1.5;
		*(uint32_t*) &_header[18] = _width  * 1.5;
		_height *= 1.5;
		_width  *= 1.5;
	}
	
	// down scaling
	else
	{
		*(uint32_t*) &_header[22] = _height / 1.5;
		*(uint32_t*) &_header[18] = _width  / 1.5;
		_height /= 1.5;
		_width  /= 1.5;
	} 
	
	// make the width a multiple of 4
	*(uint32_t*) &_header[18] = _width / 4 * 4;
	_width = _width / 4 * 4;
	
	// declare tmpImg to get new image
	char* tmpImg = new char[_height * _width * _channel];

	// update the value of file size
	*(uint32_t*) &_header[2] = _height * _width * 3 + 54;
	
	// perform bilinear interpolation according to each new pixel point and store in tmpImg
	for(int i=0;i<_height;i++)
		for(int j=0;j<_width;j++)
			for(int k=0;k<_channel;k++)
				*(uint8_t*) &tmpImg[i*_width*_channel+j*_channel+k] = bilinear(i,j,k,originalH,originalW);
	
	// resize the original image
	_image = new char[_height * _width * _channel];
	for(int i=0;i<_height*_width*_channel;i++)
		*(uint8_t*) &_image[i] = *(uint8_t*) &tmpImg[i];
	
	// update Y Cr Cb info
	RGB2YCrCb();
	
	// free memory space
	delete []tmpImg;
}

void BMP::YCrCb2RGB()
{
	for(int i=0;i<_height*_width;i++)
	{
		// to avoid B going beyond the range
		if(*(uint8_t*) &_Y[i] + 1.773 * (*(uint8_t*) &_Cb[i] - 128) >= 255)
			*(uint8_t*) &_image[i*_channel  ] = 255;
		else if(*(uint8_t*) &_Y[i] + 1.773 * (*(uint8_t*) &_Cb[i] - 128) <= 0)
			*(uint8_t*) &_image[i*_channel  ] = 0;
		else
			*(uint8_t*) &_image[i*_channel  ] = *(uint8_t*) &_Y[i] + 1.773 * (*(uint8_t*) &_Cb[i] - 128);
			
		// to avoid G going beyond the range
		if(*(uint8_t*) &_Y[i] - 0.714 * (*(uint8_t*) &_Cr[i] - 128) - 0.344 * (*(uint8_t*) &_Cb[i] - 128) >= 255)
			*(uint8_t*) &_image[i*_channel+1] = 255;
		else if(*(uint8_t*) &_Y[i] - 0.714 * (*(uint8_t*) &_Cr[i] - 128) - 0.344 * (*(uint8_t*) &_Cb[i] - 128) <= 0)
			*(uint8_t*) &_image[i*_channel+1] = 0;
		else
			*(uint8_t*) &_image[i*_channel+1] = *(uint8_t*) &_Y[i] - 0.714 * (*(uint8_t*) &_Cr[i] - 128) - 0.344 * (*(uint8_t*) &_Cb[i] - 128);
			
		// to avoid R going beyond the range
		if(*(uint8_t*) &_Y[i] + 1.403 * (*(uint8_t*) &_Cr[i] - 128) >= 255)
			*(uint8_t*) &_image[i*_channel+2] = 255;
		else if(*(uint8_t*) &_Y[i] + 1.403 * (*(uint8_t*) &_Cr[i] - 128) <= 0)
			*(uint8_t*) &_image[i*_channel+2] = 0;
		else
			*(uint8_t*) &_image[i*_channel+2] = *(uint8_t*) &_Y[i] + 1.403 * (*(uint8_t*) &_Cr[i] - 128);
	}
}

void BMP::RGB2YCrCb()
{
	_Y  = new uint8_t[_height*_width];
	_Cr = new uint8_t[_height*_width];
	_Cb = new uint8_t[_height*_width];
	
	for(int i=0;i<_height * _width;i++)
	{
		*(uint8_t*) &_Y [i] =  *(uint8_t*) &_image[i*_channel  ] * 0.114 + *(uint8_t*) &_image[i*_channel+1] * 0.587 + *(uint8_t*) &_image[i*_channel+2] * 0.299;
		*(uint8_t*) &_Cr[i] = (*(uint8_t*) &_image[i*_channel+2] - *(uint8_t*) &_Y [i]) * 0.713 + 128;
		*(uint8_t*) &_Cb[i] = (*(uint8_t*) &_image[i*_channel  ] - *(uint8_t*) &_Y [i]) * 0.564 + 128;
	}
}

void BMP::HSI2RGB(double* H, double* S, double* I)
{
	// turn HSI to RGB
	for(int i=0;i<_height*_width;i++)
	{
		double b, g, r; 
		if(0 <= H[i] && H[i] < 2.0 * M_PI / 3.0)
		{
			b = (1.0 - S[i]) * I[i];
			r = (1.0 + S[i] * cos(H[i]) / cos(M_PI / 3.0 - H[i]) ) * I[i];
			g = 3 * I[i] - (r + b);
		}
		else if(2.0 * M_PI / 3.0 <= H[i] && H[i] < 4.0 * M_PI / 3.0)
		{
			r = (1.0 - S[i]) * I[i];
			g = (1.0 + S[i] * cos(H[i] - 2.0 * M_PI / 3.0) / cos(M_PI - H[i]) ) * I[i];
			b = 3 * I[i] - (r + g);
		}
		else if(4.0 * M_PI / 3.0 <= H[i] && H[i] < 2.0 * M_PI)
		{
			g = (1.0 - S[i]) * I[i];
			b = (1.0 + S[i] * cos(H[i] - 4.0 * M_PI / 3.0) / cos(5.0 * M_PI / 3.0 - H[i]) ) * I[i];
			r = 3 * I[i] - (g + b);
		}
		
		if(b >= 255)
			b = 255;
		else if(b <= 0)
			b = 0;
		
		if(g >= 255)
			g = 255;
		else if(g <= 0)
			g = 0;
			
		if(r >= 255)
			r = 255;
		else if(r <= 0)
			r = 0.0;
		
		*(uint8_t*) &_image[i*_channel  ] = uint8_t(b);
		*(uint8_t*) &_image[i*_channel+1] = uint8_t(g);
		*(uint8_t*) &_image[i*_channel+2] = uint8_t(r);  
		
	}
}

void BMP::RGB2HSI(double* H, double* S, double* I)
{	
	// turn RGB to HSI
	for(int i=0;i<_height*_width;i++)
	{
		uint8_t b   = *(uint8_t*) &_image[i*_channel  ];
		uint8_t g   = *(uint8_t*) &_image[i*_channel+1];
		uint8_t r   = *(uint8_t*) &_image[i*_channel+2];
		
		uint8_t min = *(uint8_t*) &_image[i*_channel  ];
		if(g < min)
			min = g;
		if(r < min)
			min = r;
			
		double theta = acos( 0.5 * (r-g+r-b) / sqrt( pow(r-g,2) + (r-b) * (g-b) ) );
		
		I[i] = (r + g + b) / 3.0;
		if(I[i] == 0)
			S[i] = 0;
		else
			S[i] = 1.0 - 3.0 * min / (r + g + b);
		H[i] = (b <= g) ? theta : 2 * M_PI - theta;
	}
}

void BMP::brightnessEnhance(double alpha, double beta)
{
	// brightness enhancement
	for(int i=0;i<_height * _width;i++)
	{
		// to avoid going beyond the range
		if(*(uint8_t*) &_Y[i] * alpha + beta >= 255)
			*(uint8_t*) &_Y[i] = 255;
		else if(*(uint8_t*) &_Y[i] * alpha + beta <= 0)
			*(uint8_t*) &_Y[i] = 0;
		else
			*(uint8_t*) &_Y[i] = *(uint8_t*) &_Y[i] * alpha + beta;
	}
	
	// update image
	YCrCb2RGB();
}

void BMP::gammaCorrection(double gamma)
{
	// gamma correction
	for(int i=0;i<_height*_width;i++)
		*(uint8_t*) &_Y[i] = pow(int(*(uint8_t*) &_Y[i])/255.,gamma) * 255;

	// update image
	YCrCb2RGB();
}

void BMP::histEqual()
{
	// declare histogram array
	int hist[256] = {0};
	
	// iterate through the array, calculate the quantity of 256 levels, and store it in the histogram
	for(int i=0;i<_height*_width;i++)
		hist[int(*(uint8_t*) &_Y[i])]++;
	
	// declare probability array 
	double prob[256] = {0};

	// iterate through the array, calculate the probability of 256 levels, and store it in the histogram
	for(int i=0;i<256;i++)
		prob[i] = hist[i] / double (_height * _width) ;
	
	// declare cdf array 
	double cdf[256] = {0};

	// iterate through the array, calculate the cdf of 256 levels by probability histogram, and store it in the cdf histogram
	for(int i=0;i<256;i++)
	{
		if(i == 0)
			cdf[i] = prob[i];
		else
			cdf[i] = cdf[i-1] + prob[i];
	}
	
	// calculate the result of every pixel in the three channels
	for(int i=0;i<_height * _width;i++)
	{
		if(255 * cdf[int (*(uint8_t*) &_Y[i])] >= 255)
			*(uint8_t*) &_Y[i] = 255;
		else if(255 * cdf[int (*(uint8_t*) &_Y[i])] <= 0)
			*(uint8_t*) &_Y[i] = 0;
		else
			*(uint8_t*) &_Y[i] = (255 * cdf[int (*(uint8_t*) &_Y[i])]);
	} 
	
	// update image
	YCrCb2RGB();
}

void BMP::convMask3(int i, int j, uint8_t mask[3][3])
{
	// get nine pixels need to convolution and do replication padding
	if(i == 0 && j == 0)
	{
		mask[0][0] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[0][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[0][2] = *(uint8_t*) &_Y[ i   *_width+(j+1)];
		mask[1][0] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[1][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[1][2] = *(uint8_t*) &_Y[ i   *_width+(j+1)];
		mask[2][0] = *(uint8_t*) &_Y[(i+1)*_width+ j   ];
		mask[2][1] = *(uint8_t*) &_Y[(i+1)*_width+ j   ];
		mask[2][2] = *(uint8_t*) &_Y[(i+1)*_width+(j+1)];
	}
	else if ((i == _height - 1) && j == 0)
	{
		mask[0][0] = *(uint8_t*) &_Y[(i-1)*_width+ j   ];
		mask[0][1] = *(uint8_t*) &_Y[(i-1)*_width+ j   ];
		mask[0][2] = *(uint8_t*) &_Y[(i-1)*_width+(j+1)];
		mask[1][0] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[1][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[1][2] = *(uint8_t*) &_Y[ i   *_width+(j+1)];
		mask[2][0] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[2][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[2][2] = *(uint8_t*) &_Y[ i   *_width+(j+1)];
	}
	else if(i == 0 && (j == _width - 1))
	{
		mask[0][0] = *(uint8_t*) &_Y[ i   *_width+(j-1)];
		mask[0][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[0][2] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[1][0] = *(uint8_t*) &_Y[ i   *_width+(j-1)];
		mask[1][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[1][2] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[2][0] = *(uint8_t*) &_Y[(i+1)*_width+(j-1)];
		mask[2][1] = *(uint8_t*) &_Y[(i+1)*_width+ j   ];
		mask[2][2] = *(uint8_t*) &_Y[(i+1)*_width+ j   ];
	}
	else if((i == _height - 1) && (j == _width - 1))
	{
		mask[0][0] = *(uint8_t*) &_Y[(i-1)*_width+(j-1)];
		mask[0][1] = *(uint8_t*) &_Y[(i-1)*_width+ j   ];
		mask[0][2] = *(uint8_t*) &_Y[(i-1)*_width+ j   ];
		mask[1][0] = *(uint8_t*) &_Y[ i   *_width+(j-1)];
		mask[1][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[1][2] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[2][0] = *(uint8_t*) &_Y[ i   *_width+(j-1)];
		mask[2][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[2][2] = *(uint8_t*) &_Y[ i   *_width+ j   ];
	}
	else if (i == 0)
	{
		mask[0][0] = *(uint8_t*) &_Y[ i   *_width+(j-1)];
		mask[0][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[0][2] = *(uint8_t*) &_Y[ i   *_width+(j+1)];
		mask[1][0] = *(uint8_t*) &_Y[ i   *_width+(j-1)];
		mask[1][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[1][2] = *(uint8_t*) &_Y[ i   *_width+(j+1)];
		mask[2][0] = *(uint8_t*) &_Y[(i+1)*_width+(j-1)];
		mask[2][1] = *(uint8_t*) &_Y[(i+1)*_width+ j   ];
		mask[2][2] = *(uint8_t*) &_Y[(i+1)*_width+(j+1)];
	}
	else if (i == _height - 1)
	{
		mask[0][0] = *(uint8_t*) &_Y[(i-1)*_width+(j-1)];
		mask[0][1] = *(uint8_t*) &_Y[(i-1)*_width+ j   ];
		mask[0][2] = *(uint8_t*) &_Y[(i-1)*_width+(j+1)];
		mask[1][0] = *(uint8_t*) &_Y[ i   *_width+(j-1)];
		mask[1][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[1][2] = *(uint8_t*) &_Y[ i   *_width+(j+1)];
		mask[2][0] = *(uint8_t*) &_Y[ i   *_width+(j-1)];
		mask[2][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[2][2] = *(uint8_t*) &_Y[ i   *_width+(j+1)];
	}
	else if(j == 0)
	{
		mask[0][0] = *(uint8_t*) &_Y[(i-1)*_width+ j   ];
		mask[0][1] = *(uint8_t*) &_Y[(i-1)*_width+ j   ];
		mask[0][2] = *(uint8_t*) &_Y[(i-1)*_width+(j+1)];
		mask[1][0] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[1][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[1][2] = *(uint8_t*) &_Y[ i   *_width+(j+1)];
		mask[2][0] = *(uint8_t*) &_Y[(i+1)*_width+ j   ];
		mask[2][1] = *(uint8_t*) &_Y[(i+1)*_width+ j   ];
		mask[2][2] = *(uint8_t*) &_Y[(i+1)*_width+(j+1)];
	}
	else if (j == _width - 1)
	{
		mask[0][0] = *(uint8_t*) &_Y[(i-1)*_width+(j-1)];
		mask[0][1] = *(uint8_t*) &_Y[(i-1)*_width+ j   ];
		mask[0][2] = *(uint8_t*) &_Y[(i-1)*_width+ j   ];
		mask[1][0] = *(uint8_t*) &_Y[ i   *_width+(j-1)];
		mask[1][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[1][2] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[2][0] = *(uint8_t*) &_Y[(i+1)*_width+(j-1)];
		mask[2][1] = *(uint8_t*) &_Y[(i+1)*_width+ j   ];
		mask[2][2] = *(uint8_t*) &_Y[(i+1)*_width+ j   ];
	}
	else
	{
		mask[0][0] = *(uint8_t*) &_Y[(i-1)*_width+(j-1)];
		mask[0][1] = *(uint8_t*) &_Y[(i-1)*_width+ j   ];
		mask[0][2] = *(uint8_t*) &_Y[(i-1)*_width+(j+1)];
		mask[1][0] = *(uint8_t*) &_Y[ i   *_width+(j-1)];
		mask[1][1] = *(uint8_t*) &_Y[ i   *_width+ j   ];
		mask[1][2] = *(uint8_t*) &_Y[ i   *_width+(j+1)];
		mask[2][0] = *(uint8_t*) &_Y[(i+1)*_width+(j-1)];
		mask[2][1] = *(uint8_t*) &_Y[(i+1)*_width+ j   ];
		mask[2][2] = *(uint8_t*) &_Y[(i+1)*_width+(j+1)];
	}
}

uint8_t BMP::laplacianAFilter(uint8_t mask[3][3])
{
	// Laplacian A filter
	if		(	mask[0][0] * ( 0) + mask[0][1] * (-1) + mask[0][2] * ( 0) + 
				mask[1][0] * (-1) + mask[1][1] * ( 5) + mask[1][2] * (-1) +
				mask[2][0] * ( 0) + mask[2][1] * (-1) + mask[2][2] * ( 0) >= 255)
		return 255;
	else if	(	mask[0][0] * ( 0) + mask[0][1] * (-1) + mask[0][2] * ( 0) + 
				mask[1][0] * (-1) + mask[1][1] * ( 5) + mask[1][2] * (-1) +
				mask[2][0] * ( 0) + mask[2][1] * (-1) + mask[2][2] * ( 0) <= 0)
		return 0;
	else
		return	mask[0][0] * ( 0) + mask[0][1] * (-1) + mask[0][2] * ( 0) + 
				mask[1][0] * (-1) + mask[1][1] * ( 5) + mask[1][2] * (-1) +
				mask[2][0] * ( 0) + mask[2][1] * (-1) + mask[2][2] * ( 0) ;
}

uint8_t BMP::laplacianBFilter(uint8_t mask[3][3])
{
	// Laplacian B filter
	if		(	mask[0][0] * (-1) + mask[0][1] * (-1) + mask[0][2] * (-1) + 
				mask[1][0] * (-1) + mask[1][1] * ( 9) + mask[1][2] * (-1) +
				mask[2][0] * (-1) + mask[2][1] * (-1) + mask[2][2] * (-1) >= 255)
		return 255;
	else if	(	mask[0][0] * (-1) + mask[0][1] * (-1) + mask[0][2] * (-1) + 
				mask[1][0] * (-1) + mask[1][1] * ( 9) + mask[1][2] * (-1) +
				mask[2][0] * (-1) + mask[2][1] * (-1) + mask[2][2] * (-1) <= 0)
		return 0;
	else
		return	mask[0][0] * (-1) + mask[0][1] * (-1) + mask[0][2] * (-1) + 
				mask[1][0] * (-1) + mask[1][1] * ( 9) + mask[1][2] * (-1) +
				mask[2][0] * (-1) + mask[2][1] * (-1) + mask[2][2] * (-1) ;
}

uint8_t BMP::prewittXFilter(uint8_t mask[3][3])
{
	// Prewitt X filter
	if		(	mask[0][0] * ( 1) + mask[0][1] * ( 1) + mask[0][2] * ( 1) + 
				mask[1][0] * ( 0) + mask[1][1] * ( 1) + mask[1][2] * ( 0) +
				mask[2][0] * (-1) + mask[2][1] * (-1) + mask[2][2] * (-1) >= 255)
		return 255;
	else if	(	mask[0][0] * ( 1) + mask[0][1] * ( 1) + mask[0][2] * ( 1) + 
				mask[1][0] * ( 0) + mask[1][1] * ( 1) + mask[1][2] * ( 0) +
				mask[2][0] * (-1) + mask[2][1] * (-1) + mask[2][2] * (-1) <= 0)
		return 0;
	else
		return	mask[0][0] * ( 1) + mask[0][1] * ( 1) + mask[0][2] * ( 1) + 
				mask[1][0] * ( 0) + mask[1][1] * ( 1) + mask[1][2] * ( 0) +
				mask[2][0] * (-1) + mask[2][1] * (-1) + mask[2][2] * (-1) ;
}

uint8_t BMP::prewittYFilter(uint8_t mask[3][3])
{
	// Prewitt Y filter
	if		(	mask[0][0] * (-1) + mask[0][1] * ( 0) + mask[0][2] * ( 1) + 
				mask[1][0] * (-1) + mask[1][1] * ( 1) + mask[1][2] * ( 1) +
				mask[2][0] * (-1) + mask[2][1] * ( 0) + mask[2][2] * ( 1) >= 255)
		return 255;
	else if	(	mask[0][0] * (-1) + mask[0][1] * ( 0) + mask[0][2] * ( 1) + 
				mask[1][0] * (-1) + mask[1][1] * ( 1) + mask[1][2] * ( 1) +
				mask[2][0] * (-1) + mask[2][1] * ( 0) + mask[2][2] * ( 1) <= 0)
		return 0;
	else
		return	mask[0][0] * (-1) + mask[0][1] * ( 0) + mask[0][2] * ( 1) + 
				mask[1][0] * (-1) + mask[1][1] * ( 1) + mask[1][2] * ( 1) +
				mask[2][0] * (-1) + mask[2][1] * ( 0) + mask[2][2] * ( 1) ;
}

uint8_t BMP::sobelXFilter(uint8_t mask[3][3])
{
	// Sobel X filter
	if		(	mask[0][0] * ( 1) + mask[0][1] * ( 2) + mask[0][2] * ( 1) + 
				mask[1][0] * ( 0) + mask[1][1] * ( 1) + mask[1][2] * ( 0) +
				mask[2][0] * (-1) + mask[2][1] * (-2) + mask[2][2] * (-1) >= 255)
		return 255;
	else if	(	mask[0][0] * ( 1) + mask[0][1] * ( 2) + mask[0][2] * ( 1) + 
				mask[1][0] * ( 0) + mask[1][1] * ( 1) + mask[1][2] * ( 0) +
				mask[2][0] * (-1) + mask[2][1] * (-2) + mask[2][2] * (-1) <= 0)
		return 0;
	else
		return	mask[0][0] * ( 1) + mask[0][1] * ( 2) + mask[0][2] * ( 1) + 
				mask[1][0] * ( 0) + mask[1][1] * ( 1) + mask[1][2] * ( 0) +
				mask[2][0] * (-1) + mask[2][1] * (-2) + mask[2][2] * (-1) ;
}

uint8_t BMP::sobelYFilter(uint8_t mask[3][3])
{
	// Sobel Y filter
	if		(	mask[0][0] * (-1) + mask[0][1] * ( 0) + mask[0][2] * ( 1) + 
				mask[1][0] * (-2) + mask[1][1] * ( 1) + mask[1][2] * ( 2) +
				mask[2][0] * (-1) + mask[2][1] * ( 0) + mask[2][2] * ( 1) >= 255)
		return 255;
	else if	(	mask[0][0] * (-1) + mask[0][1] * ( 0) + mask[0][2] * ( 1) + 
				mask[1][0] * (-2) + mask[1][1] * ( 1) + mask[1][2] * ( 2) +
				mask[2][0] * (-1) + mask[2][1] * ( 0) + mask[2][2] * ( 1) <= 0)
		return 0;
	else
		return	mask[0][0] * (-1) + mask[0][1] * ( 0) + mask[0][2] * ( 1) + 
				mask[1][0] * (-2) + mask[1][1] * ( 1) + mask[1][2] * ( 2) +
				mask[2][0] * (-1) + mask[2][1] * ( 0) + mask[2][2] * ( 1) ;
}

uint8_t BMP::gaussianFilter(uint8_t mask[3][3])
{
	// Gaussian filter
	return (mask[0][0] * (1) + mask[0][1] * (2) + mask[0][2] * (1) + 
			mask[1][0] * (2) + mask[1][1] * (4) + mask[1][2] * (2) +
			mask[2][0] * (1) + mask[2][1] * (2) + mask[2][2] * (1) ) / 16;
}

uint8_t BMP::medianFilter(uint8_t mask[3][3])
{
	// store value of nine pixels in a vector
	vector <uint8_t> tmp;
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			tmp.push_back(mask[i][j]);
			
	// sort nine pixels
	sort(tmp.begin(),tmp.end());
	
	// get median
	return tmp[4];	
}

uint8_t BMP::meanFilter(uint8_t mask[3][3])
{
	// Mean filter 
	return (mask[0][0] * (1) + mask[0][1] * (1) + mask[0][2] * (1) + 
			mask[1][0] * (1) + mask[1][1] * (1) + mask[1][2] * (1) +
			mask[2][0] * (1) + mask[2][1] * (1) + mask[2][2] * (1) ) / 9;
}

void BMP::sharp(int mode)
{
	// mask is used to store the nine pixels for the current convolution operation
	uint8_t mask[3][3];
	
	// declare tmpY to caluculate new result of Y
	uint8_t* tmpY = new uint8_t[_height*_width];
	
	// convolution
	for(int i=0;i<_height;i++)
	{
		for(int j=0;j<_width;j++)
		{
			
			// get the 3*3 convolution mask
			convMask3(i,j,mask);
			
			// convolution
			
			// use different filters in different modes
			if(mode == 1)
				// Laplacian A filter
				tmpY[i*_width+j] = laplacianAFilter(mask);
			else if(mode == 2)
				// Laplacian B filter
				tmpY[i*_width+j] = laplacianBFilter(mask);
			else if(mode == 3)
				// Prewitt X filter
				tmpY[i*_width+j] = prewittXFilter(mask);
			else if(mode == 4)
				// Prewitt Y filter
				tmpY[i*_width+j] = prewittYFilter(mask);
			else if(mode == 5)
				// Sobel X filter
				tmpY[i*_width+j] = sobelXFilter(mask);
			else if(mode == 6)
				// Sobel Y filter
				tmpY[i*_width+j] = sobelYFilter(mask);
		}
	}
	
	// update Y
	for(int i=0;i<_height*_width;i++)
		*(uint8_t*) &_Y[i] = *(uint8_t*) &tmpY[i];
	
	// update image
	YCrCb2RGB();
	
	// free memory space
	delete []tmpY;
}

void BMP::denoise(int mode)
{
	// mask is used to store the nine pixels for the current convolution operation
	uint8_t mask[3][3];
	
	// declare tmpY to caluculate new result of Y
	uint8_t* tmpY = new uint8_t[_height*_width];
	
	// convolution
	for(int i=0;i<_height;i++)
	{
		for(int j=0;j<_width;j++)
		{
			// get the 3*3 convolution mask
			convMask3(i,j,mask);
			
			// convolution
			
			// use different filters in different modes
			if(mode == 1)
				// Gaussian filter
				tmpY[i*_width+j] = gaussianFilter(mask);
			else if(mode == 2)
				// Median filter
				tmpY[i*_width+j] = medianFilter(mask);
			else if(mode == 3)
				// Mean filter
				tmpY[i*_width+j] = meanFilter(mask);
		}
	}
	
	// update Y
	for(int i=0;i<_height*_width;i++)
		*(uint8_t*) &_Y[i] = *(uint8_t*) &tmpY[i];
	
	// update image
	YCrCb2RGB();
	
	// free memory space
	delete []tmpY;
}

void BMP::maxRGB()
{
	// declare maxB, maxG and maxR to store the maximum level of the three channel
	uint8_t maxB = 0;
	uint8_t maxG = 0;
	uint8_t maxR = 0;
	
	// find value of maxB, maxG amd maxR
	for(int i=0;i<_height*_width*_channel;i++)
	{
		if(i%3 == 0)
			maxB = *(uint8_t*) &_image[i] > maxB ? *(uint8_t*) &_image[i] : maxB;
		else if(i%3 == 1)
			maxG = *(uint8_t*) &_image[i] > maxG ? *(uint8_t*) &_image[i] : maxG;
		else
			maxR = *(uint8_t*) &_image[i] > maxR ? *(uint8_t*) &_image[i] : maxR;
	}
	
	// output information
	cout << "maxB: " << int (maxB) << "\t\tmaxG: " << int (maxG) << "\t\tmaxR: " << int (maxR) << '\n';
	
	// declare ratioB, ratioG, ratioR and calculate their value
	double ratioB, ratioG, ratioR;
	ratioB = 255. / int (maxB);
	ratioG = 255. / int (maxG);
	ratioR = 255. / int (maxR);
	
	// output information
	cout << "ratioB: " << ratioB << "\t\tratioG: " << ratioG << "\t\tratioR: " << ratioR << '\n';
	
	// multiply every pixel in the three channel with their corresponding ratio (maxRGB method)
	for(int i=0;i<_height*_width*_channel;i++)
	{
		if(i%3 == 0)
			*(uint8_t*) &_image[i] = *(uint8_t*) &_image[i] * ratioB;
		else if(i%3 == 1)
			*(uint8_t*) &_image[i] = *(uint8_t*) &_image[i] * ratioG;
		else
			*(uint8_t*) &_image[i] = *(uint8_t*) &_image[i] * ratioR;
	}	
	cout << "---------------------------------\n";
	
	// update Y Cr Cb info
	RGB2YCrCb();
}

void BMP::grayWorld()
{
	// declare sumB, sumG and sumR to calculate mean in the three channels
	double sumB = 0;
	double sumG = 0;
	double sumR = 0;
	
	// calculate sumB, sumG and sumR
	for(int i=0;i<_height*_width*_channel;i++)
	{
		if(i%3 == 0)
			sumB += *(uint8_t*) &_image[i];
		else if(i%3 == 1)
			sumG += *(uint8_t*) &_image[i];
		else
			sumR += *(uint8_t*) &_image[i];
	}
	
	// output information
	cout << "sumB: " << sumB << "\tsumG: " << sumG << "\tsumR: " << sumR << '\n';
	
	// declare mB, mG and mR, and calculate mean value of the mean in the three channels
	double mB, mG, mR;
	mB = sumB * 1.0 / (_height * _width);
	mG = sumG * 1.0 / (_height * _width);
	mR = sumR * 1.0 / (_height * _width);
	
	// output information
	cout << "mB: " << mB << "\t\tmG: " << mG << "\t\tmR: " << mR << '\n';
	
	// declare ratioB, ratioG, ratioR and calculate their value
	double ratioB, ratioG, ratioR;
	ratioB = (mB + mG + mR) / (3.0 * mB);
	ratioG = (mB + mG + mR) / (3.0 * mG);
	ratioR = (mB + mG + mR) / (3.0 * mR);
	
	// output information
	cout << "ratioB: " << ratioB << "\t\tratioG: " << ratioG << "\t\tratioR: " << ratioR << '\n';
	
	// multiply every pixel in the three channel with their corresponding ratio (Grey World method)
	for(int i=0;i<_height * _width * _channel;i++)
	{
		if(i%3 == 0)
			*(uint8_t*) &_image[i] = *(uint8_t*) &_image[i] * ratioB >= 255 ? 255 : *(uint8_t*) &_image[i] * ratioB;
		else if(i%3 == 1)
			*(uint8_t*) &_image[i] = *(uint8_t*) &_image[i] * ratioG >= 255 ? 255 : *(uint8_t*) &_image[i] * ratioG;
		else
			*(uint8_t*) &_image[i] = *(uint8_t*) &_image[i] * ratioR >= 255 ? 255 : *(uint8_t*) &_image[i] * ratioR;
	}	
	cout << "---------------------------------\n";
	
	// update Y Cr Cb info
	RGB2YCrCb();
}

void BMP::saturationEnhance(double gamma)
{
	// declare H, S, I
	double* H = new double[_height*_width];
	double* S = new double[_height*_width];
	double* I = new double[_height*_width];
	
	// transform sRGB to HSI 
	RGB2HSI(H,S,I);
	
	// enhance Saturation by power-law transformation
	for(int i=0;i<_height*_width;i++)
		S[i] = min(1.0, pow(S[i] , gamma));
	
	// update image
	HSI2RGB(H,S,I);
	RGB2YCrCb(); 
}

