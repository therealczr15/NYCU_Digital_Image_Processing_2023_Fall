#include "BMP.h"

using namespace std;

int main()
{
	BMP bmp;
	
	// input3.bmp
	bmp.read("input3.bmp");
	bmp.show();	
	
	// iteration 10 times
	for(int i=0;i<10;i++)
		bmp.denoise();
	bmp.write("output3_1.bmp");
	
	// input3.bmp
	bmp.read("input3.bmp");
	bmp.show();	
	
	// iteration 50 times
	for(int i=0;i<50;i++)
		bmp.denoise();
	bmp.write("output3_2.bmp");
	return 0;
}
