#include "BMP.h"

using namespace std;

int main()
{
	BMP bmp;
	
	// input1.bmp
	bmp.read("input1.bmp");
	bmp.grayWorld();
	bmp.write("output1_1.bmp");

	// input2.bmp
	bmp.read("input2.bmp");
	bmp.maxRGB();
	bmp.write("output2_1.bmp");

	// input3.bmp
	bmp.read("input3.bmp");
	bmp.maxRGB();
	bmp.write("output3_1.bmp");
	
	// input4.bmp
	bmp.read("input4.bmp");
	bmp.grayWorld();
	bmp.write("output4_1.bmp");
	
	return 0;
}
