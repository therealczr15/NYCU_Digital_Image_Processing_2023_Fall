#include "BMP.h"

using namespace std;

int main()
{
	BMP bmp;

	// input1.bmp
	bmp.read("output1_1.bmp");
	bmp.histEqual();
	bmp.saturationEnhance(0.6);
	bmp.gammaCorrection(0.8);
	bmp.write("output1_2.bmp");
	
	// input2.bmp
	bmp.read("output2_1.bmp");
	bmp.histEqual();
	bmp.gammaCorrection(0.95);
	bmp.write("output2_2.bmp");
	
	// input3.bmp
	bmp.read("output3_1.bmp");
	bmp.saturationEnhance(0.8);
	bmp.gammaCorrection(0.9);
	bmp.write("output3_2.bmp");
	
	// input4.bmp
	bmp.read("output4_1.bmp");
	bmp.histEqual();
	bmp.gammaCorrection(1.05);
	bmp.saturationEnhance(0.8);
	bmp.write("output4_2.bmp");
	return 0;
}
