#include "BMP.h"

using namespace std;

int main()
{
	BMP bmp;
	
	// input1.bmp
	bmp.read("input1.bmp");
	bmp.show();
	
	// brute force method	
	//bmp.bruteForce(20);
	
	// gamma correction method, gamma = 0.75
	bmp.gammaCorrection(0.75);
	bmp.write("output1_1.bmp");
	
	// input1.bmp
	bmp.read("input1.bmp");
	bmp.show();	
	
	// brute force method
	//bmp.bruteForce(40);
	
	// gamma correction method, gamma = 0.5
	bmp.gammaCorrection(0.5);
	bmp.write("output1_2.bmp");
	return 0;
}
