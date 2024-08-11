#include "BMP.h"

using namespace std;

int main()
{
	BMP bmp;
	
	// input1.bmp
	bmp.read("input1.bmp");
	bmp.show();	
	bmp.quan(2);
	bmp.write("output1_1.bmp");
	bmp.quan(4);
	bmp.write("output1_2.bmp");
	bmp.quan(6);
	bmp.write("output1_3.bmp");
	
	// input2.bmp
	bmp.read("input2.bmp");
	bmp.show();	
	bmp.quan(2);
	bmp.write("output2_1.bmp");
	bmp.quan(4);
	bmp.write("output2_2.bmp");
	bmp.quan(6);
	bmp.write("output2_3.bmp");
	return 0;
}
