#include "BMP.h"

using namespace std;

int main()
{
	BMP bmp;
	
	// input1.bmp
	bmp.read("input1.bmp");	
	bmp.show();	
	bmp.flip();
	bmp.write("output1_flip.bmp");
	
	// input2.bmp
	bmp.read("input2.bmp");	
	bmp.show();	
	bmp.flip();
	bmp.write("output2_flip.bmp");
	return 0;
}
