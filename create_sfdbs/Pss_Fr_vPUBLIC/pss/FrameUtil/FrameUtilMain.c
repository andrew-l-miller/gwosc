#include "frame_util.h"

//int main()
int main(int argc, char *argv[])
// argc[0] batch file
// if arg is present, batch operation
{
	if (argc > 1) frame_util_bat(argv[1]);

	else frame_util();

	return 0;
//	frame_util();
}
