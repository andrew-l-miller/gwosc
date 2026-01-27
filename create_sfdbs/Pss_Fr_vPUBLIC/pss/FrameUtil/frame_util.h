#include <stdio.h>

void frame_util();

int frame_util_bat(char *filbat);

void print_dump_std(char *dumpfile, char *filename);

void print_dump(char *dumpfile, char *filename);

void to_sds(char *filename, char *direcout, char *channel, char *antacr, char *datatype, char *sub_typ);

FILE *new_sdsfile(double t0, double t0gps, long ULeapS, char *direcout, char *antacr, char *datatype, char *channel, char *sub_typ, double dt);
