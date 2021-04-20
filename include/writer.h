#include "parser.h"
#ifndef WRITER_H
#define WRITER_H

// -----------------------------------------|
// Functions to write the results in Jmol   |
// -----------------------------------------|

const char* FCT_atm_name(double norm);
void write(Mesh &mesh,darray &disp,dvector &VM);
void writeJmol(Mesh &mesh,darray &disp,dvector &VM);

#endif