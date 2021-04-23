#include "mesh.h"
#ifndef WRITER_H
#define WRITER_H

// ----------------------------------|
// Functions available in the file   |
// ----------------------------------|

const char* FCT_atm_name(double norm);
void graph(Mesh &mesh,darray &disp,int step);
void write(Mesh &mesh,darray &disp,dvector &VM);
void writeJmol(Mesh &mesh,darray &disp,dvector &VM);

#endif