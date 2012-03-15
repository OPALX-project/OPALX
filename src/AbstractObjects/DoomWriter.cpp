// ------------------------------------------------------------------------
// $RCSfile: DoomWriter.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: DoomWriter
//   This class allows construction of an object to be written to the DOOM
//   data base and writes it out automatically at destruction time.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/12/01 10:37:09 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/DoomWriter.h"
#include "doom.h"
#include "doomex.h"
#include <algorithm>
#include <string.h>


// Class DoomWriter
//   NOTE: This class must use malloc()/realloc(), since the objects
//         allocated will be deleted by the DOOM data base using free().
// ------------------------------------------------------------------------

DoomWriter::DoomWriter()
{
  doomStruct = make_obj("", 0, 0, 0, 0);
}

DoomWriter::DoomWriter(const string &name)
{
  doomStruct = make_obj("", 0, 0, 0, 0);
  truncateName(name, doomStruct->key);
}


DoomWriter::DoomWriter(const string &name, const int keyList[])
{
  doomStruct = make_obj("", 0, 0, 0, 0);
  make_key(name.c_str(), keyList, doomStruct->key);
}


DoomWriter::~DoomWriter()
{
  // Avoid saving an unnamed object.
  if (doomStruct->key[0] != 0) {
    // Determine the length for the a_char array.
    int c_char = 0;
    int c_obj = doomStruct->c_obj;


    // Tue Apr 18 ada: rules out the case
    //                 in which doomStruct->key[0] != 0
    // but the string is NULL and strlen crashes on the
    // linux systems.
    // May the problem is deeper into OPAL9 and the DOOM stuff?

    for (int i = 0; i < c_obj; ++i) {
      if (doomStruct->names[i] != 0 )
	c_char += strlen(doomStruct->names[i]) + 1;
      else
	c_char++;
    }

    // Build the a_char array.
    if (c_char > 0) {
      int c = 0;
      doomStruct->a_char = (char *) malloc(c_char * sizeof(char));
      for (int i = 0; i < c_obj; ++i) {
	if (doomStruct->names[i] != 0) {
	  strcpy(doomStruct->a_char + c, doomStruct->names[i]);
	  c += strlen(doomStruct->names[i]);
	  free(doomStruct->names[i]);
	  doomStruct->names[i] = 0;
	}
	
	doomStruct->a_char[c++] = '|';
      }
      doomStruct->l_char = doomStruct->c_char = c_char;
    }

    // Build the a_obj array.
    doomStruct->p_obj =
      (struct object **) calloc(sizeof(struct object *), c_obj);

    // Write out the DOOM data base object.
    // Do not destroy, the object will be owned by the data base.
    doom_save(doomStruct);
  }
}


void DoomWriter::putInt(int index, const int value)
{
  int l_int = doomStruct->l_int;

  if (index >= l_int) {
    int nl_int = 2;
    do {
      nl_int = nl_int + nl_int;
    } while (index >= nl_int);

    int *a_int = doomStruct->a_int;
    doomStruct->a_int = (int *) calloc(sizeof(int), nl_int);
    doomStruct->l_int = nl_int;

    if (a_int != 0) {
      std::copy(a_int, a_int + l_int, doomStruct->a_int);
      free(a_int);
    }
  }

  doomStruct->a_int[index] = value;
  if (index >= doomStruct->c_int) doomStruct->c_int = index + 1;
}


void DoomWriter::putReal(int index, double value)
{
  int l_dble = doomStruct->l_dble;

  if (index >= l_dble) {
    int nl_dble = 2;
    do {
      nl_dble = nl_dble + nl_dble;
    } while (index >= nl_dble);

    double *a_dble = doomStruct->a_dble;
    doomStruct->a_dble = (double *) calloc(sizeof(double), nl_dble);
    doomStruct->l_dble = nl_dble;
    if (a_dble != 0) {
      std::copy(a_dble, a_dble + l_dble, doomStruct->a_dble);
      free(a_dble);
    }
  }

  doomStruct->a_dble[index] = value;
  if (index >= doomStruct->c_dble) doomStruct->c_dble = index + 1;
}


void DoomWriter::putString(int index, const string &value)
{
  int l_obj = doomStruct->l_obj;

  if (index >= l_obj) {
    int nl_obj = 2;
    do {
      nl_obj = nl_obj + nl_obj;
    } while (index >= nl_obj);

    char **names = doomStruct->names;
    doomStruct->names = (char **) calloc(sizeof(char *), nl_obj);
    doomStruct->l_obj = nl_obj;
    if (names != 0) {
      std::copy(names, names + l_obj, doomStruct->names);
      free(names);
    }
  }

  if (doomStruct->names[index] != 0)
    free(doomStruct->names[index]);

  doomStruct->names[index] = (char *) malloc(value.length() + 1);
  strcpy(doomStruct->names[index], value.c_str());
  if (index >= doomStruct->c_obj) doomStruct->c_obj = index + 1;
}


void DoomWriter::setObjectName(const string &name)
{
  truncateName(name, doomStruct->key);
}


void DoomWriter::setObjectName(const string &name, const int keyList[])
{
  make_key(name.c_str(), keyList, doomStruct->key);
}


void DoomWriter::setBaseName(const string &name)
{
  truncateName(name, doomStruct->base_name);
}


void DoomWriter::setParentName(const string &name)
{
  truncateName(name, doomStruct->par_name);
}


void DoomWriter::setTypeName(const string &name)
{
  truncateName(name, doomStruct->obj_type);
}


void DoomWriter::truncateName(const string &name, char nam[24])
{
  int len = name.length();
  if (len >= 24) len = 23;
  strncpy(nam, name.data(), len);
  nam[len] = '\0';
}
