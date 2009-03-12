/*
    comet - an actin-based motility simulator
    Copyright (C) 2005 Mark J Dayel

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>
*/

#include "stdafx.h"

#include "vect.h"
// DUMMY Header file for systems not linking to vtk      

class CometVtkVis {
 public:
  CometVtkVis(bool,bool) {};
  ~CometVtkVis() {};
  
  void buildVTK(int, vect&, vect&);
};

void CometVtkVis::buildVTK(int, vect&, vect&) 
{
  cout << "(vtk not linked)" << endl;
};

