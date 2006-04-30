/*
  Copyright (C) 2005 Mark J Dayel

  You may not distribute, deploy, provide copies of, rent, loan, lease, 
  transfer, or grant any rights in the software or derivative works thereof 
  in any form to any person.  Reproduction, adaptation, or translation of 
  this program without prior written permission from the author is 
  prohibited.  Proofs of infringement of copyright violations include but 
  not limited to similar code style and structures, similar layout and 
  design, similar algorithm design, and containing parts of the original 
  software source code.  Copyright notice must remain intact and cannot be 
  removed without prior written permission from the author.
*/

#include "stdafx.h"

// DUMMY Header file for systems not linking to vtk

class CometVtkVis {
 public:
  CometVtkVis() {};
  ~CometVtkVis() {};
  
  void buildVTK(int);
};

void CometVtkVis::buildVTK(int) 
{
  cout << "(vtk not linked)" << endl;
};

