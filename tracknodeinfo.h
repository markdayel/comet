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


#ifndef tracknodeinfo_H
#define tracknodeinfo_H


class tracknodeinfo
{
public:
    tracknodeinfo(){};
    ~tracknodeinfo(){};
    //tracknodeinfo(int const & setnodenum, const vect & setposn, const vect & setnucposn, const rotationmatrix & setrotation, int const & setframe ):
    tracknodeinfo(int const & setnodenum, const double &setx, const double &sety, const double &setz,
        const int & setBMPx, const int & setBMPy, int const & setframe ):
    nodenum(setnodenum),frame(setframe),x(setx), y(sety), z(setz), BMPx(setBMPx), BMPy(setBMPy)
    {
        //posn = setposn;
        //nucposn = setnucposn;
        //rotation = setrotation;
    }  

    int nodenum;
    //vect posn, nucposn;
    //rotationmatrix rotation;
    int frame;
    double x,y,z;
    int BMPx,BMPy;

    friend ostream &operator<<(ostream &stm, tracknodeinfo const &trk){
	stm << trk.nodenum << " " << trk.frame << " "
        << trk.x << " " << trk.y << " " << trk.z << " "
        << trk.BMPx << " " << trk.BMPy;
	return stm;
    }
    
    friend istream &operator>>(istream &stm, tracknodeinfo &trk){
	stm >> setprecision(SAVE_DATA_PRECISION) 
        >> trk.nodenum >> trk.frame 
        >> trk.x >> trk.y >> trk.z
        >> trk.BMPx >> trk.BMPy;
	return stm;
    }
    
};

#endif

