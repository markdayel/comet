
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

