/*
    comet - a actin-based motility simulator
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

#ifndef mytimer_H
#define mytimer_H

#include <ostream>




class mytimer
{
public:
    mytimer(void)
    {
        setprecision(2);
        start();
    };
    ~mytimer(void){};
    double starttime;
    unsigned int precision;

    void start()
    {
        starttime = getsectime();
    }
    
    void setprecision(unsigned int p)
    {  
        precision = p;

#ifndef TIMERCPUTIME  
        precision = 0;  // force 0 precision for walltime since it's in seconds
#endif

    }

#ifndef TIMERCPUTIME
    double getsectime() const
    {   // return current time in seconds
        return (double) time(NULL) ; // wall time
    }
#else
    double getsectime() const
    {   // return current time in seconds
        return (double) clock()  / (double) CLOCKS_PER_SEC; // processor time
    }

#endif
    
};

ostream& operator << (ostream &out, const mytimer& timer)
{
    out << setprecision(timer.precision) <<  timer.getsectime() - timer.starttime ;
    return out;
}


#endif

