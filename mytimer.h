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

