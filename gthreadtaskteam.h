// thread_test.cpp
//
#ifndef THREADTASKTEAM_H_
#define THREADTASKTEAM_H_

#ifdef _WIN32
#include <wherever\timeb.h\is\>
#else
#include <sys/timeb.h>
#endif

#include "pthread.h"
#include <vector>
#include <queue>
#include <string>
#include <map>

class taskteam 
{
private:
    std::queue<void*> tasks;
    std::vector<pthread_t> workers;
    pthread_mutex_t mutex;
    pthread_cond_t start;
    pthread_cond_t finish;
    bool ready;
    size_t ntasks_todo;
    void* (*taskfcn)(void *args, pthread_mutex_t *mutex);
    
    // timing 
    std::string team_id;
    // relies on p_thread being int, which it might not be.    
    std::map<pthread_t, double> cumulativetimes;
    std::map<pthread_t, int> cumulativeruns;
    timeb starttime;
    timeb endtime;
    int ntaskteam_runs;
    
    void lockteam();
    void unlockteam();
    static void* wrapper_do_thread_work(void* wrapper_arg);
    void* do_thread_work(void* thread_arg);
    double elapsedtime(timeb start, timeb end);
    
public:
    taskteam();
    ~taskteam();
    void  set_taskfcn(void* (*ataskfcn)(void *args, pthread_mutex_t *mutex));
    void  create_team(const int num_workerthreads);
    void  add_task(void* arg);
    void  clear_tasks();
    int   do_work();
    // timing
    void showtimings();
    void set_teamid();
    std::string get_teamid();
};

#endif // THREADTASKTEAM_H_
