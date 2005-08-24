// thread_test.cpp
//
#ifndef THREADTASKTEAM_H_
#define THREADTASKTEAM_H_

#include "pthread.h"
#include <vector>
#include <queue>

class taskteam 
{
protected:
    std::vector<pthread_t> workers;
    std::queue<void*> tasks;

    pthread_mutex_t mutex;
    pthread_cond_t start;
    pthread_cond_t finish;
    bool ready;
    size_t ntasks_todo;
    void* (*taskfcn)(void *args, pthread_mutex_t *mutex);
      
    void lockteam();
    void unlockteam();
    static void* wrapper_do_thread_work(void* wrapper_arg);
    void* do_thread_work(void* thread_arg);
    
public:
    taskteam();
    virtual ~taskteam() {};
    void  set_taskfcn(void* (*ataskfcn)(void *args, pthread_mutex_t *mutex));
    void  create_team(const int num_workerthreads);
    void  add_task(void* arg);
    void  clear_tasks();
    int   do_work();
};

#endif // THREADTASKTEAM_H_
