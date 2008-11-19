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

// thread_test.cpp
//
#ifndef THREADEDTASKQUEUE_H_
#define THREADEDTASKQUEUE_H_

#include "pthread.h"
#include <vector>
#include <queue>


class TaskQueue // sorry couldn't look at lower case user defined types any longer
{
protected:
    typedef void* (*TaskFcn_ptr)(void *args);//, pthread_mutex_t *mutex);
    struct TaskDetails {
      TaskFcn_ptr fcn;
      void* args;
      
      TaskDetails(const TaskFcn_ptr &the_fcn, void* the_args)
	: fcn(the_fcn), args(the_args) {}
    };
    std::vector<pthread_t>  workers;
    std::queue<TaskDetails> tasks;
    
    pthread_mutex_t mutex;
    pthread_mutex_t queue_mtx;
    pthread_cond_t start;
    pthread_cond_t finish;
    bool ready, finished_tasks;
    int ntasks_todo;
    size_t ntasks_inprogress;
    
    void lockteam();
    void unlockteam();
    void do_task();
    
 public:
    TaskQueue();
    virtual ~TaskQueue() {};
    virtual void* activate_workers();
    // TODO override copy and assignment cstrs
    virtual void create_threads(const int num_workerthreads);
    void queue_task(TaskFcn_ptr task_fcn, void* task_args);
    void clear_tasks();
    virtual void start_tasks();
    virtual void wait_to_finish_tasks();
    virtual void complete_queued_tasks();
};

#endif // THREADEDTASKQUEUE_H_
