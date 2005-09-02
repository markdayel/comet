
#include "threadedtaskqueue.h"

#include <assert.h>
#include <iostream>

// static members are not *guaranteed* to have C linkage
extern "C" void* wrapper_activate_workers(void *my_this)
{
  TaskQueue* the_team = (TaskQueue*) my_this;    
  return the_team->activate_workers();
}

TaskQueue::TaskQueue()
{
  int err;
  
  // set predicates
  err = pthread_mutex_init(&mutex, 0);
  if(err != 0){
    std::cerr << "mutex initialisation error" << std::endl;
  }
  err = pthread_cond_init(&start, 0);
  if(err != 0){
    std::cerr << "cond initialisation error" << std::endl;
  }
  err = pthread_cond_init(&finish, 0);
  if(err != 0){
    std::cerr << "cond initialisation error" << std::endl;
  }
  ready = false;
}

void TaskQueue::lockteam()
{
  int err;
  err = pthread_mutex_lock(&mutex);
  // FIXME: Throw an exception here
  if(err != 0){
    std::cerr << "lock team mutex error" << std::endl;
  }
}

void TaskQueue::unlockteam()
{
  int err;
  err = pthread_mutex_unlock(&mutex);
  // FIXME: Throw an exception here
  if(err != 0){
    std::cerr << "unlock team mutex error" << std::endl;
  }
}

void TaskQueue::queue_task(TaskFcn_ptr task_fcn, void* task_args)
{
  assert(task_fcn !=0);
  assert(task_args!=0);
  
  TaskDetails task(task_fcn, task_args);
  
  lockteam();
  tasks.push(task); // add the function and args for this task
  unlockteam();
}

void TaskQueue::clear_tasks()
{
  lockteam();
  // REVISIT: what to do if tasks are running?
  while(!tasks.empty()) {
    tasks.pop();
  }
  unlockteam();
}

void TaskQueue::create_threads(const int num_workerthreads)
{
  lockteam();
  
  workers.resize(num_workerthreads);

  for(int i=0; i<num_workerthreads; ++i){
    pthread_t thread;
    int err = pthread_create(&thread, 0,
			     wrapper_activate_workers,
			     this);
    if(err != 0){
      std::cerr << "Error creating pthread" << std::endl;
    }
    workers[i] = thread;
  }
  
  unlockteam();
}

void TaskQueue::start_tasks()
{
  // ! ! ! ! REVISIT: add a wait if there are threads working ??  
  assert(ready==false);

  lockteam();
  
  int err = pthread_cond_broadcast(&start);
  if(err != 0){
    std::cerr << "error starting the team" << std::endl;
  }
  ntasks_todo = tasks.size();;
  ready = true;
  
  unlockteam();
}

void TaskQueue::wait_to_finish_tasks()
{
  int err;
  lockteam();
  
  while( ntasks_todo != 0  ) {
    err = pthread_cond_wait(&finish, &mutex);
    if(err != 0){
      std::cerr << "error waiting for work" << std::endl;
    }
  }    
  ready = false;
  
  assert( tasks.empty() );
  assert( ntasks_todo == 0 );
  
  unlockteam();
}

void TaskQueue::complete_current_tasks()
{
  start_tasks();
  wait_to_finish_tasks();
}

// Treat as a private/protected member fcn
void* TaskQueue::activate_workers()
{
  int err;
  
  while(true) { // FIXME allow thread termination?
    
    lockteam();		
    // wait until there is some work to do
    while( !ready || tasks.empty() ) {
      err = pthread_cond_wait(&start, &mutex);
      if(err != 0){
	std::cerr << "start error" << std::endl;
      }
    }
    
    // get the next work unit
    TaskDetails task = tasks.front();
    tasks.pop();
    unlockteam();
    
    // Call the fcn to perform the task
    task.fcn(task.args, &mutex);

    lockteam();
    ntasks_todo--;
    // check if all the work is done
    // not same as tasks.empty(), as tasks can be in process
    // of executing.
    if( ntasks_todo == 0 ) {
      // indicate that we are done
      err = pthread_cond_broadcast(&finish);
      if(err != 0){
	std::cerr << "finish error" << std::endl;
      }
    }
    unlockteam();
  }
  // UNREACHABLE ...
  return 0;
}

