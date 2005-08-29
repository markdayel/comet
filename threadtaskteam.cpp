
#include "threadtaskteam.h"

#include "stdafx.h"
#include <assert.h>
#include <cmath>
#include <iostream>

taskteam::taskteam()
{
    int err;
    
    // set predicates
    err = pthread_mutex_init(&mutex, NULL);
    if(err != 0){
		std::cerr << "mutex initialisation error" << std::endl;
    }
    err = pthread_cond_init(&start, NULL);
    if(err != 0){
		std::cerr << "cond initialisation error" << std::endl;
    }
    err = pthread_cond_init(&finish, NULL);
    if(err != 0){
		std::cerr << "cond initialisation error" << std::endl;
    }
    ready = false;
    taskfcn = NULL;
}

void taskteam::lockteam()
{
    int err;
    err = pthread_mutex_lock(&mutex);
    // FIXME: Throw an exception here
    if(err != 0){
		std::cerr << "lock team mutex error" << std::endl;
    }
}

void taskteam::unlockteam()
{
    int err;
    err = pthread_mutex_unlock(&mutex);
    // FIXME: Throw an exception here
    if(err != 0){
		std::cerr << "unlock team mutex error" << std::endl;
    }
}

void taskteam::set_taskfcn(void* (*ataskfcn)(void *args, pthread_mutex_t *mutex))
{
	assert(ataskfcn!=NULL);
	taskfcn = ataskfcn;
}

void taskteam::add_task(void* arg)
{
    lockteam();
    assert(arg!=NULL);
    tasks.push(arg); // add the relevant args
    unlockteam();
}

void taskteam::clear_tasks()
{
    lockteam();
    while(!tasks.empty()) {
    	tasks.pop();
    }
    unlockteam();
}

void taskteam::create_team(const int num_workerthreads)
{
    lockteam();

    workers.resize(num_workerthreads);
    
    for(int i=0; i<num_workerthreads; ++i){
    	pthread_t thread;
		int err = pthread_create(&thread, NULL,
				 				 taskteam::wrapper_do_thread_work,
				 				 this);
		if(err != 0){
	    	std::cerr << "Error creating pthread" << std::endl;
		}
		workers[i] = thread;
    }

    unlockteam();
}

int taskteam::do_work()
{
    // ! ! ! ! wait if there are threads working ??
    lockteam();
	assert(taskfcn != NULL);
	
    //std::cout << "telling threads to start" << std::endl;
    ntasks_todo = tasks.size();;
    ready = true;
    int err = pthread_cond_broadcast(&start);
    if(err != 0){
		std::cerr << "error starting the team" << std::endl;
    }
    
    //std::cout << "waiting for finish" << std::endl;
    while( !tasks.empty() ) {
		err = pthread_cond_wait(&finish, &mutex);
		if(err != 0){
	    	std::cerr << "error waiting for work" << std::endl;
		}
    }
    
    //std::cout << "-- Finished" << std::endl;
    ready = false;
    unlockteam();

    assert( tasks.empty() );
    assert( ntasks_todo == 0 );
    
    return 0;
}

void* taskteam::wrapper_do_thread_work(void *my_this)
{
    taskteam* the_team = (taskteam*) my_this;    
    the_team->do_thread_work(NULL);
    
    return NULL;
}

void* taskteam::do_thread_work(void* arg)
{
    int err;
    
    while(GRASS_IS_GREEN) 
    {
	
		lockteam();	
		// wait until there is some work to do
		while( !ready || tasks.empty() ) 
        {
	    	// cout << this << " waiting for start " << threads_ready << std::endl;
	    	err = pthread_cond_wait(&start, &mutex);
	    	if(err != 0)
            {
				std::cerr << "start error" << std::endl;
	    	}
		}
	
		// get the next work unit
		void * task_arg = tasks.front();
		tasks.pop();
		unlockteam();
		
		// Call the fcn to perform the task
		taskfcn(task_arg, &mutex);
		
		lockteam();
		ntasks_todo--;
		// check if all the work is done
		if( ntasks_todo == 0 ) 
        {
	    	// indicate that we are done
	    	err = pthread_cond_broadcast(&finish);
	    	if(err != 0)
            {
				std::cerr << "finish error" << std::endl;
	    	}
		}

		unlockteam();	
    }
    
    return NULL;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// delay function
/*
double f(double x) {
    // do some meaningless work
    for (int i=0;i<5000;++i)
	x = pow(cos(x/(i*2.0)),3);
    return x;
}

// Task A
typedef struct task_a_arg_t {
	int thread_id;
	int bgn_indx;
	int end_indx;
	std::vector<int> *array;
} task_a_arg;


void* task_fcn_a(void* arg)
{
    // cast arg
    task_a_arg *data;
    data = (task_a_arg*) arg;

    // - - - DO WORK
    // cout << "Work A] Thread:" << data->thread_id << " started "
    // 	    <<"(begin, end) :("
    //      << data->bgn_indx << ","
    //      << data->end_indx << ")"
    //      << std::endl;
    
    // ignore this, just simulate some work.
    // loop with some make-work
    for(int i=data->bgn_indx; i<data->end_indx; ++i){
		(*data->array)[i] += 5;
		f(3);
    }
    
    return NULL;
}

// Task B
typedef struct task_b_arg_t {
	int thread_id;
	int bgn_indx;
	int end_indx;
	std::vector<int> *array;
} task_b_arg;

void* task_fcn_b(void* arg)
{
    // cast arg
    task_a_arg *data;
    data = (task_a_arg*) arg;
    
    // - - - DO WORK
    // cout << "Work B] Thread:" << data->thread_id << " started "
    //	 <<"(begin, end) :("
    //   << data->bgn_indx << ","
    //   << data->end_indx << ")"
    //   << std::endl;
    
    // ignore this, just simulate some work.
    // loop with some make-work
    for(int i=data->bgn_indx; i<data->end_indx; ++i){
	(*data->array)[i] += 7;
	f(3);
    }
    
    return NULL;
}
 
// - MAIN - - - - - - - - - - - - 
void sout_array(std::vector<int> & array){
    // write output to show some work was done
    int nout = 5;
    for(int i=0; i<nout; ++i){
	
	int tst_indx = (int)(array.size() / nout)*i + 2;
	std::cout << "array (" << tst_indx << ") "
	     << array[tst_indx] << std::endl;
    }    
}

int main(int argc, char* argv[])
{	
    if (argc < 2) {
		std::cout << "Usage:" << argv[0] << " n_threads" << std::endl;
		return 1;
    }
    int num_workers = atoi(argv[1]);
    
    taskteam tteam_a;
    tteam_a.set_taskfcn(&task_fcn_a);
    
    const int array_size = 1500;
    std::vector<int> array(array_size);
    for(int i=0; i<array_size; ++i){
		array[i] = i;
    }
    
    std::cout << "Thread test using " << num_workers 
    		  << " thread(s)." << std::endl;
    int workunit = array_size / num_workers;
    
    // init the thread args
    task_a_arg targs_a[num_workers];
    for(int i=0; i<num_workers; ++i){
		targs_a[i].bgn_indx = i * workunit;
		targs_a[i].end_indx = (i+1) * workunit; // NB < end_indx
	
		targs_a[i].thread_id = (i);
		targs_a[i].array = &array;
	
		if( i==num_workers-1 ) {
	    	targs_a[i].end_indx = array.size();
		}			
    }
    
    std::cout << "Creating threads." << std::endl;
    
    tteam_a.create_team(num_workers);
    
    for (int i=0; i<num_workers; ++i){
		tteam_a.add_task(&targs_a[i]);
    }    
    tteam_a.do_work();
    std::cout << "a) done." << std::endl;
    sout_array(array);

    // do B
    taskteam tteam_b;
    tteam_b.set_taskfcn(task_fcn_b);
    std::cout << "Creating B threads." << std::endl;
    tteam_b.create_team(num_workers);
    
    // init the thread args
    task_b_arg targs_b[num_workers];
    for(int i=0; i<num_workers; ++i){
		targs_b[i].bgn_indx = i * workunit;
		targs_b[i].end_indx = (i+1) * workunit; // NB < end_indx
	
		targs_b[i].thread_id = (i);
		targs_b[i].array = &array;
	
		if( i==num_workers-1 ) {
	    	targs_b[i].end_indx = array.size();
		}			
		
		tteam_b.add_task(&targs_b[i]);
    }
        
    tteam_b.do_work();
    std::cout << "b) done." << std::endl;
    sout_array(array);
    
    std::cout << "All done." << std::endl;
    return 0;
}
*/
