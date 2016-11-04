/**
 * This class defines a FIFO cache (First In First Out) scheduler
 **/
#ifndef FIFO_CACHE_NC_SCHEDULER_HPP
#define FIFO_CACHE_NC_SCHEDULER_HPP

#include <queue>
#include <map>
#include <cmath>
#include <cassert>

#include <graphlab/logger/logger.hpp>
#include <graphlab/graph/graph.hpp>
#include <graphlab/scope/iscope.hpp>
#include <graphlab/tasks/update_task.hpp>
#include <graphlab/schedulers/ischeduler.hpp>
#include <graphlab/parallel/pthread_tools.hpp>
#include <graphlab/schedulers/support/vertex_task_set.hpp>
#include <graphlab/schedulers/support/direct_callback.hpp>
//#include <util/shared_termination.hpp>
#include <graphlab/util/task_count_termination.hpp>

// #include <bitmagic/bm.h>


#include <graphlab/macros_def.hpp>

namespace graphlab {


 
  template<typename Graph>
  class fifo_cache_nc_scheduler: public ischeduler<Graph> {
  public:
    typedef Graph graph_type;
    typedef ischeduler<Graph> base;

    typedef typename base::iengine_type iengine_type;
    typedef typename base::update_task_type update_task_type;
    typedef typename base::update_function_type update_function_type;
    typedef typename base::callback_type callback_type;
    typedef typename base::monitor_type monitor_type;
     
    
    
  private:
    using base::monitor;

  public:

    fifo_cache_nc_scheduler(iengine_type* engine,
                   Graph& g, 
                   size_t ncpus)  : 
      callbacks(ncpus, direct_callback<Graph>(this, engine)), 
      vertex_tasks(g.num_vertices()), g(g) {
      numvertices = g.num_vertices();
      decisions = 0;
      last_scheduled = new vertex_id_t[ncpus];
      for (size_t id = 0; id < ncpus; id++)
	last_scheduled[id] = 0xffffffff;
    }

    ~fifo_cache_nc_scheduler() {
      delete last_scheduled;
    }

    callback_type& get_callback(size_t cpuid) {
      return callbacks[cpuid];
    }

    /** Get the next element in the queue */
    sched_status::status_enum get_next_task(size_t cpuid, update_task_type &ret_task) {    
      if (terminator.finish()) {
        if (cpuid == 0)
                logger(LOG_INFO, "#total decisions. %d\n", decisions);
	return sched_status::COMPLETE;
      }
      bool success(false);
      queue_lock.lock();
      if(!task_set.empty()) {
	      int max_cached = 0;
	      vertex_id_t cur = -1;
	      if (last_scheduled[cpuid] != 0xffffffff) {
		      vertex_id_t last = last_scheduled[cpuid];
		      if (task_set.find(last) != task_set.end()) {
			      g.query_cache(&g.vertex_data(last), &max_cached);
			      cur = last;
		      }
		      /*
			 const std::vector<edge_id_t> & in_edges = g.in_edge_ids(last);
			 std::vector<edge_id_t>::const_iterator itr = in_edges.begin();
			 for (; itr != in_edges.end(); itr++) {
			 edge_id_t edge = *itr;
			 vertex_id_t v = g.source(edge);
			 if (task_set.find(v) == task_set.end())
			 continue;
			 int cached = 0;
			 g.query_cache(&g.vertex_data(v), &cached);
			 if (cached > max_cached) {
			 max_cached = cached;
			 cur = v;
			 }
			 }	
			 const std::vector<edge_id_t> & out_edges = g.out_edge_ids(last);
			 itr = out_edges.begin();
			 for (;itr != out_edges.end(); itr++) {
			 edge_id_t edge = *itr;
			 vertex_id_t v = g.target(edge);
			 if (task_set.find(v) == task_set.end())
			 continue;
			 int cached = 0;
			 g.query_cache(&g.vertex_data(v), &cached);
			 if (cached > max_cached) {
			 max_cached = cached;
			 cur = v;
			 }
			 }
		       */
	      }

	      if (cur == -1) {
		      ret_task = task_set.begin()->second;
		      task_set.erase(task_set.begin());
	      } else {
		      ret_task = task_set[cur];
		      task_set.erase(cur);
	      }

	      last_scheduled[cpuid] = ret_task.vertex();
	      decisions++;
	      success = true;
      }
      queue_lock.unlock();
      
      if(success) {
        if (monitor != NULL) {
          double priority = vertex_tasks.top_priority(ret_task.vertex());
          monitor->scheduler_task_scheduled(ret_task, priority);
        }
        vertex_tasks.remove(ret_task);
        return sched_status::NEWTASK;
      } else {
        return sched_status::WAITING;
      }
    } // end of get_next_task
    

    void add_task(update_task_type task, double priority) {
      if (vertex_tasks.add(task)) {
        terminator.new_job();
        queue_lock.lock();
	task_set[task.vertex()] = task;
        queue_lock.unlock();
        if (monitor != NULL) 
          monitor->scheduler_task_added(task, priority);
      } else {
        if (monitor != NULL) 
          monitor->scheduler_task_pruned(task);
      }
    } // end of add_task

    void add_tasks(const std::vector<vertex_id_t> &vertices,
                   update_function_type func,
                   double priority) {
      foreach(vertex_id_t vertex, vertices) {
        add_task(update_task_type(vertex, func), priority);
      }
    } // end of add_tasks

    void add_task_to_all(update_function_type func, double priority) {
      for (vertex_id_t vertex = 0; vertex < numvertices; ++vertex){
        add_task(update_task_type(vertex, func), priority);
      }
    } // end of add_task_to_all

  
    void update_state(size_t cpuid,
                      const std::vector<vertex_id_t> &updated_vertices,
                      const std::vector<edge_id_t>& updatededges) {};

    void scoped_modifications(size_t cpuid, vertex_id_t rootvertex,
                              const std::vector<edge_id_t>& updatededges){}

    void completed_task(size_t cpuid, const update_task_type &task) {
      terminator.completed_job();
    }

    void abort() { terminator.abort(); }
 
    void restart() { terminator.restart(); }

  private:
    size_t numvertices; /// Remember the number of vertices in the graph
  
    std::map<vertex_id_t, update_task_type> task_set; /// The actual task queue
    
    spinlock queue_lock; // The lock to get an element from the queue

    /// The callbacks pre-created for each cpuid
    std::vector< direct_callback<Graph> > callbacks; 

    // Task set for task pruning
    vertex_task_set<Graph> vertex_tasks;
  
    task_count_termination terminator;

    Graph& g;

    vertex_id_t* last_scheduled;

    long decisions;
  }; 
} // end of namespace graphlab
#include <graphlab/macros_undef.hpp>

#endif
