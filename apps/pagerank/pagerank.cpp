/*
 *  PageRank tutorial application.
 *  pagerank.cpp
 *
 *  For description of the PageRank algorithm, see Wikipedia article
 *  http://en.wikipedia.org/wiki/Pagerank
 */

#include <string>
#include <vector>

#include <graphlab.hpp>
#include <graphlab/macros_def.hpp>

// Constants for the algorithm. Better way would be to
// pass them in the shared_data to the update function, but for
// the simplicity of this example, we simply define them here.

#define termination_bound 1e-5
#define damping_factor 0.85   // PageRank damping factor

/**
 * Edge data represents the weight as well as the weight times the
 * last value of the source vertex when the target value was computed.
 */
struct edge_data {
  float weight;
  float old_source_value;
  edge_data(float weight) :
    weight(weight), old_source_value(0) { } 
    
}; // End of edge data

static int global = 0;
/**
 * Stores the value and the self weight
 */
struct vertex_data {
  float value;
  float self_weight; // GraphLab does not support edges from vertex to itself, so
  				     // we save weight of vertex's self-edge in the vertex data
  std::vector<double>* diff;
  int label;
  vertex_data(float value = 1) : value(value), self_weight(0) {
    label = __sync_fetch_and_add(&global, 1);
    diff = new std::vector<double>();
  }

  vertex_data(const vertex_data& copy) {
    this->value = copy.value;
    this->self_weight = copy.self_weight;
    this->label = copy.label;
    this->diff = new std::vector<double>(*copy.diff);
  }

  ~vertex_data() {
    delete diff;
  }

}; // End of vertex data

static int wasted = 0;
static int total = 0;

//! The type of graph used in this program
typedef graphlab::graph<vertex_data, edge_data> pagerank_graph;


/**
 * The collection of graphlab types restricted to the graph type used
 * in this program.
 */
typedef graphlab::types<pagerank_graph> gl_types;




/**
 * The Page rank update function
 */
void pagerank_update(gl_types::iscope &scope,
                     gl_types::icallback &scheduler,
                     gl_types::ishared_data* shared_data) {
  //printf("vertex %d: update fn starts\n", scope.vertex_data().label); 
  
                       
  // Get the data associated with the vertex
  vertex_data& vdata = scope.vertex_data();
  
  // Sum the incoming weights; start by adding the 
  // contribution from a self-link.
  float sum = vdata.value * vdata.self_weight;
  std::vector<graphlab::edge_id_t> in_edges = scope.in_edge_ids();
  
  int size = scope.in_edge_ids().size();
  int count = 0;
  int i = 0;
  foreach(graphlab::edge_id_t eid, scope.in_edge_ids()) {
    // Get the neighobr vertex value
    const vertex_data& neighbor_vdata =
      scope.const_neighbor_vertex_data(scope.source(eid));
    double neighbor_value = neighbor_vdata.value;
   
    if ((*vdata.diff)[i] < 0)
    	(*vdata.diff)[i++] = neighbor_value;
    else {
	if(std::fabs((*vdata.diff)[i] - neighbor_value) < termination_bound) {
	  count++;
	  i++;
	}
	else
	  (*vdata.diff)[i++] = neighbor_value;
    }

    // Get the edge data for the neighbor
    edge_data& edata = scope.edge_data(eid);
    // Compute the contribution of the neighbor
    double contribution = edata.weight * neighbor_value;
    
    // Add the contribution to the sum
    sum += contribution;
    
    // Remember this value as last read from the neighbor
    edata.old_source_value = neighbor_value;
  }

  //printf("vertex %d wasted %d, total %d, ratio %.2f\n", 
  //              scope.vertex_data().label, count, size, count*1.0f/size);
  __sync_fetch_and_add(&wasted, count);
  __sync_fetch_and_add(&total, size);

  // compute the jumpweight
  sum = (1-damping_factor)/scope.num_vertices() + damping_factor*sum;
  vdata.value = sum;
  
  // Schedule the neighbors as needed
  foreach(graphlab::edge_id_t eid, scope.out_edge_ids()) {
    edge_data& outedgedata = scope.edge_data(eid);
    
    // Compute edge-specific residual by comparing the new value of this
    // vertex to the previous value seen by the neighbor vertex.
    double residual =
      outedgedata.weight *
      std::fabs(outedgedata.old_source_value - vdata.value);
    // If the neighbor changed sufficiently add to scheduler.
    if(residual > termination_bound) {
      gl_types::update_task task(scope.target(eid), pagerank_update);
      scheduler.add_task(task, residual);
    }
  }

  //printf("vertex %d: update fn ends\n", scope.vertex_data().label); 
} // end of pagerank update function





// Creates simple 5x5 graph
void create_graph(pagerank_graph& graph, std::string filename) {
	// if this is a text file
	if (filename.substr(filename.length()-3,3)=="txt") {
		std::ifstream fin(filename.c_str());
		size_t edgecount = 0;
		while(fin.good()) {
			// read a line. if it begins with '#' it is a comment
			std::string line;
			std::getline(fin, line);
			if (line[0] == '#') {
				std::cout << line << std::endl;
				continue;
			}

			// read the vertices
			std::stringstream s(line);
			size_t srcv, destv;
			s >> srcv >> destv;
			
			// make sure we have all the vertices we need
			while (graph.num_vertices() <= srcv || graph.num_vertices() <= destv) {
				vertex_data v;
				v.value = 1.0;
				v.self_weight = 0.0;
				graph.add_vertex(v);
			}
			
			// for graph that has self edges, the weight of self edges are stored in the
			// vertice's self_weight field.
			if (srcv == destv) {
				// selfedge
				vertex_data& v = graph.vertex_data(srcv);
				// we use selfweight temporarily to store the number of self edges
				v.self_weight += 1.0;
			} else {
				// check if the edge already exists
				std::pair<bool, graphlab::edge_id_t> edgefind = graph.find(srcv, destv);
				if (edgefind.first) {
					edge_data& edata = graph.edge_data(edgefind.first);
					// if it does, increment the weight
					edata.weight += 1.0;
				}
				else {
					// we use weight temporarily as a counter
					edge_data e(1.0);
					graph.add_edge(srcv, destv, e);
				}
			}
			++edgecount;
			if (edgecount % 1000000 == 0)
				std::cout << edgecount << " Edges inserted" << std::endl;
      		}
	} else {
		printf("graph types other than .txt not supported right now.\n");
		assert(0);
	}
	
	// now we have to go through the edges again and figure out the
	// weights
	for (graphlab::vertex_id_t i = 0;i < graph.num_vertices(); ++i) {
		vertex_data& v = graph.vertex_data(i);
		// count the effective number of out edges
		double weight = 0;
		foreach(graphlab::edge_id_t e, graph.out_edge_ids(i)){
			edge_data edata = graph.edge_data(e);
			weight += edata.weight;
		}
		// remember to count selfweights
		weight += v.self_weight;
		if (weight != 0) {
			// now the weights should all be scaled to 1
			foreach(graphlab::edge_id_t e, graph.out_edge_ids(i)) {
				edge_data& edata = graph.edge_data(e);
				edata.weight /= weight;
			}
			v.self_weight /= weight;
		}
		//update neighbors counts
		v.diff->reserve(graph.in_edge_ids(i).size());
		for (int j = 0; j < graph.in_edge_ids(i).size(); j++) {
			(*v.diff).push_back(-1.0);
		}
	}
        
	graph.finalize();
}

FILE * fp;
#include <pthread.h>

int main(int argc, char** argv) {
  global_logger().set_log_level(LOG_INFO);
  global_logger().set_log_to_console(true);
  logger(LOG_INFO, "PageRank starting\n");
  
  int numCPU = sysconf(_SC_NPROCESSORS_ONLN); 
  // Setup the parser
  graphlab::command_line_options
    clopts("Run the PageRank algorithm.", numCPU);
 
  std::string input_filename; 
  clopts.attach_option("infile", &input_filename,
		  "PageRank input file. In src, dest, weight format.");

  // Parse the command line input
  if(!clopts.parse(argc, argv)) {
	  std::cout << "Error in parsing input." << std::endl;
	  return EXIT_FAILURE;
  }
  if(!clopts.is_set("infile")) {
	  std::cout << "Input file no provided!" << std::endl;
	  clopts.print_description();
	  return EXIT_FAILURE;
  }

  // Create a graphlab core
  gl_types::core core;

  // Set the engine options
  core.set_engine_options(clopts);
  
  // Create a synthetic graph
  create_graph(core.graph(), input_filename);

  // Schedule all vertices to run pagerank update on the
  // first round.
  core.add_task_to_all(pagerank_update, 100.0);
  
  // Run the engine
  double runtime = core.start();

  // We are done, now output results.
  std::cout << "Graphlab finished, runtime: " << runtime << " seconds." << std::endl;
  
  std::cout << "wasted ratio: " << wasted*(1.0)/total << std::endl; 

  fp = fopen("result.txt", "a");
  fprintf(fp, "%s: wasted ratio: %.2f%%\n", input_filename.c_str(), wasted*(100.0)/total);
  fclose(fp);
  // First we need to compute a normalizer. This could be
  // done with the sync facility, but for simplicity, we do
  // it by hand.
  double norm = 0.0;
  for(graphlab::vertex_id_t vid=0; vid<core.graph().num_vertices(); vid++) {
  	 norm += core.graph().vertex_data(vid).value;
  }
  
  // And output 5 first vertices pagerank after dividing their value
  // with the norm.
  for(graphlab::vertex_id_t vid=0; vid<5 && vid<core.graph().num_vertices(); vid++) {
  	 std::cout << "Page " << vid << " pagerank = " <<
  	 	core.graph().vertex_data(vid).value / norm << std::endl;
  }
  
	  
  return EXIT_SUCCESS;
} // End of main


  // Configuration information
// std::string filename;

//   clopts.attach_option("infile", &filename,
//                        "PageRank input file. In src, dest, weight format.");
// 
//   // Parse the command line input
//   if(!clopts.parse(argc, argv)) {
//     std::cout << "Error in parsing input." << std::endl;
//     return EXIT_FAILURE;
//   }
//   if(!clopts.is_set("infile")) {
//     std::cout << "Input file no provided!" << std::endl;
//     clopts.print_description();
//     return EXIT_FAILURE;
//   }

// Load the graph
//   if(!load_graph(filename, core.graph())) {
//     std::cout << "Error in parsing graph!" << std::endl;
//     return EXIT_FAILURE;
//   }



/**
 * Load a graph file specified in the format:
 *
 *   source_id, target_id, weight
 *   source_id, target_id, weight
 *   source_id, target_id, weight
 *               ....
 *
 * The file should not contain repeated edges.
 */
bool load_graph(const std::string& filename,
                pagerank_graph& graph) {
  std::ifstream fin(filename.c_str());
  if(!fin.good()) return false;
  // Loop through file reading each line
  while(fin.good()) {
    size_t source = 0;
    size_t target = 0;
    float weight = -1;
    fin >> source;
    fin.ignore(1); // skip comma
    fin >> target;
    fin.ignore(1); // skip comma
    fin >> weight;
    //    fin.ignore(1); // skip comma
    // Ensure that the number of vertices is correct
    if(source >= graph.num_vertices() ||
       target >= graph.num_vertices())
      graph.resize(std::max(source, target) + 1);
    if(source != target) {
      // Add the edge
      edge_data edata(weight);
      graph.add_edge(source, target, weight);
    } else {
      // add the self edge by updating the vertex weight
      graph.vertex_data(source).self_weight = weight;
    }       
  }
  graph.finalize();
  return true;
} // end of load graph
