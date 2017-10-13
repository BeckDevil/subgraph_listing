#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "graph.h"
#include "wtime.h"

#define CLR_BLACK 0
#define CLR_GRAY 1
#define CLR_WHITE 2

#define UNVISITED 0
#define VISITED 1

struct TCB{
	unsigned long long int embeddings;
	unsigned long long int recursive_calls;
	uint8_t *visited;
};

typedef std::vector<int> Gpsi;
typedef std::vector<std::vector<int> > Intermediate;

class Psgl{
	public:
		int num_thrds;
		graph *data_graph, *query_graph;
		const char* data_graph_filename, *query_graph_filename, *preset_filename;
		bool break_automorph;
	
		unsigned long long int total_recursive_calls;
		unsigned long long int total_embedding_found;

		uint8_t *matching_order;
		uint8_t *matching_order_map;
		uint8_t *automorph_group_id;
		uint8_t *color;
		uint8_t *visited;

		TCB *myTCB;

		std::vector<std::vector<int> > all_gpsi;
	public:
		Psgl(const char *data_graph_file, const char *query_graph, bool break_auto, const char* preset_file);
		void generic_query_proc();
		void expand_instance(Gpsi &my_gpsi);
		bool process_neighbor(std::vector<int> my_gpsi, int vp_prime, int current_dnode, std::vector<int>& cand_list);
		bool get_candidate_set(std::vector<int> my_gpsi, int vp_prime, int current_dnode, std::vector<int>& candidates);
		std::vector<int> intersection(std::vector<int> &v1, std::vector<int> &v2);
		void read_presets();
};
