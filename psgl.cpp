#include<stdio.h>
#include<vector>
#include<map>
#include<limits.h>
#include<stdlib.h>
#include<float.h>
#include<math.h>
#include<mpi.h>
#include<omp.h>
#include<assert.h>
#include<string.h>
#include<fstream>

#include"graph.h"
#include"psgl.h"
#include"wtime.h"

Psgl::Psgl(const char* data_graph_file, const char* query_graph_file, bool break_auto, const char* preset_file)
{
	//read parameters into the file
	data_graph_filename 	= data_graph_file;
	query_graph_filename 	= query_graph_file;
	preset_filename 		= preset_file;
	break_automorph 		= break_auto;

	//create dat and query graphs
	data_graph = new graph(data_graph_filename);
	query_graph = new graph(query_graph_filename);

	// create adjacency lists for 
	data_graph->buildLabelVertexList();
	data_graph->build_sorted_csr();

	query_graph->buildLabelVertexList();
	query_graph->build_sorted_csr();

	total_recursive_calls = 0;
	total_embedding_found = 0;

	matching_order 		= new uint8_t[query_graph->vert_count];
	matching_order_map 	= new uint8_t[query_graph->vert_count];
	automorph_group_id 	= new uint8_t[query_graph->vert_count];
	read_presets();
//	for(int i = 0; i < query_graph->vert_count; ++i){
//		printf("%d %d %d\n", matching_order[i], matching_order_map[i], automorph_group_id[i]);
//	}

	num_thrds = 32;
	myTCB = new TCB[num_thrds];
	for(int thid = 0; thid < num_thrds; ++thid){
		myTCB[thid].embeddings = 0;
		myTCB[thid].recursive_calls = 0;
		myTCB[thid].visited = new uint8_t[query_graph->vert_count];
		memset(myTCB[thid].visited, UNVISITED, query_graph->vert_count);
	}

	color			= new uint8_t[query_graph->vert_count];
	memset(color, CLR_WHITE , query_graph->vert_count);
	
	visited 		= new uint8_t[query_graph->vert_count];
	memset(visited, UNVISITED, query_graph->vert_count);
	// Vector of GPSIs
	MPI_Barrier(MPI_COMM_WORLD);
}

void Psgl::generic_query_proc()
{
	// TODO: Initialization
	double time;
	int current_qnode = matching_order[0];
	for(int i = 0; i < data_graph->vert_count; ++i){
		// If degree is not satisfied, skip
		if(data_graph->degree[i] < query_graph->degree[current_qnode]){
			continue;
		}
		//create a new Gpsi for satisfying candidates
		Gpsi my_gpsi;
		my_gpsi.push_back(i);
		all_gpsi.push_back(my_gpsi);
	}
	printf("Total initial Gpsis is %d\n", all_gpsi.size());
	// Each data vertex creates a Gpsi, contains 1 match 
	// Initial set of Gpsi, starting nodes
	
	while(!all_gpsi.empty()){
		time = wtime();
		std::vector< std::vector<int> > next_gpsi(all_gpsi.begin(), all_gpsi.end());
		all_gpsi.clear();
#pragma omp parallel num_threads(num_thrds)
		{
			double xtime = omp_get_wtime();
	//		printf("Hello, this is new level expansion\n");
#pragma omp for schedule(dynamic)
			for(int itr = 0; itr < next_gpsi.size(); ++itr){
				int thid = omp_get_thread_num();
//				printf("Thread %d\n", thid);
				std::vector<int> my_gpsi = next_gpsi[itr];
//				for(int i = 0; i < my_gpsi.size(); ++i){
//					printf("%d ", my_gpsi[i]);
//				}
//				printf("\n");
				expand_instance(my_gpsi);
			}
		}
	}

	for(int i = 0; i< num_thrds; ++i){
		total_embedding_found += myTCB[i].embeddings;
	}
	printf("Total %d Embeddings are found in %f seconds\n", total_embedding_found, wtime()-time);
}

// Individual expansion instances 
void Psgl::expand_instance(Gpsi& my_gpsi)
{
	// Get currently expanding v_p and v_d from Gpsii
	Intermediate my_frontier;
	int thid = omp_get_thread_num();
	int current_qnode = matching_order[my_gpsi.size()-1];
	color[current_qnode] = CLR_BLACK;
	myTCB[thid].visited[current_qnode] = VISITED;
	//int current_dnode = my_gpsi[my_gpsi.size()-1];

//	printf("The size of embedding is %d\n", my_gpsi.size());

	// set candlist to NULL, color of current qnode to black
	std::vector<int> cand_list;
	//This is lame
	for(int i = 0; i< data_graph->vert_count; ++i){
		cand_list.push_back(i);
	}
	
//	printf("candidate list size is %d \n", cand_list.size());
	int next_qnode = matching_order[my_gpsi.size()];
//	printf("next node to match is %d\n", next_qnode);
	for(int i = query_graph->beg_pos[next_qnode]; i < query_graph->beg_pos[next_qnode+1]; ++i){
		int vp_prime = query_graph->csr[i];
		if(myTCB[thid].visited[vp_prime] == VISITED){
			int map_vp_prime = my_gpsi[matching_order_map[vp_prime]];
//			printf("parent %d\n", my_gpsi[matching_order_map[vp_prime]]);
			std::vector<int> neighbors = data_graph->sorted_csr[map_vp_prime];
			cand_list = intersection(neighbors, cand_list);
		}
	}

	for(int citr = 0; citr < cand_list.size(); ++citr){
		int match = cand_list[citr];
//		printf("match is %d\n", match);
		if(data_graph->degree[match] < query_graph->degree[next_qnode]){
			cand_list[citr] = -1;
			continue;
		}
		for(std::vector<int>::iterator mitr = my_gpsi.begin(); mitr != my_gpsi.end(); ++mitr){
			if(match == *mitr){
				cand_list[citr] = -1;
				break;
			}
		}
		//apply  automorphism breaking rule here
		for(int itr = 0; itr <  my_gpsi.size(); ++itr){
			int neighbor = matching_order[itr];
			if(automorph_group_id[neighbor] == automorph_group_id[next_qnode]){
				if((neighbor < next_qnode && my_gpsi[matching_order_map[neighbor]] > match) || (neighbor > next_qnode && my_gpsi[matching_order_map[neighbor]] < match)){
					cand_list[citr] = -1;
					break;
				}
			}
		}
	}

//	printf("The size of candidate list is %d", cand_list.size());

	if(my_gpsi.size() == query_graph->vert_count-1){
		// append the new number at the end and print the output
//		printf("new embedding: %d %d\n", my_gpsi[0], my_gpsi[1]);
		for(std::vector<int>::iterator cnode = cand_list.begin(); cnode != cand_list.end(); ++cnode){
			if(*cnode == -1)
				continue;
			my_gpsi[query_graph->vert_count - 1] = *cnode;
		printf("new embedding: %d %d %d\n", my_gpsi[0], my_gpsi[1], my_gpsi[2]);
//			printf("%d\n", my_gpsi[2]);
			myTCB[thid].embeddings++;
		}
	}
	else{
		int x = 0;
		int cnode = 0;
		for(x = 0; x < cand_list.size()-1; ++x){
			cnode = cand_list[x];
			if(cnode == -1)
				continue;
			std::vector<int> new_gpsi(my_gpsi.begin(), my_gpsi.end());
			new_gpsi.push_back(cnode);
			my_frontier.push_back(new_gpsi);
		}
		cnode = cand_list[x];
		if(cnode != -1){
			my_gpsi.push_back(cnode);
			my_frontier.push_back(my_gpsi);
		}
	}
#pragma omp critical
	all_gpsi.insert(all_gpsi.end(), my_frontier.begin(), my_frontier.end());

}

/*
// To make sure that each neighbor have a matching candidate in sight
bool Psgl::process_neighbor(std::vector<int> my_gpsi, int vp_prime, int current_dnode, std::vector<int>& cand_list)
{
	// This function need to return the updated cand_list
	// GRAY means the neighbor is next in frontier
	if(visited[vp_prime] == VISITED){
		// check if map(vp_prime) is in N(current_dnode)
		int map_vp_prime = my_gpsi[matching_order_map[vp_prime]];
		if(!data_graph->isEdge(map_vp_prime, target) || current_dnode == map_vp_prime){
			return false;
		}
		std::vector<int> neighbors = data_graph->sorted_csr[map_vp_prime];
		cand_list = intersection(cand_list, neighbors);
		// cand_list == intersection 
	}
	return true;
}

bool Psgl::get_candidate_set(std::vector<int> my_gpsi, int vp_prime, int current_dnode, std::vector<int>& candidates)
{
	for(int i = data_graph->beg_pos[current_dnode]; i < data_graph->beg_pos[current_dnode+1]; ++i){
		int vd_prime = data_graph->csr[i];
		// Pruning rule 1
		if(data_graph->degree[vd_prime] < query_graph->degree[vp_prime]){// TODO: partial order restriction
			continue;
		}
		// Pruning rule 2
		bool valid = true;
		for(int i = query_graph->beg_pos[vp_prime]; i < query_graph->beg_pos[vp_prime+1]; ++i){
			int vp_pprime = query_graph->csr[i];
			if((color[vp_pprime] == CLR_GRAY) && !data_graph->isEdge(vd_prime, my_gpsi[vp_pprime])){
				valid = false;
				break;
			}
		}
		if(valid == false){
			continue;
		}
		else{
			candidates.push_back(vd_prime);
		}
	}
	if(candidates.empty()){
		return false;
	}
	return true;
}*/

std::vector<int> Psgl::intersection(std::vector<int> &v1, std::vector<int> &v2)
{
	std::vector<int> v3;
    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());
    set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v3));
    return v3;
}

// Read the preset file to hardcode the presets
void Psgl::read_presets()
{
	std::ifstream ifile;
	ifile.open(preset_filename);
	if(!ifile){
		printf("Error in reading file\n");
	}
	int x;

	for(int i = 0; i< query_graph->vert_count; ++i){
		ifile >> x;
		matching_order[i] = x;
	}
	
	for(int i = 0; i < query_graph->vert_count; ++i){
		ifile >> x;
		matching_order_map[i] = x;
	}

	for(int i = 0; i < query_graph->vert_count; ++i){
		ifile >> x;
		automorph_group_id[i] = x;
	}
	ifile.close();
}

