#ifndef __GRAPH_H__
#define __GRAPH_H__
#include <iostream>
#include <vector>
#include <map>

#include "util.h"
		
struct adjLabelFrequency{
	std::map<int, std::vector<int> > labelCount;
};


class graph
{
    public:
        index_t *beg_pos;
        vertex_t *csr;
        index_t *beg_pos_label;
        vertex_t *csr_label;
        index_t vert_count;
        index_t vert_max;
        vertex_t edge_count;
        index_t label_count;
        index_t *degree;
		//The data structures for storing the inverse label index
		std::map<int, std::vector<int> > labelVertexList;
		std::vector<std::vector<int> > sorted_csr;

		adjLabelFrequency* adjVertices;
    public:
        graph(const char *graphFile);
        ~graph();
        void test();

		//function for generating inverse label indexes
		void buildLabelVertexList();
		void buildVertexLabelVertexList();
		void build_sorted_csr();
		bool isEdge(int x, int y);

		std::map<int, std::vector<int> >* getLabelVertexList(){
			return &labelVertexList;
		}

		void testIndex(int label){
			std::map<int, std::vector<int> >::iterator itr = adjVertices[4].labelCount.begin();
			for(; itr != adjVertices[4].labelCount.end(); itr++){
				std::vector<int>::iterator j = itr->second.begin();
				std::cout << itr->first <<" " << itr->second.size() << "\n";
				for(; j != itr->second.end(); ++j){
					std::cout << itr->second[*j] << "\t";
				}
				std::cout << "\n";
			}
		}
};

#endif

