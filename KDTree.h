#ifndef __kdtree__
#define __kdtree__
#include <stack>
#include <string>
#include <vector>

class KDNode {
    private:
        double lat;         // latitude
        double lon;         // longitude
        std::string name;   // description
        KDNode* left;       // pointer to left child
        KDNode* right;      // pointer to right child
        
        /* provided by instructor */
        // calculate the distance in miles between this object and the location
        // given by _la and _lo
        double distance(double _la, double _lo);
    
    public:
        // constructor / destructor
        KDNode(double la, double lo, const std::string &desc);
        ~KDNode();
    
    friend class KDTree;
};

class KDTree {
    private:
        unsigned int size;  // number of nodes in the tree
        KDNode* root;       // pointer to the root of the tree
        
        // method called by the destructor (deletes all nodes)
        // postorder traversal
        void destroy(KDNode* p);
        
        // recursive helper method k-dimensionally inserts a node into the tree
        KDNode* insert(KDNode* p, int depth, double la, double lo, const std::string &desc);
        
        // recursive helper method prunes the partitions of the tree for
        // all nodes under a distance 'rad' from 'la,lo'
        void rangeQuery(KDNode* p, std::stack<KDNode*>* query, unsigned int* count, unsigned int* comparisons, int depth, double la, double lo, double rad, const std::string &filter);
    	
    	// graham scan algorithm
    	void grahamScan(std::vector<KDNode*> &A, std::stack<KDNode*> &convex_points);
    	// find hull
    	void findHull(std::vector<KDNode*> &A);
    	// determine ccw turn
    	int ccw(KDNode* a, KDNode* b, KDNode* c);
    	    	
    	// merge sort sub-methods
    	void r_mergesort(std::vector<KDNode*> &A, std::vector<KDNode*> &aux, unsigned int lo, unsigned int hi);
    	void merge(std::vector<KDNode*> &A, std::vector<KDNode*> &aux, unsigned int lo, unsigned int mid, unsigned int hi);
    	// merge sort algorithm
        void mergesort(std::vector<KDNode*> &A, unsigned int n);
    
    public:
        // constructor / destructor
        KDTree();
        ~KDTree();
        
        // insert a new node at the beginning of the list
        void insert(double la, double lo, const std::string &desc);
        // print all the nodes under a distance 'rad' from 'la,lo' and where
        // filter is a non-empty substring of their description
        std::vector<KDNode*> printNeighbors(double la, double lo, double rad, const std::string &filter);
        // print all the nodes on the convex hull
        // filter is a non-empty substring of their description
        unsigned int printConvexHull(double la, double lo, double rad, const std::string &filter);
        // return the number of nodes in the list
        unsigned int getSize();
};

#endif
