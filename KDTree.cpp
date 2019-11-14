#include "KDTree.h"
#include <math.h>
#include <iostream>
#include <fstream>

KDNode::KDNode(double la, double lo, const std::string &desc) {
    left = NULL;
    right = NULL;
    name = desc;
    lat = la;
    lon = lo;
}

KDNode::~KDNode() {
}

/* provided by instructor */
double KDNode::distance(double _la, double _lo) {
    double param = M_PI / 180.0; 	// required for conversion from degrees to radians
    double rad = 3956.0; 		// radius of earth in miles
    double d_lat = (_la - lat) * param;
    double d_lon = (_lo - lon) * param;
    double dist = sin(d_lat/2) * sin(d_lat/2) + cos(lat*param) * cos(_la*param) * sin(d_lon/2) * sin(d_lon/2);
    dist = 2.0 * atan2(sqrt(dist), sqrt(1-dist));
    return rad * dist;
}

KDTree::KDTree() {
    size = 0;
    root = NULL;
}

KDTree::~KDTree() {
    destroy(root);
}

void KDTree::destroy(KDNode* p) {
    if (p) {    // postorder traversal
        destroy(p->left);
        destroy(p->right);
        delete p;
    }
}

// return the number of nodes in the tree
unsigned int KDTree::getSize() {
    return size;
}

// insert a new node into the tree
void KDTree::insert(double la, double lo, const std::string &desc) {
    // pointer to root of tree with an initial depth of 0
    root = insert(root, 0, la, lo, desc);
}

// recursively insert a new node k-dimensionally (partitioning) into the tree
KDNode* KDTree::insert(KDNode* p, int depth, double la, double lo, const std::string &desc) {
    if (!p) {                               // node is a null
        size++;                             // increment number of nodes in the tree
        return new KDNode(la, lo, desc);    // return a new node
    }
    
    // node has a depth of 0 -- latitude partitioning
    // node has a depth of 1 -- longitude partitioning
    switch (depth) {
        case 0:
            if (la >= p->lat) {     // new latitude is greater than or equal to node->latitude
                p->right = insert(p->right, 1, la, lo, desc);
            } else {                // new latitude is less than node->latitude
                p->left = insert(p->left, 1, la, lo, desc);
            }
            break;
        case 1:
            if (lo >= p->lon) {     // new longitude is greater than or equal to node->longitude
                p->right = insert(p->right, 0, la, lo, desc);
            } else {                // new longitude is less than node->longitude
                p->left = insert(p->left, 0, la, lo, desc);
            }
            break;
    }
    return p;   // return node
}

// print all the nodes under a distance 'rad' from 'la,lo' and where
// filter is a non-empty substring of their description
std::vector<KDNode*> KDTree::printNeighbors(double la, double lo, double rad, const std::string &filter) {
    std::stack<KDNode*> query;
    std::vector<KDNode*> q_arr;
    std::ofstream markers_js;
    unsigned int count = 0, comparisons = 0;
    
    rangeQuery(root, &query, &count, &comparisons, 0, la, lo, rad, filter);
    
    markers_js.open("data/markers.js");
    
    markers_js << "var markers = [\n";
    markers_js << "\t[\"" << "CENTER" << "\", " << la << ", " << lo << "],\n";
    
    q_arr = std::vector<KDNode*>(count);
    for (unsigned int i = 0; i < count; i++) {
    	q_arr[i] = query.top();
    	markers_js << "\t[\"" << q_arr[i]->name << "\", " << q_arr[i]->lat << ", " << q_arr[i]->lon << "],\n";
    	query.pop();
    }

    markers_js << "];\n";
    std::cerr << comparisons << " comparisons were made\n";
    
    return q_arr;
}

// recursively prune and print a node under a distance 'rad' from 'la, lo'
void KDTree::rangeQuery(KDNode* p, std::stack<KDNode*>* query, unsigned int*count, unsigned int*comparisons, int depth, double la, double lo, double rad, const std::string &filter) {
    if (!p) {   // node is null
        return;
    }
    // if the distance between the node and the center is less than the radius, it lies within the radius
    else if (p->distance(la,lo) < rad && p->name.find(filter) != std::string::npos) {
        query->push(p);
        *count += 1;    // increment number of records fetched
    }
    
    *comparisons += 1;  // increment number of comparisions made
    
    // node has a depth of 0 -- latitude partitioning
    // node has a depth of 1 -- longitude partitioning
    switch (depth) {    // ...pruning...
        case 0:
            if (p->distance(la, p->lon) > rad && la < p->lat)           // if farthest radial point east of center is less than node->latitude
                rangeQuery(p->left, query, count, comparisons, 1, la, lo, rad, filter);
            else if (p->distance(la, p->lon) >= rad && la > p->lat)     // if farthest radial point west of center is greater than or equal to node->latitude
                rangeQuery(p->right, query, count, comparisons, 1, la, lo, rad, filter);
            else {                                                      // circular perimeter lies on both sides of partition
                rangeQuery(p->left, query, count, comparisons, 1, la, lo, rad, filter);
                rangeQuery(p->right, query, count, comparisons, 1, la, lo, rad, filter);
            }
            break;
        case 1:
            if (p->distance(p->lat, lo) > rad && lo < p->lon)           // if farthest radial point south of center is less than node->longitude
                rangeQuery(p->left, query, count, comparisons, 0, la, lo, rad, filter);
            else if (p->distance(p->lat, lo) >= rad && lo > p->lon)     // if farthest radial point north of center is greater than or equal to node->longitude
                rangeQuery(p->right, query, count, comparisons, 0, la, lo, rad, filter);
            else {                                                      // circular perimeter lies on both sides of partition
                rangeQuery(p->left, query, count, comparisons, 0, la, lo, rad, filter);
                rangeQuery(p->right, query, count, comparisons, 0, la, lo, rad, filter);
            }
            break;
    }
}

// print all the nodes on the convex hull and where
// filter is a non-empty substring of their description
unsigned int KDTree::printConvexHull(double la, double lo, double rad, const std::string &filter) {
    std::vector<KDNode*> q_arr;					// store queries found in printNeighbors
    std::stack<KDNode*> convex_points;				// store nodes on the convex hull
    std::ofstream convex_js;					// file output
    q_arr = printNeighbors(la, lo, rad, filter);
    
	convex_js.open("data/convex.js");			// create convex.js in ./data/
    convex_js << "var convex = [\n";
    
    // q_arr has enough points to form a polygon/convex hull
    if (q_arr.size() >= 3) {
    	grahamScan(q_arr, convex_points);
    	
    	// print nodes and angles to verify sorting order
    	for (int i = 0; i < q_arr.size(); i++)
    		std::cerr << "\t[\"" << q_arr[i]->name << "\", " << q_arr[i]->lat << ", " << q_arr[i]->lon << "],\t"
    			<< atan2(q_arr[i]->lat - q_arr[0]->lat, q_arr[i]->lon - q_arr[0]->lon) * 180 / M_PI << "\n";
    	
    	// output all convex nodes until stack is empty
    	while (!convex_points.empty()) {
    		KDNode* node = convex_points.top();
    		convex_points.pop();
    	
    		convex_js << "\t{lat: " << node->lat << ", lng: "<< node->lon << "},\n";
    	}
	}
    convex_js << "];\n";
    convex_js.close();
    
    return q_arr.size();
}

// store convex hull nodes from query array into the provided stack
void KDTree::grahamScan(std::vector<KDNode*> &q_arr, std::stack<KDNode*> &convex_points) {
	findHull(q_arr);					// hull is stored at q_arr[0]
	mergesort(q_arr, q_arr.size());		// sorts query by polar angles--O(nlog(n))
	
	convex_points.push(q_arr[0]);		// first two nodes are stored
    convex_points.push(q_arr[1]);
    
    for (unsigned int i = 2; i < q_arr.size(); i++) {
    	KDNode* point_b = convex_points.top();
    	convex_points.pop();
    	while (ccw(convex_points.top(), point_b, q_arr[i]) <= 0) {	// if collinear or
    		point_b = convex_points.top();				// cw, remove point b
    		convex_points.pop();					// ccw(a, b, c)
    	}
    	convex_points.push(point_b);	// restore point b
    	convex_points.push(q_arr[i]);   // insert ccw point
    }
}

// swap location of node with lowest latitude with
// the first node in the query
void KDTree::findHull(std::vector<KDNode*> &q_arr) {
    KDNode* temp = q_arr[0];
    int index_of_anchor = 0;
    // searches for index of node with lowest latitude
    for (int i = 1; i < q_arr.size(); i++) {
    	if (temp->lat > q_arr[i]->lat) {
    		temp = q_arr[i];
			index_of_anchor = i;
		}
		else if (temp->lat == q_arr[i]->lat) {	// break by rightmost lattitude
			if (temp->lon > q_arr[i]->lon) {
				temp = q_arr[i];
				index_of_anchor = i;
			}
		}
	}
	
	temp = q_arr[0];					// temporarily store first result
	q_arr[0] = q_arr[index_of_anchor];	// assign result with lowest latitude to first element
	q_arr[index_of_anchor] = temp;		// return first result to original index of lowest latitude
}

/* Author: Robert Sedgewick & Kevin Wayne
 * Algorithms, 4th Edition by Robert Sedgewick and Kevin Wayne
 * Modified by: Osto Vargas
 */
int KDTree::ccw(KDNode* a, KDNode* b, KDNode* c) {
	// cross product
	double area2 = (b->lon - a->lon)*(c->lat - a->lat) - (b->lat - a->lat)*(c->lon - a->lon);
	if 	(area2 < 0) 	return -1;	// clockwise
	else if (area2 > 0)	return 1;	// counter-clockwise
	else 			return 0;	// collinear
}

/********************************************************************************************************/
// merge sort
void KDTree::merge(std::vector<KDNode*> &q_arr, std::vector<KDNode*> &aux, unsigned int lo, unsigned int mid, unsigned int hi) {
    // copy array
    std::memcpy(&aux[0] + lo, &q_arr[0] + lo, (hi-lo+1)*sizeof(KDNode*));

    // merge
    unsigned int i = lo, j = mid + 1;
    for (int k = lo; k <= hi; k++) {
        if (i > mid) q_arr[k] = aux[j++];
        else if (j > hi) q_arr[k] = aux[i++];
        else if (
        	atan2(aux[j]->lat - q_arr[0]->lat, aux[j]->lon - q_arr[0]->lon) * 180 / M_PI <
        	atan2(aux[i]->lat - q_arr[0]->lat, aux[i]->lon - q_arr[0]->lon) * 180 / M_PI
        	) q_arr[k] = aux[j++];
        else q_arr[k] = aux[i++];
    }
}

void KDTree::r_mergesort(std::vector<KDNode*> &q_arr, std::vector<KDNode*> &aux, unsigned int lo, unsigned int hi) {
    // base case
    if (hi <= lo) return;
    // divide
    unsigned int mid = lo + (hi - lo) / 2;
    // recursively sort halves
    r_mergesort(q_arr, aux, lo, mid);
    r_mergesort(q_arr, aux, mid+1, hi);
    // merge results
    merge(q_arr, aux, lo, mid, hi);
}

void KDTree::mergesort(std::vector<KDNode*> &q_arr, unsigned int n) {
    // allocate space for aux
    std::vector<KDNode*> aux (n);
    // call recursive mergesort
    r_mergesort(q_arr, aux, 1, n - 1);
    // free memory
}
