#include "kdtreeKISS.hpp"

const double TWO = 2.0;
const int SAMPLESIZE = 100;

bool nodeComparison::operator()(const node& a, const node& b) {
    return (a.coordinates[myAxis] < b.coordinates[myAxis]);
}

bool nodeGIDComparison::operator()(const node& a, const node& b) {
    return (a.GID < b.GID);
}

bool nodeRangeComparison::operator()(const node& a, const node& b){
		return (abs(a.coordinates[myAxis] - b.coordinates[myAxis]) < myHorizon);
}

// Adapted from the partition method associated with the quicksort algorithm article from wikipedia
// After calling this method, the container nodeer to by the iterators becomes: [less than median ... median ... greater than median]

std::vector<node >::iterator Kdtree::partitionNodeList( std::vector<node >::iterator begining,  std::vector<node >::iterator ending,  std::vector<node >::iterator pivotIndex, nodeComparison& comparator) {
    // The estimated median node is placed at the end of the vector
    std::swap(*(pivotIndex), *(ending-1));
     std::vector<node >::iterator storeIndex = begining;

    // Compare remaining nodes against pivot
    for( std::vector<node >::iterator it = begining; it != (ending-1); ++it) {

        if(comparator(*it, *(ending-1))) {
            std::swap(*it, *storeIndex);
            storeIndex++;
        }
    }

    //move pivot to final place roughly in the middle of the vector.
    std::swap(*storeIndex, *(ending-1));
    return storeIndex;
}

// This method samples some values along the selected axis from the node vector and finds the location of the median of those samples
// estimating that it is similar to the median of the whole vector.

std::vector<node >::iterator Kdtree::estimateMedian( std::vector<node >::iterator begining,  std::vector<node >::iterator ending, nodeComparison& comparisonFunctionObject) {
    int size = std::distance(begining, ending);
    int tempIndex = 0;
     std::vector<node > temp(SAMPLESIZE);
     std::vector<int> indices(SAMPLESIZE);

    for( std::vector<node >::iterator it = temp.begin(); it != temp.end(); ++it, ++tempIndex) {
        indices[tempIndex] = rand()%size;
        temp[tempIndex] = *(begining + indices[tempIndex]);
    }
    tempIndex = 0;

    std::sort (temp.begin(), temp.end(), comparisonFunctionObject);

    for( std::vector<node >::iterator it = temp.begin(); it != temp.end(); ++it, ++tempIndex) {
        // I dont wanna write another comparator, so use the one from earlier creatively.
        if(comparisonFunctionObject(temp[SAMPLESIZE/2], *(begining + indices[tempIndex]))) {
            asm("nop");
        }
        else if(comparisonFunctionObject(*(begining + indices[tempIndex]), temp[SAMPLESIZE/2])) {
            asm("nop");
        }
        else
            return( begining + indices[tempIndex] );
    }
    return begining; //Error
}

std::vector<node >::iterator Kdtree::computeMedian( std::vector<node >::iterator begining,  std::vector<node >::iterator ending, nodeComparison& comparisonFunctionObject) {
    std::sort (begining, ending, comparisonFunctionObject);
		std::vector<node>::iterator median = begining;
		std::advance(median, std::distance(begining, ending)/2);
		return median;
}

int Kdtree::buildTree( std::vector<node >::iterator begining,  std::vector<node >::iterator ending, unsigned int depth) {

    // split the node set we've been given along axis
    int axis = depth % 3; //This is class is meant for 3 spatial dimensions (but many non spatial dimensions are ok) only partially to make development easier

    // because we have AoS, and we're only comparing one coordinate, we need a fancy comparison function
    // since we're limited to binary operators
    nodeComparison comparisonFunctionObject(axis);

    //estimate median by sampling
    //std::vector<node >::iterator median = Kdtree::estimateMedian(begining, ending, comparisonFunctionObject);

    // Heavyweight method
    //median = Kdtree::partitionNodeList(begining, ending, median, comparisonFunctionObject);
		//
		std::vector<node>::iterator median = Kdtree::computeMedian(begining, ending, comparisonFunctionObject);

    // ID the children if any of the node
		node ORIGIN;
		ORIGIN.coordinates[0] = 0.0;
		ORIGIN.coordinates[1] = 0.0;
		ORIGIN.coordinates[2] = 0.0;

    int mylNeighb, myrNeighb;
    median->mySplitAxis = axis;
		median->distanceFromOrigin = Kdtree::euclidianDistance(*median, ORIGIN);

		// I give nodes the GID for their neighbors because a test for NULL does not work
    if(std::distance(begining, median) > 0) {
        median->lnode = Kdtree::buildTree( begining, median, depth+1);
    }
    else
        median->lnode = -99;

    if(std::distance(median+1, ending) > 0) {
        median->rnode = Kdtree::buildTree( median+1, ending, depth+1);
    }
    else
        median->rnode = -99;

    return median->GID;
}

void Kdtree::finalizeTree( std::vector<node >::iterator begin,  std::vector<node >::iterator end) {
    for( std::vector<node >::iterator it = begin; it != end; ++it) {
        if(it->lnode != -99) {
            for( std::vector<node >::iterator lnode = begin; lnode != end; ++lnode) {
                if(it->lnode == lnode->GID) {
                    it->myLesserNeighbor = &*(lnode);
                }
            }
        }

        if(it->rnode != -99) {
            for( std::vector<node >::iterator rnode = begin; rnode != end; ++rnode) {
                if(it->rnode == rnode->GID) {
                    it->myGreaterNeighbor = &*(rnode);
                }
            }
        }
    }
}

void Kdtree::printTreeVasyaNovikov(std::string prefix, const node& current, bool isTail) {
    std::cout << prefix + (isTail ? "└── " : "├── ") + "" << current.GID << " " << current.mySplitAxis << " (" << current.coordinates[0] << ", " << current.coordinates[1] << ", " << current.coordinates[2] << ")" << std::endl;

    if(current.lnode != -99)
        if(current.rnode != -99) {
            Kdtree::printTreeVasyaNovikov(prefix + (isTail ? "    " : "│   "), *(current.myLesserNeighbor), false);
            Kdtree::printTreeVasyaNovikov(prefix + (isTail ? "    " : "│   "), *(current.myGreaterNeighbor), true);
        }
        else
            Kdtree::printTreeVasyaNovikov(prefix + (isTail ? "    " : "│   "), *(current.myLesserNeighbor), true);
    else if(current.rnode != -99) {
        Kdtree::printTreeVasyaNovikov(prefix + (isTail ? "    " : "│   "), *(current.myGreaterNeighbor), true);
    }
    //don't do anything if node has no children
} 


double Kdtree::euclidianDistance(const node& nodeA, const node& nodeB) {
// Accumulates the squared differences of elements from the locations of the nodes and takes the sqaure root.
    return pow(nodeA.coordinates[0] - nodeB.coordinates[0], TWO) + 
					pow(nodeA.coordinates[1] - nodeB.coordinates[1], TWO) + 
					pow(nodeA.coordinates[2] - nodeB.coordinates[2], TWO);
}

void Kdtree::searchTree(const node& focus, const node& current, std::vector<int>& neighborhood, const double horizon) {
    double splitAxisVector = focus.coordinates[current.mySplitAxis] - current.coordinates[current.mySplitAxis];
    // If the lhs of the splitting axis isn't empty and could possibly intersect the neighborhood sphere search the lhs
    if(current.lnode != -99){
			if(splitAxisVector < horizon){
    		Kdtree::searchTree(focus, *current.myLesserNeighbor, neighborhood, horizon);
			}
			else
				asm("nop");
		}
		else
			asm("nop");

    // If the rhs of the splitting axis isn't empty and could possibly intersect the neighborhood sphere search the rhs
    if(current.rnode != -99){
			if(splitAxisVector > -horizon){
       Kdtree::searchTree(focus, *current.myGreaterNeighbor, neighborhood, horizon);
			}
			else
				asm("nop");
		}
		else 
			asm("nop");

		if(Kdtree::euclidianDistance(focus, current) < horizon*horizon){
   		neighborhood.push_back(current.GID);
		}
				
}

void Kdtree::buildNeighborhoods(const node& root, std::vector<node>& nodes, std::vector<std::vector<int> >& neighborhoods, const double horizon) {
    for(std::vector<node>::iterator focus = nodes.begin(); focus != nodes.end(); ++focus) {
        Kdtree::searchTree(*focus, root, neighborhoods[std::distance(nodes.begin(), focus)], horizon);
    }
}

