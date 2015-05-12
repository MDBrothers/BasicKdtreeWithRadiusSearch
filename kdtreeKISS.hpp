#include <algorithm>
#include <iterator>
#include <vector>
#include <vector>
#include <iostream>
#include <string>
#include <math.h>
#include <cstring>

const int CPYSIZE = 3*sizeof(double);

struct node {
    double coordinates[3];
		double distanceFromOrigin;
    int GID, mySplitAxis, rnode, lnode, myBinNum;
    node* myLesserNeighbor, *myGreaterNeighbor;
} temp;

namespace std{
inline	void swap(node& a, node& b){
		temp.GID = b.GID;
		temp.mySplitAxis = b.mySplitAxis;
		temp.rnode = b.rnode;
		temp.lnode = b.lnode;
		temp.distanceFromOrigin = b.distanceFromOrigin;
	//	temp.myLesserNeighbor = b.myLesserNeighbor;
//		temp.myGreaterNeighbor = b.myGreaterNeighbor;
		std::memcpy(temp.coordinates, b.coordinates, CPYSIZE);

		b.GID = a.GID;
		b.mySplitAxis = a.mySplitAxis;
		b.rnode = a.rnode;
		b.lnode = a.lnode;
		b.distanceFromOrigin = a.distanceFromOrigin;
	//	b.myLesserNeighbor = a.myLesserNeighbor;
//		b.myGreaterNeighbor = a.myGreaterNeighbor;
		std::memcpy(b.coordinates, a.coordinates, CPYSIZE);

		a.GID = temp.GID;
		a.mySplitAxis = temp.mySplitAxis;
		a.rnode = temp.rnode;
		a.lnode = temp.lnode;
		a.distanceFromOrigin = temp.distanceFromOrigin;
//		a.myLesserNeighbor = temp.myLesserNeighbor;
//		a.myGreaterNeighbor = temp.myGreaterNeighbor;
		std::memcpy(a.coordinates, temp.coordinates, CPYSIZE);
	}
}

class nodeComparison
{
	private:
		int myAxis;
		double myHorizon;
  public:	
		nodeComparison(int axis){myAxis = axis;};
		bool operator() (const node& a, const node& b);	
};

class nodeGIDComparison{
	public:
		bool operator()(const node& a, const node& b);
};

class nodeRangeComparison{
		private:
			int myAxis;
			double myHorizon;
		public:
			nodeRangeComparison(int axis, double horizon){myAxis = axis; myHorizon=horizon;};
			bool operator()(const node& a, const node& b);
};

namespace Kdtree {
    std::vector<node>::iterator estimateMedian(std::vector<node>::iterator begining,  std::vector<node>::iterator ending, nodeComparison& comparisonFunctionObject);

    std::vector<node>::iterator computeMedian(std::vector<node>::iterator begining,  std::vector<node>::iterator ending, nodeComparison& comparisonFunctionObject);

	  std::vector<node >::iterator partitionNodeList( std::vector<node >::iterator begining,  std::vector<node >::iterator ending,  std::vector<node >::iterator pivotIndex, nodeComparison& comparator);

    int buildTree( std::vector<node>::iterator begin,  std::vector<node>::iterator end, unsigned int depth);

		void finalizeTree(std::vector<node>::iterator begin,  std::vector<node>::iterator end);
    // Based on a solution from stack overflow user named in method
    void printTreeVasyaNovikov(std::string prefix, const node& current, bool isTail);

    double euclidianDistance(const node& nodeA, const node& nodeB);

		void searchTree(const node& focus, const node& current, std::vector<int>& neighborhood, const double horizon); 

    void buildNeighborhoods(const node& root, std::vector<node>& nodes, std::vector<std::vector<int> >& neighborhoods, const double horizon);
}

