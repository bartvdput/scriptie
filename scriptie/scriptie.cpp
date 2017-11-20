// CheckCriteria.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include "stdafx.h"

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/GraphCopyAttributes.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/layered/SugiyamaLayout.h>
#include <ogdf/layered/OptimalRanking.h>
#include <ogdf/layered/MedianHeuristic.h>
#include <ogdf/layered/OptimalHierarchyLayout.h>
#include <ogdf/orthogonal/OrthoLayout.h>
#include <ogdf/planarity/EmbedderMinDepthMaxFaceLayers.h>
#include <ogdf/planarity/PlanarSubgraphFast.h>
#include <ogdf/planarity/PlanarizationLayout.h>
#include <ogdf/planarity/SubgraphPlanarizer.h>
#include <ogdf/planarity/VariableEmbeddingInserter.h>
#include <ogdf/planarity/PlanarizationGridLayout.h>
#include <ogdf/planarity/PlanRep.h>
#include <ogdf/planarity/PlanarSubgraphBoyerMyrvold.h>
#include <ogdf/module/CrossingMinimizationModule.h>
#include <ogdf/planarity/MaximumPlanarSubgraph.h>

using namespace ogdf;
using namespace std;

void CreateGraph(Graph&, GraphAttributes&);
void CreateGraphTwo(Graph&, GraphAttributes&);

void BendPromotion(Graph&, GraphAttributes&);
void SetGraphLayout(Graph&, GraphAttributes&);
double BendCriterium(Graph&, GraphAttributes&);
double CrossingCriterium(Graph&, GraphAttributes&, double);
double NodeOrthogonalityCriterium(Graph&, GraphAttributes&);
double EdgeOrthogonalityCriterium(Graph&, GraphAttributes&);
void OrthogonalLayout(Graph&, GraphAttributes&);
void PlanarRepresentation(Graph&, GraphAttributes&);
void test();

const double NODE_SIZE = 20.0;
const double PI = 3.141592653589793238463;

// now hardcoded for test graph, need to be able to get this value from my own layout-algorithm
const double CROSSINGS = 2;

int main() {
	//test();

	Graph test;
	GraphAttributes testAttributes(test, GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::nodeLabel | GraphAttributes::nodeStyle | GraphAttributes::edgeType | GraphAttributes::edgeArrow | GraphAttributes::edgeStyle);

	//CreateGraphTwo(test, testAttributes);
	randomSimpleGraph(test, 20, 40);
	SetGraphLayout(test, testAttributes);
	PlanarRepresentation(test, testAttributes);

	// construct graph
	Graph graph;
	GraphAttributes graphAttributes(graph, GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics |	GraphAttributes::nodeLabel | GraphAttributes::nodeStyle |	GraphAttributes::edgeType |	GraphAttributes::edgeArrow | GraphAttributes::edgeStyle	);

	Graph graphCopy;
	GraphAttributes graphAttributesCopy(graphCopy, GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::nodeLabel | GraphAttributes::nodeStyle | GraphAttributes::edgeType | GraphAttributes::edgeArrow | GraphAttributes::edgeStyle);

	// set layouts for graph 
	CreateGraph(graph, graphAttributes);
	SetGraphLayout(graph, graphAttributes);

	CreateGraph(graphCopy, graphAttributesCopy);
	SetGraphLayout(graphCopy, graphAttributesCopy);

	// test criteria
	double Nb = BendCriterium(graphCopy, graphAttributesCopy);
	double Nc = CrossingCriterium(graphCopy, graphAttributesCopy, CROSSINGS);
	double Nno = NodeOrthogonalityCriterium(graphCopy, graphAttributesCopy);
	double Neo = EdgeOrthogonalityCriterium(graphCopy, graphAttributesCopy);

	cout << "Criteria" << endl << endl;
	cout << "Crossings (N_c): " << Nc << endl;
	cout << "Bends (N_b): " << Nb << endl;
	cout << "Node Ortho (N_no): " << Nno << endl;
	cout << "Edge Ortho (N_eo): " << Neo << endl;
	cout << endl << endl;

	//draw graphs
	GraphIO::write(graphAttributes, "C:\\Users\\Bart\\Desktop\\Output.svg", GraphIO::drawSVG);
	GraphIO::write(graphAttributesCopy, "C:\\Users\\Bart\\Desktop\\Output2.svg", GraphIO::drawSVG);

	return 0;
}

void test() {
	Graph G;
	GraphAttributes GA(G, GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::nodeLabel | GraphAttributes::nodeStyle | GraphAttributes::edgeType | GraphAttributes::edgeArrow | GraphAttributes::edgeStyle);

	randomSimpleGraph(G, 100, 150);
	SubgraphPlanarizer SP;
	SP.setSubgraph(new PlanarSubgraphFast<int>);
	SP.setInserter(new VariableEmbeddingInserter);
	int crossNum;
	PlanRep PR(G);
	SP.call(PR, 0, crossNum);
	cout << crossNum << " crossings" << endl;
	cout << isPlanar(G) << " planar" << endl;
	GraphIO::write(PR, "C:\\Users\\Bart\\Desktop\\test.gml", GraphIO::writeGML);
}

void PlanarRepresentation(Graph& G, GraphAttributes& GA) {
	cout << "Edges before: " << G.edges.size() << endl;

	//PlanarSubgraphBoyerMyrvold ps;
	List<edge> deletedEdges;
	MaximumPlanarSubgraph ps;
	ps.call(G, deletedEdges);
	
	cout << "Deleted edges: " << endl;
	for (edge e: deletedEdges) {
		cout << "edge: " << e->source() << " --- " << e->target() << endl;
		G.delEdge(e);
	}
	cout << endl << endl;
	cout << "Edges after: " << G.edges.size() << endl;

	cout << endl << endl;
	cout << endl << endl;
	cout << "IS PLANAR? " << isPlanar(G) << endl;

	PlanarizationGridLayout pgl;

	pgl.call(GA);
	GraphIO::write(GA, "C:\\Users\\Bart\\Desktop\\Output4.svg", GraphIO::drawSVG);

	cout << "Number of crossings: " << pgl.numberOfCrossings() << endl;
}

void OrthogonalLayout(Graph& G, GraphAttributes& GA) {
	PlanarizationGridLayout pgl;
	pgl.call(GA);
	GraphIO::write(GA, "C:\\Users\\Bart\\Desktop\\Output3.svg", GraphIO::drawSVG);
}

void BendPromotion(Graph& G, GraphAttributes& GA) {
	List<edge> edges;
	G.allEdges(edges);

	while (!edges.empty()) {
		edge e = edges.popFrontRet();
		DPolyline bends_e = GA.bends(e);
		node s = e->source();
		node t = e->target();

		//check if an edge has bendpoints
		if (!bends_e.empty()) {
			while (!bends_e.empty()) {
				DPoint p = bends_e.front();

				//insert new node
				node n = G.newNode();
				GA.x(n) = p.m_x;
				GA.y(n) = p.m_y;

				edge e_ = G.newEdge(s, n);
				GA.arrowType(e_) = ogdf::EdgeArrow::None;
				GA.strokeColor(e_) = Color("#0000FF");
				s = n;

				bends_e.popFront();
			}
			edge e_ = G.newEdge(s, t);
			GA.arrowType(e_) = ogdf::EdgeArrow::None;
			GA.strokeColor(e_) = Color("#0000FF");

			G.delEdge(e);
		}
	}
}

/*
 * calculate bend criterium
 */
double BendCriterium(Graph& G, GraphAttributes& GA) {
	// original amount of edges
	double m = G.edges.size();

	/*
	 * Maybe get this out of this function, because this is not the only criterium that uses bendpromotion
	 */
	BendPromotion(G, GA);

	// amount of edges after bend promotion
	double m2 = G.edges.size();

	double b_avg = (m2 - m) / m2;
	double Nb = 1 - b_avg;

	return Nb;
}

/*
 * calculate crossing criterium
 */
double CrossingCriterium(Graph& G, GraphAttributes& GA, double c) {
	// calculate all (theoretically) possible edge crossings
	double m2 = G.edges.size();
	double c_all = (m2 * (m2 - 1)) / 2;

	// calculate imopssible edge crossings
	double c_impossible = 0;

	for (node n : G.nodes) {
		double sum = n->degree() * (n->degree() - 1);
		c_impossible += sum;
	}

	c_impossible = c_impossible / 2;

	// calculate real max possible crossings
	double c_max = c_all - c_impossible;

	// calculate crossings criterium
	double Nc = 1;
	if (c_max > 0) {
		Nc -= (c / c_max);
	}

	return Nc;
}

/*
 * calculate edge orthogonality criterium
 */
double EdgeOrthogonalityCriterium(Graph& G, GraphAttributes& GA) {
	double sum = 0;

	for (edge e : G.edges) {
		node s = e->source();
		node t = e->target();

		double s_x = GA.x(s);
		double s_y = GA.y(s);

		double t_x = GA.x(t);
		double t_y = GA.y(t);

		double angle = abs(atan2(s_y - t_y, t_x - s_x) * 180 / PI);
		Array<double> values = { angle, abs(90.0 - angle), abs(180.0 - angle) };
		double *deviation = min_element(begin(values), end(values));
		*deviation = *deviation / 45;

		sum += *deviation;
	}
	
	double Neo = 1.0 - ((1.0 / G.edges.size()) * sum);

	return Neo;
}

/*
 * calculate node orthogonality criterium
 *
 * TODO: calculate value '2 * NODE_SIZE' from: 
 * GCD of the set of vertical and horizontal (pixel) differences between all geometrically adjacent nodes.
 */
double NodeOrthogonalityCriterium(Graph& G, GraphAttributes& GA) {
	double width = 0, height = 0;
	for (node n : G.nodes) {
		double nodeWidth = GA.x(n) / (2 * NODE_SIZE);
		double nodeHeight = GA.y(n) / (2 * NODE_SIZE);

		if (nodeWidth > width)
			width = nodeWidth;

		if (nodeHeight > height)
			height = nodeHeight;
	}

	double A = (1 + width)*(1 + height);
	double Nno = G.nodes.size() / A;

	return Nno;
}

// create testGraph to test criteria imlementations
void CreateGraph(Graph& G, GraphAttributes& GA) {
	// add nodes
	node zero = G.newNode();
	node one = G.newNode();
	node two = G.newNode();
	node three = G.newNode();
	node four = G.newNode();

	// set node positions
	GA.x(zero) = 4 * NODE_SIZE;
	GA.y(zero) = 0;

	GA.x(one) = 4 * NODE_SIZE;
	GA.y(one) = 4 * NODE_SIZE;

	GA.x(two) = 0;
	GA.y(two) = 2 * NODE_SIZE;

	GA.x(three) = 4 * NODE_SIZE;
	GA.y(three) = 8 * NODE_SIZE;

	GA.x(four) = 0;
	GA.y(four) = 8 * NODE_SIZE;

	// add edges
	edge zero_one = G.newEdge(zero, one);
	edge zero_three = G.newEdge(zero, three);
	edge zero_four = G.newEdge(zero, four);
	edge one_two = G.newEdge(one, two);
	edge one_three = G.newEdge(one, three);
	edge two_three = G.newEdge(two, three);

	DPolyline &p = GA.bends(zero_three);
	p.pushBack(DPoint(6 * NODE_SIZE, 2 * NODE_SIZE));
	p.pushBack(DPoint(6 * NODE_SIZE, 6 * NODE_SIZE));
}

void SetGraphLayout(Graph& G, GraphAttributes& GA) {
	for (node v : G.nodes) {
		GA.height(v) = NODE_SIZE; // set the height to 20.0
		GA.width(v) = NODE_SIZE; // set the width to 40.0
		GA.shape(v) = ogdf::Shape::Rect;
		//GA.label(v) = v->index();
	}

	for (edge e : G.edges) {// set default edge color and type
		GA.arrowType(e) = ogdf::EdgeArrow::None;
		GA.strokeColor(e) = Color("#bababa");
	}
}

// create testGraph to test criteria imlementations
void CreateGraphTwo(Graph& graph, GraphAttributes& GA) {
	// add nodes
	node A = graph.newNode();	
	node B = graph.newNode();	
	node C = graph.newNode();	
	node D = graph.newNode();	
	node E = graph.newNode();	
	node F = graph.newNode();	
	node G = graph.newNode();	
	node S = graph.newNode();
	node T = graph.newNode();

	GA.label(A) = "A";
	GA.label(B) = "B";
	GA.label(C) = "C";
	GA.label(D) = "D";
	GA.label(E) = "E";
	GA.label(F) = "F";
	GA.label(G) = "G";
	GA.label(S) = "S";
	GA.label(T) = "T";

	// add edgraphes
	edge AB = graph.newEdge(A, B);
	edge AC = graph.newEdge(A, C);
	edge AG = graph.newEdge(A, G);
	edge BC = graph.newEdge(B, C);
	edge BG = graph.newEdge(B, G);
	edge CS = graph.newEdge(C, S);
	edge CT = graph.newEdge(C, T);
	edge DE = graph.newEdge(D, E);
	edge DT = graph.newEdge(D, T);
	edge ET = graph.newEdge(E, T);
	edge FS = graph.newEdge(F, S);
	edge GS = graph.newEdge(G, S);

	for (edge e : graph.edges) {// set default edge color and type
		GA.arrowType(e) = ogdf::EdgeArrow::None;
		GA.strokeColor(e) = Color("#0000FF");
	}
}