// CheckCriteria.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <utility>
#include <iterator>
#include "stdafx.h"
#include "json.hpp"

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
#include <ogdf/graphalg/Clusterer.h>
#include <ogdf/orthogonal/EdgeRouter.h>
#include <ogdf/planarity/VariableEmbeddingInserter.h>
#include <ogdf/planarity/PlanarSubgraphFast.h>
#include <ogdf/planarlayout/MixedModelLayout.h>
#include <ogdf/planarlayout/MMCBDoubleGrid.h>
#include <ogdf/planarlayout/MMCBLocalStretch.h>
#include <ogdf/planarity/FixedEmbeddingInserter.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/energybased/SpringEmbedderGridVariant.h>
#include <ogdf/planarity/PlanarSubgraphCactus.h>
#include <ogdf/planarlayout/SchnyderLayout.h>
#include <ogdf/planarlayout/PlanarDrawLayout.h>
#include <ogdf/planarlayout/PlanarStraightLayout.h>

using namespace ogdf;
using namespace std;
using json = nlohmann::json;

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
void GetDegrees(Graph&, GraphAttributes&);
void CriteriaTesting();
void SubGraphPlan(Graph&, GraphAttributes&);
void addRelations(Graph&, GraphAttributes&);
void ERLayoutAlgorithm(Graph&, GraphAttributes&);
void Planarize(Graph&, GraphAttributes&);
void Embed(Graph&, GraphAttributes&);
void Orthogonalize(Graph&, GraphAttributes&);
int getBiconnectedComponents(Graph&);
void createGraphFromJson(Graph&, GraphAttributes&, string);
void CreateGraphThree(Graph&, GraphAttributes&);

const double NODE_WIDTH = 40.0;
const double NODE_HEIGHT = 10.0;
const float EDGE_STROKEWIDTH = 4;
const double PI = 3.141592653589793238463;

// now hardcoded for test graph, need to be able to get this value from my own layout-algorithm
const double CROSSINGS = 2;

int main() {
	Graph test;
	GraphAttributes testAttributes(test, GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::nodeLabel | GraphAttributes::nodeStyle | GraphAttributes::edgeType | GraphAttributes::edgeArrow | GraphAttributes::edgeStyle | GraphAttributes::edgeLabel);
	
	string file = "entities-big_system.json";	
	createGraphFromJson(test, testAttributes, file);
	SetGraphLayout(test, testAttributes);

	ERLayoutAlgorithm(test, testAttributes);
	//PlanarRepresentation(test, testAttributes);

	//CriteriaTesting();
	return 0;
}

void createGraphFromJson(Graph& G, GraphAttributes& GA, string file) {
	// Read JSON file
	ifstream i(file);
	json js;
	i >> js;

	// map to be able to find nodes with name
	map<string, node> nodes;
	map<string, node>::iterator map_it;

	// create all nodes
	for (size_t i = 0; i < js.size(); i++) {
		string name = js[i]["name"];
		node n = G.newNode();
		GA.label(n) = name;
		nodes.insert(pair<string, node>(name, n));
	}	

	int count = 0;
	int count2 = 0;
	
	// create all edges
	for (size_t i = 0; i < js.size(); i++) {
		// walk through node members
		for (size_t j = 0; j < js[i]["members"].size(); j++) {
			string type = js[i]["members"][j]["relation"];
			// check if edge/relation is found
			if (type != "NONE") {
				// find source node
				string source = js[i]["name"];
				map_it = nodes.find(source);

				// if source node is found, continue
				if (map_it != nodes.end()) {
					// get source node from map
					node s = map_it->second;

					// find target node
					string target = js[i]["members"][j]["type"]["name"];
					map_it = nodes.find(target);

					// if target node is found, continue
					if (map_it != nodes.end()) {
						// get target node from map
						node t = map_it->second;

						if (G.searchEdge(t, s) != 0 && GA.label(s) != GA.label(t)) {
							cout << "Duplicate edge: " << type << endl;
							count2++;
						}

						// check for double edges and self-loops
						if (G.searchEdge(t, s) == 0 && GA.label(s) != GA.label(t)) {
							// make new edge
							edge e = G.newEdge(s, t);
							GA.strokeWidth(e) = 0.5;

							if (type == "UNI_TO_ONE") {
								GA.strokeType(e) = ogdf::StrokeType::Solid;
								GA.arrowType(e) = ogdf::EdgeArrow::Last;
							} else if (type == "BI_MANY_TO_ONE") {
								GA.strokeType(e) = ogdf::StrokeType::Dash;
								GA.arrowType(e) = ogdf::EdgeArrow::First;
							} else if (type == "BI_ONE_TO_MANY") {
								GA.strokeType(e) = ogdf::StrokeType::Dash;
								GA.arrowType(e) = ogdf::EdgeArrow::Last;
							} else if (type == "BI_MANY_TO_MANY") {
								GA.strokeType(e) = ogdf::StrokeType::Dash;
								GA.arrowType(e) = ogdf::EdgeArrow::Both;
							} else if (type == "BI_ONE_TO_ONE") {
								GA.strokeType(e) = ogdf::StrokeType::Dash;
								GA.arrowType(e) = ogdf::EdgeArrow::None;
							}								
						}
					}					
				}
								
			}
		}
	}

	cout << "Count: " << count << endl;
	cout << "Count2: " << count2 << endl;

	// check degree and delete non-connected nodes
	for (map_it = nodes.begin(); map_it != nodes.end(); map_it++) {
		node n = map_it->second;
		if (n->degree() == 0) {
			G.delNode(n);
		}
	}	

	/*
	// List all nodes
	for (node n : G.nodes) {
		cout << "FINAL NODES: " << GA.label(n) << endl;
	}
	
	// List all edges
	cout << endl << endl;
	for (edge e : G.edges) {
		cout << "FINAL EDGES: " << GA.label(e->source()) << " -- " << GA.label(e->target()) << endl;
	}
	*/
}

void CriteriaTesting() {
	// construct graph
	Graph graph;
	GraphAttributes graphAttributes(graph, GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::nodeLabel | GraphAttributes::nodeStyle | GraphAttributes::edgeType | GraphAttributes::edgeArrow | GraphAttributes::edgeStyle);

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
}

int getBiconnectedComponents(Graph& G) {
	EdgeArray<int>& component = EdgeArray<int>(G);
	int count = 0;
	count = biconnectedComponents(G, component);	

	int edgeIndex = 0;
	for (int i = 0; i < count; i++) {
		for (int c : component) {
			if (c == i) {
				cout << "Component: " << i << "  --  " << edgeIndex << endl;
			}
			edgeIndex++;
		}
		edgeIndex = 0;
	}
	
	return count;
}

void ERLayoutAlgorithm(Graph& G, GraphAttributes& GA) {
	PlanarizationGridLayout pl;
	SubgraphPlanarizer *crossMin = new SubgraphPlanarizer;

	// Get a planar subgraph using Boyer Myrvold Algorithm
	PlanarSubgraphBoyerMyrvold *ps = new PlanarSubgraphBoyerMyrvold;
	VariableEmbeddingInserter *ves = new VariableEmbeddingInserter;

	crossMin->setSubgraph(ps);
	crossMin->setInserter(ves);

	MixedModelLayout *ol = new MixedModelLayout;
	MMCBLocalStretch *cb = new MMCBLocalStretch;
	//ol->separation(180.0);

	ol->setCrossingsBeautifier(cb);

	pl.setPlanarLayouter(ol);

	pl.call(GA);

	cout << "number of crossings: " << pl.numberOfCrossings() << endl;

	GraphIO::SVGSettings settings;
	settings.fontSize(2);

	GraphIO::drawSVG(GA, "C:\\Users\\Bart\\Desktop\\Output6.svg", settings);
}

void GetDegrees(Graph&G, GraphAttributes& GA) {
	for (node n : G.nodes) {
		if (GA.shape(n) == Shape::Rect) {
			cout << "Node " << n->index() << " degree: " << n->degree() << endl;
		}
	}
}

void addRelations(Graph& G, GraphAttributes& GA) {
	List<edge> originalEdges;
	G.allEdges(originalEdges);
		
	for (edge& e : originalEdges) {		
		node s = e->source();
		node t = e->target();

		node R = G.newNode();

		G.newEdge(s, R);
		G.newEdge(R, t);

		GA.shape(R) = Shape::Rhomb;

		G.delEdge(e);
	}
}

void PlanarRepresentation(Graph& G, GraphAttributes& GA) {
	PlanarSubgraphBoyerMyrvold ps = PlanarSubgraphBoyerMyrvold();
	List<edge> deletedEdges;

	GraphCopy GC = GraphCopy(G);
	ps.call(G, deletedEdges);

	for (edge e : deletedEdges) {
		G.delEdge(e);
	}

	cout << endl << "deleted: " << deletedEdges.size() << endl;

	//SchnyderLayout *sl = new SchnyderLayout;
	//sl->call(GA);

	//PlanarDrawLayout *pdl = new PlanarDrawLayout;
	//pdl->call(GA);

	PlanarStraightLayout *psl = new PlanarStraightLayout;
	psl->call(GA);

	// set fontsize
	GraphIO::SVGSettings settings;
	settings.fontSize(2);
	// draw graph
	GraphIO::drawSVG(GA, "C:\\Users\\Bart\\Desktop\\Output6.svg", settings);
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
				GA.strokeColor(e_) = Color("#bababa");
				s = n;

				bends_e.popFront();
			}
			edge e_ = G.newEdge(s, t);
			GA.arrowType(e_) = ogdf::EdgeArrow::None;
			GA.strokeColor(e_) = Color("#bababa");

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
		double nodeWidth = GA.x(n) / (2 * NODE_WIDTH);
		double nodeHeight = GA.y(n) / (2 * NODE_HEIGHT);

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
	GA.x(zero) = 4 * NODE_WIDTH;
	GA.y(zero) = 0;

	GA.x(one) = 4 * NODE_WIDTH;
	GA.y(one) = 4 * NODE_HEIGHT;

	GA.x(two) = 0;
	GA.y(two) = 2 * NODE_HEIGHT;

	GA.x(three) = 4 * NODE_WIDTH;
	GA.y(three) = 8 * NODE_HEIGHT;

	GA.x(four) = 0;
	GA.y(four) = 8 * NODE_HEIGHT;

	// add edges
	edge zero_one = G.newEdge(zero, one);
	edge zero_three = G.newEdge(zero, three);
	edge zero_four = G.newEdge(zero, four);
	edge one_two = G.newEdge(one, two);
	edge one_three = G.newEdge(one, three);
	edge two_three = G.newEdge(two, three);

	DPolyline &p = GA.bends(zero_three);
	p.pushBack(DPoint(6 * NODE_WIDTH, 2 * NODE_HEIGHT));
	p.pushBack(DPoint(6 * NODE_WIDTH, 6 * NODE_HEIGHT));
}

void SetGraphLayout(Graph& G, GraphAttributes& GA) {
	GA.setAllHeight(NODE_HEIGHT);
	GA.setAllWidth(NODE_WIDTH);
}

// create testGraph to test criteria imlementations
void CreateGraphTwo(Graph& graph, GraphAttributes& GA) {
	// add nodes
	node Adresses = graph.newNode();	
	node Schools = graph.newNode();	
	node Subjects = graph.newNode();	
	node Parent_Adresses = graph.newNode();	
	node Student_Adresses = graph.newNode();	
	node Parents = graph.newNode();	
	node Student_Parents = graph.newNode();	
	node Teachers = graph.newNode();
	node Classes = graph.newNode();
	node Family_Members = graph.newNode();
	node Students = graph.newNode();
	node Student_Classes = graph.newNode();
	node Families = graph.newNode();
	node Homework = graph.newNode();
	node Reports = graph.newNode();

	GA.label(Adresses) = "Adresses";
	GA.label(Schools) = "Schools";
	GA.label(Subjects) = "Subjects";
	GA.label(Parent_Adresses) = "Parent_Adresses";
	GA.label(Student_Adresses) = "Student_Adresses";
	GA.label(Parents) = "Parents";
	GA.label(Student_Parents) = "Student_Parents";
	GA.label(Teachers) = "Teachers";
	GA.label(Classes) = "Classes";
	GA.label(Family_Members) = "Family_Members";
	GA.label(Students) = "Students";
	GA.label(Student_Classes) = "Student_Classes";
	GA.label(Families) = "Families";
	GA.label(Homework) = "Homework";
	GA.label(Reports) = "Reports";

	// add edgraphes
	edge SchoolsToAdresses = graph.newEdge(Schools, Adresses);
	edge Parent_AdressesToAdresses = graph.newEdge(Parent_Adresses, Adresses);
	edge Parent_AdressesToParents = graph.newEdge(Parent_Adresses, Parents);
	edge Student_AdressesToAdresses = graph.newEdge(Student_Adresses, Adresses);
	edge Student_AdressesToStudents = graph.newEdge(Student_Adresses, Students);
	edge Student_ParentsToParents = graph.newEdge(Student_Parents, Parents);
	edge Student_ParentsToStudents = graph.newEdge(Student_Parents, Students);
	edge TeachersToSchools = graph.newEdge(Teachers, Schools);
	edge ClassesToSubjects = graph.newEdge(Classes, Subjects);
	edge ClassesToTeachers = graph.newEdge(Classes, Teachers);
	edge Family_MembersToParents = graph.newEdge(Family_Members, Parents);
	edge Family_MembersToFamilies = graph.newEdge(Family_Members, Families);
	edge Family_MembersToStudents = graph.newEdge(Family_Members, Students);
	edge Student_ClassesToStudents = graph.newEdge(Student_Classes, Students);
	edge Student_ClassesToClasses = graph.newEdge(Student_Classes, Classes);
	edge FamiliesToParents = graph.newEdge(Families, Parents);
	edge HomeworkToStudents = graph.newEdge(Homework, Students);
	edge ReportsToStudents = graph.newEdge(Reports, Students);	

	for (edge e : graph.edges) {// set default edge color and type
		GA.arrowType(e) = ogdf::EdgeArrow::Last;
		GA.strokeType(e) = ogdf::StrokeType::Solid;
		GA.strokeColor(e) = Color("#bababa");
	}
}

void CreateGraphThree(Graph& G, GraphAttributes& GA) {
	// add nodes
	node ROLES = G.newNode();
	node USERS = G.newNode();
	node MESSAGES = G.newNode();
	node ORDERS = G.newNode();
	node PRODUCTS = G.newNode();
	node ORDER_LINES = G.newNode();

	// add edges
	edge a = G.newEdge(ROLES, USERS);
	edge b = G.newEdge(MESSAGES, USERS);
	edge c = G.newEdge(ORDERS, USERS);
	edge d = G.newEdge(ORDERS, ORDER_LINES);
	edge e = G.newEdge(ORDER_LINES, PRODUCTS);
	//edge f = G.newEdge(ORDER_LINES, ORDERS);
	//edge g = G.newEdge(USERS, ROLES);
}