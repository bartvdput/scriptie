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
double BendCriterion(Graph&, GraphAttributes&);
double CrossingCriterion(Graph&, GraphAttributes&, double);
double AngleCriterion(Graph&, GraphAttributes&);
double NodeOrthogonalityCriterion(Graph&, GraphAttributes&);
double EdgeOrthogonalityCriterion(Graph&, GraphAttributes&);
void CriteriaTesting(Graph&, GraphAttributes&, int);
void addRelations(Graph&, GraphAttributes&);
int ERLayoutAlgorithm(Graph&, GraphAttributes&);
int getBiconnectedComponents(Graph&);
void createGraphFromJson(Graph&, GraphAttributes&, string);

const double NODE_WIDTH = 40.0;
const double NODE_HEIGHT = 20.0;
const float EDGE_STROKEWIDTH = 4;
const double PI = 3.141592653589793238463;

int main() {
	Graph test;
	GraphAttributes testAttributes(test, GraphAttributes::nodeType | GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::nodeLabel | GraphAttributes::nodeStyle | GraphAttributes::edgeType | GraphAttributes::edgeArrow | GraphAttributes::edgeStyle | GraphAttributes::edgeLabel);
	
	// read file
	string file = "entities-big_system.json";
	//cout << "Enter file name to read in: " << endl;	
	//cin >> file;
	createGraphFromJson(test, testAttributes, file);

	// set some layout properties
	SetGraphLayout(test, testAttributes);

	// calculate layout and return number of crossings
	int CROSSINGS = ERLayoutAlgorithm(test, testAttributes);

	// draw graph to svg file
	GraphIO::SVGSettings settings;
	settings.fontSize(2);
	GraphIO::drawSVG(testAttributes, "C:\\Users\\Bart\\Desktop\\ERD.svg", settings);

	cout << "nr of crossings: " << CROSSINGS << endl;

	// calculate criteria value
	CriteriaTesting(test, testAttributes, CROSSINGS);

	// draw after bend promotion
	GraphIO::drawSVG(testAttributes, "C:\\Users\\Bart\\Desktop\\ERD2.svg", settings);
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

	//map<edge, string> relTypes;
	//map<edge, string>::iterator map_it2;

	// create all nodes
	for (size_t i = 0; i < js.size(); i++) {
		string name = js[i]["name"];
		node n = G.newNode();
		GA.label(n) = name;
		GA.fillColor(n) = Color::Name::Aquamarine;
		nodes.insert(pair<string, node>(name, n));
	}	
	
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

						/*
						edge ed = G.searchEdge(t, s);
						if (ed != 0) {
							map_it2 = relTypes.find(ed);
							cout << "edge: " << GA.label(s) << " -- " << GA.label(t) << " type: " << type << endl;
							cout << "edge: " << GA.label(t) << " -- " << GA.label(s) << " type: " << map_it2->second << endl << endl;
						}
						*/					

						// check for double edges and self-loops
						if (G.searchEdge(t, s) == 0 && GA.label(s) != GA.label(t)) {
							// make new edge
							edge e = G.newEdge(s, t);
							//relTypes.insert(pair<edge, string>(e, type));
							GA.strokeWidth(e) = 0.5;

							if (type == "UNI_TO_ONE") {
								GA.strokeType(e) = ogdf::StrokeType::Solid;
								GA.arrowType(e) = ogdf::EdgeArrow::None;
								GA.fillColor(s) = Color::Name::White;
								GA.fillColor(t) = Color::Name::White;
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
								GA.fillColor(s) = Color::Name::White;
								GA.fillColor(t) = Color::Name::White;
							}								
						}
					}					
				}
								
			}
		}
	}

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

int ERLayoutAlgorithm(Graph& G, GraphAttributes& GA) {
	int CROSSINGS = 0;

	PlanarizationLayout pl;
	SubgraphPlanarizer *crossMin = new SubgraphPlanarizer;

	PlanarSubgraphFast<int> *ps = new PlanarSubgraphFast<int>;
	VariableEmbeddingInserter *ves = new VariableEmbeddingInserter;

	crossMin->setSubgraph(ps);
	crossMin->setInserter(ves);

	OrthoLayout *ol = new OrthoLayout;
	//ol->separation(20.0);
	//ol->cOverhang(0.4);

	pl.setPlanarLayouter(ol);	

	pl.call(GA);

	CROSSINGS = pl.numberOfCrossings();
	
	return CROSSINGS;
}

void CriteriaTesting(Graph& G, GraphAttributes& GA, int CROSSINGS) {
	// test criteria
	double Nb = BendCriterion(G, GA);
	double Nc = CrossingCriterion(G, GA, CROSSINGS);
	double Nno = NodeOrthogonalityCriterion(G, GA);
	double Neo = EdgeOrthogonalityCriterion(G, GA);

	cout << "Criteria" << endl << endl;
	cout << "Crossings (N_c): " << Nc << endl;
	cout << "Bends (N_b): " << Nb << endl;
	cout << "Node Ortho (N_no): " << Nno << endl;
	cout << "Edge Ortho (N_eo): " << Neo << endl;
	cout << endl << endl;
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
 * calculate bend criterion
 */
double BendCriterion(Graph& G, GraphAttributes& GA) {
	// original amount of edges
	double m = G.edges.size();

	/*
	 * Maybe get this out of this function, because this is not the only criterion that uses bendpromotion
	 */
	BendPromotion(G, GA);

	// amount of edges after bend promotion
	double m2 = G.edges.size();

	double b_avg = (m2 - m) / m2;
	double Nb = 1 - b_avg;

	return Nb;
}

/*
 * calculate crossing criterion
 */
double CrossingCriterion(Graph& G, GraphAttributes& GA, double c) {
	// calculate all (theoretically) possible edge crossings
	double m2 = G.edges.size();
	double c_all = (m2 * (m2 - 1)) / 2;

	// calculate impossible edge crossings
	double c_impossible = 0;

	for (node n : G.nodes) {
		double sum = n->degree() * (n->degree() - 1);
		c_impossible += sum;
	}

	c_impossible = c_impossible / 2;

	// calculate real max possible crossings
	double c_max = c_all - c_impossible;

	// calculate crossings criterion
	double Nc = 1;
	if (c_max > 0) {
		Nc -= (c / c_max);
	}

	return Nc;
}

/*
 * calculate edge orthogonality criterion
 */
double EdgeOrthogonalityCriterion(Graph& G, GraphAttributes& GA) {
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
 * calculate node orthogonality criterion
 *
 * TODO: calculate value '2 * NODE_SIZE' from: 
 * GCD of the set of vertical and horizontal (pixel) differences between all geometrically adjacent nodes.
 */
double NodeOrthogonalityCriterion(Graph& G, GraphAttributes& GA) {
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

double AngleCriterion(Graph& G, GraphAttributes& GA) {
	double maxAngle = 0.0;
	List<edge> edges;
	for (node n : G.nodes) {
		maxAngle = 360 / n->degree();
		n->adjEdges(edges);
	}

	return maxAngle;
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