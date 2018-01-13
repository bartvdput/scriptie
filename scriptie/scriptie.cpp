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
#include <ogdf/graphalg/PageRank.h>

using namespace ogdf;
using namespace std;
using json = nlohmann::json;

void CreateGraph(Graph&, GraphAttributes&);
void Model1(Graph&, GraphAttributes&);
void Model2(Graph&, GraphAttributes&);
void Model3(Graph&, GraphAttributes&);
void Model4(Graph&, GraphAttributes&);
void Model6(Graph&, GraphAttributes&);
void Model7(Graph&, GraphAttributes&);
void Model8(Graph&, GraphAttributes&);
void Model9(Graph&, GraphAttributes&);
void Model10(Graph&, GraphAttributes&);

void BendPromotion(Graph&, GraphAttributes&);
void SetGraphLayout(Graph&, GraphAttributes&);
double BendCriterion(Graph&, GraphAttributes&);
double CrossingCriterion(Graph&, double);
double UniformEdgeCriterion(Graph&, GraphAttributes&);
double NodeOrthogonalityCriterion(Graph&, GraphAttributes&);
double EdgeOrthogonalityCriterion(Graph&, GraphAttributes&);
void CriteriaTesting(Graph&, GraphAttributes&, int);
void addRelations(Graph&, GraphAttributes&);
int ERLayoutAlgorithm(Graph&, GraphAttributes&);
int getBiconnectedComponents(Graph&);
void CreateGraphFromJson(Graph&, GraphAttributes&, string);
void FindImportantNodes(Graph&, GraphAttributes&);
void findSimplePaths(Graph&, GraphAttributes&, node, node, int);
void DFSUtil(Graph&, GraphAttributes&, node, node, int, bool[], int&, node[], Graph&, bool[]);
void GenerateSubgraph(Graph&, GraphAttributes&, int, Graph&, GraphAttributes&);
double edgeLength(edge&, GraphAttributes&);

const double NODE_WIDTH = 40.0;
const double NODE_HEIGHT = 20.0;
const float STROKEWIDTH = 0.4f;
const double PI = 3.141592653589793238463;
const Color MAJOR_NODES = Color::Name::Blue;
const Color SECONDARY_NODES = Color::Name::Lightblue;
const Color MINOR_NODES = Color::Name::White;
const Color BEND_COLOR = Color::Name::Black;
const Color PINK_COLOR = Color::Name::Pink;

int main() {
	Graph G;
	GraphAttributes GA(G, GraphAttributes::nodeType | GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::nodeLabel | GraphAttributes::nodeStyle | GraphAttributes::edgeType | GraphAttributes::edgeArrow | GraphAttributes::edgeStyle | GraphAttributes::edgeLabel);
	
	Graph subGraph;
	GraphAttributes subGA(subGraph, GraphAttributes::nodeType | GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::nodeLabel | GraphAttributes::nodeStyle | GraphAttributes::edgeType | GraphAttributes::edgeArrow | GraphAttributes::edgeStyle | GraphAttributes::edgeLabel);

	GraphIO::SVGSettings settings;
	settings.fontSize(4);

	// read file
	string file = "entities-big_system.json";
	//cout << "Enter file name to read in: " << endl;	
	//cin >> file;
	//CreateGraphFromJson(G, GA, file);
	Model4(G, GA);

	int maxIndex = 0;
	for (node n : G.nodes) {
		if (n->index() > maxIndex)
			maxIndex = n->index();
	}

	// set some layout properties
	SetGraphLayout(G, GA);
	SetGraphLayout(subGraph, subGA);

	//FindImportantNodes(G, GA);
	//GenerateSubgraph(G, GA, maxIndex, subGraph, subGA);

	// calculate layout and return number of crossings
	int CROSSINGS = ERLayoutAlgorithm(G, GA);

	// draw graph to svg file
	GraphIO::drawSVG(GA, "C:\\Users\\Bart\\Desktop\\model4_2.svg", settings);

	cout << "nr of crossings: " << CROSSINGS << endl;

	// calculate criteria value
	CriteriaTesting(G, GA, CROSSINGS);

	/*
	// draw graph to svg file
	GraphIO::drawSVG(GA, "C:\\Users\\Bart\\Desktop\\ERD3.svg", settings);

	cout << "----------------------------" << endl << endl << endl;

	// calculate layout and return number of crossings
	int CROSSINGS_subgraph = ERLayoutAlgorithm(subGraph, subGA);

	// draw graph to svg file
	GraphIO::drawSVG(subGA, "C:\\Users\\Bart\\Desktop\\ERD2.svg", settings);

	cout << "nr of crossings: " << CROSSINGS_subgraph << endl;

	// calculate criteria value
	CriteriaTesting(subGraph, subGA, CROSSINGS_subgraph);
	*/	
	
	return 0;
}

void CreateGraphFromJson(Graph& G, GraphAttributes& GA, string file) {
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

						// check for double edges and self-loops
						if (G.searchEdge(t, s) == 0 && GA.label(s) != GA.label(t)) {
							// make new edge
							edge e = G.newEdge(s, t);
							//relTypes.insert(pair<edge, string>(e, type));
							GA.strokeWidth(e) = 0.5;

							if (type == "UNI_TO_ONE") {
								GA.strokeType(e) = ogdf::StrokeType::Solid;
								GA.arrowType(e) = ogdf::EdgeArrow::None;
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

	// check degree and delete non-connected nodes
	for (map_it = nodes.begin(); map_it != nodes.end(); map_it++) {
		node n = map_it->second;
		if (n->degree() == 0) {
			G.delNode(n);
		}
	}	
}

/*
 * pagerank algorithm to find major/minor nodes
 */
void FindImportantNodes(Graph& G, GraphAttributes& GA) {
	EdgeArray<double> edgeWeight = EdgeArray<double>(G);
	NodeArray<double> pageRank = NodeArray<double>(G);
	NodeArray<double>::iterator it;

	BasicPageRank pr;
	pr.call(G, edgeWeight, pageRank);

	for (it = pageRank.begin(); it != pageRank.end(); it++) {
		double rank = it.value();
		node n = it.key();
		if (rank > 0.5) {
			GA.fillColor(n) = MAJOR_NODES;
		}
		else if (rank > 0.2) {
			GA.fillColor(n) = SECONDARY_NODES;
		}
		//cout << GA.label(it.key()) << ": " << it.value() << endl;
	}
}

/*
 * the actual layout algorithm
 */
int ERLayoutAlgorithm(Graph& G, GraphAttributes& GA) {
	cout << "Node count: " << G.nodes.size() << endl;
	cout << "Edge count: " << G.edges.size() << endl;
	int CROSSINGS = 0;

	PlanarizationLayout pl;
	SubgraphPlanarizer *crossMin = new SubgraphPlanarizer;

	PlanarSubgraphFast<int> *ps = new PlanarSubgraphFast<int>;
	//PlanarSubgraphBoyerMyrvold *ps = new PlanarSubgraphBoyerMyrvold;
	VariableEmbeddingInserter *ves = new VariableEmbeddingInserter;

	crossMin->setSubgraph(ps);
	crossMin->setInserter(ves);

	OrthoLayout *ol = new OrthoLayout;
	ol->separation(NODE_WIDTH);
	ol->cOverhang(1.0);

	pl.setPlanarLayouter(ol);	

	pl.call(GA);

	CROSSINGS = pl.numberOfCrossings();
	
	return CROSSINGS;
}

void CriteriaTesting(Graph& G, GraphAttributes& GA, int CROSSINGS) {
	double Nue = UniformEdgeCriterion(G, GA);
	double Nno_nodes = NodeOrthogonalityCriterion(G, GA);
	double Nb = BendCriterion(G, GA);
	double Nc = CrossingCriterion(G, CROSSINGS);
	double Neo = EdgeOrthogonalityCriterion(G, GA);	

	cout << "Criteria" << endl << endl << endl;
	cout << "Crossings (N_c): " << Nc << endl;
	cout << "Bends (N_b): " << Nb << endl;
	cout << "Node Ortho (N_no): " << Nno_nodes << endl;
	cout << "Edge Ortho (N_eo): " << Neo << endl;
	cout << "Uniform edge lengths (N_ue): " << Nue << endl;
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

void BendPromotion(Graph& G, GraphAttributes& GA) {
	List<edge> edges;
	G.allEdges(edges);

	while (!edges.empty()) {
		edge e = edges.popFrontRet();
		DPolyline bends_e = GA.bends(e);		
		node s = e->source();
		node t = e->target();

		//check if an edge has bendpoints
		while (!bends_e.empty()) {
			DPoint p = bends_e.front();

			//insert new node
			node n = G.newNode();
			GA.width(n) = 1;
			GA.height(n) = 1;
			GA.fillColor(n) = BEND_COLOR;			
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

/*
 * calculate bend criterion
 */
double BendCriterion(Graph& G, GraphAttributes& GA) {
	// original amount of edges
	double m = G.edges.size();

	BendPromotion(G, GA);

	// amount of edges after bend promotion
	double m2 = G.edges.size() - (m*2); // - (m*2) because bend points are created where edges connect to nodes.
	double b = (m2 - m);
	double b_avg = b / m2;

	return 1 - b_avg;
}

/*
 * calculate crossing criterion
 */
double CrossingCriterion(Graph& G, double c) {
	// calculate all (theoretically) possible edge crossings
	double m = G.edges.size();
	double c_all = (m * (m - 1)) / 2;

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
	} else Nc = 0;

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

		if (GA.fillColor(s) == BEND_COLOR && GA.fillColor(t) == BEND_COLOR) {
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
	}
	
	double avg = ((1.0 / G.edges.size()) * sum);
	return  1.0 - avg;
}

/*
 * calculate node orthogonality criterion
 */
double NodeOrthogonalityCriterion(Graph& G, GraphAttributes& GA) {
	List<int> differencesX;
	List<int>::iterator itX;
	int countX = 0;

	List<int> differencesY;
	List<int>::iterator itY;
	int countY = 0;

	// find all distances between nodes
	for (node n : G.nodes) {
		for (node m : G.nodes) {
			if (n != m) {
				if (GA.x(n) == GA.x(m)) {
					int dif = abs(GA.y(n) - GA.y(m));
					if (dif != 0) {
						differencesY.pushBack(dif);
						countY++;
					}
						
				} else if (GA.y(n) == GA.y(m)) {
					int dif = abs(GA.x(n) - GA.x(m));
					if (dif != 0) {
						differencesX.pushBack(dif);
						countX++;
					}						
				}
			}								
		}					
	}

	// use distances to find gcd of distance between nodes
	itX = differencesX.get(0);
	int gcdResultX = *itX;
	for (int i = 1; i<countX; i++) {
		itX = differencesX.get(i);
		gcdResultX = Math::gcd(gcdResultX, *itX);
	}

	int gcdResultY = 1;
	if (differencesY.size() > 0) {
		itY = differencesY.get(0);
		gcdResultY = *itY;
		for (int i = 1; i<countY; i++) {
			itY = differencesY.get(i);
			gcdResultY = Math::gcd(gcdResultY, *itY);
		}
	}

	cout << "gcdResultX: " << gcdResultX << endl;
	cout << "gcdResultY: " << gcdResultY << endl;
	double minX = DBL_MAX, minY = DBL_MAX, maxX = 0, maxY = 0;

	for (node n : G.nodes) {
		if (GA.x(n) < minX)
			minX = GA.x(n);

		if (GA.x(n) > maxX)
			maxX = GA.x(n);

		if (GA.y(n) < minY)
			minY = GA.y(n);

		if (GA.y(n) > maxY)
			maxY = GA.y(n);
	}

	cout << "min: " << minX << ", " << minY << endl;
	cout << "max: " << maxX << ", " << maxY << endl;

	double width = (maxX - minX) / gcdResultX + 1;
	double height = (maxY - minY)/ gcdResultY + 1;

	// calculate actual criterion value
	cout << "width: " << width << endl;
	cout << "height: " << height << endl;
	double A = (width)*(height);

	double Nno = G.nodes.size() / A;

	return Nno;
}

double UniformEdgeCriterion(Graph& G, GraphAttributes& GA) {
	double total_length = 0;
	for (edge e : G.edges) {
		total_length += edgeLength(e, GA);
	}
	double avg_length = total_length / G.edges.size();

	double total_deviation = 0;
	for (edge e : G.edges) {
		total_deviation += abs(edgeLength(e, GA) - avg_length);
	}
	double avg_deviation = total_deviation / G.edges.size();

	cout << "l " << avg_length << endl;
	cout << "d " << avg_deviation << endl;
	
	if (avg_deviation < avg_length)
		return 1 - (avg_deviation / avg_length);
	else {
		cout << endl << "oi m8" << endl << endl;
		return 0;
	}
}

double edgeLength(edge& e, GraphAttributes& GA) {
	double length = 0;

	DPolyline bends_e = GA.bends(e);
	DPoint s = bends_e.popFrontRet();
	DPoint t = bends_e.popBackRet();

	//check if an edge has bendpoints
	while (!bends_e.empty()) {
		DPoint p = bends_e.popFrontRet();

		if (s.m_x == p.m_x)
			length += abs(s.m_y - p.m_y);
		else
			length += abs(s.m_x - p.m_x);
		
		s = p;
	}

	if (s.m_x == t.m_x)
		length += abs(s.m_y - t.m_y);
	else
		length += abs(s.m_x - t.m_x);

	return length;
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
	for (edge e : G.edges) {
		GA.strokeWidth(e) = STROKEWIDTH;
		GA.strokeColor(e) = Color::Name::Black;
		GA.strokeType(e) = ogdf::StrokeType::Solid;
		GA.arrowType(e) = ogdf::EdgeArrow::None;
	}

	for (node n : G.nodes) {
		GA.fillColor(n) = Color::Name::Lightgray;
		//GA.strokeWidth(n) = STROKEWIDTH;
	}
}

/*
 * find all simple paths between two nodes
 */
void findSimplePaths(Graph& G, GraphAttributes& GA, node n, node t, int nodeCount) {
	Graph subGraph;
	GraphAttributes subGA(subGraph, GraphAttributes::nodeType | GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::nodeLabel | GraphAttributes::nodeStyle | GraphAttributes::edgeType | GraphAttributes::edgeArrow | GraphAttributes::edgeStyle | GraphAttributes::edgeLabel);
	
	// keep track of visited nodes
	bool *visited = new bool[nodeCount];
	bool *nodeInNewGraph = new bool[nodeCount];

	// keep track of current path
	int *path = new int[nodeCount];
	node *nodePath = new node[nodeCount];
	int pathIndex = 0;
	
	// initialize all nodes as NOT visited AND NOT in new graph
	for (int i = 0; i < nodeCount; i++) {
		visited[i] = false;
		nodeInNewGraph[i] = false;
	}

	int index = n->index();

	// call helper function
	DFSUtil(G, GA, n, t, index, visited, pathIndex, nodePath, subGraph, nodeInNewGraph);

	cout << endl;

	ERLayoutAlgorithm(subGraph, subGA);
	// draw graph to svg file
	GraphIO::SVGSettings settings;
	settings.fontSize(2);
	GraphIO::drawSVG(subGA, "C:\\Users\\Bart\\Desktop\\ERD2.svg", settings);
}

/*
 * helper function for finding all simple paths
 */
void DFSUtil(Graph& G, GraphAttributes& GA, node n, node t, int index, bool visited[], int &pathIndex, node nodePath[], Graph& subGraph, bool nodeInNewGraph[]) {
	// set node to visited and add node to path
	visited[index] = true;
	nodePath[pathIndex] = n;
	pathIndex++;

	List<edge> adj;
	n->adjEdges(adj);
	node m;

	// check with all adjacent nodes if it is the target node, else recursively call the helper function on an adjacent node
	for (edge e : adj) {
		if (e->target() != n)
			m = e->target();
		else m = e->source();

		int newIndex = m->index();

		if (m == t) {
			for (int i = 0; i < pathIndex; i++) {				
				int tempIndex = nodePath[i]->index();
				if (!nodeInNewGraph[tempIndex]) {
					node newNode = subGraph.newNode();
					nodeInNewGraph[tempIndex] = true;
				}				
			}				
		} else if (!visited[newIndex] && !(GA.fillColor(m) == MINOR_NODES)){
			DFSUtil(G, GA, m, t, newIndex, visited, pathIndex, nodePath, subGraph, nodeInNewGraph);
		}
	}

	// remove current node from path and set to not visited
	pathIndex--;
	visited[index] = false;
}

void GenerateSubgraph(Graph& G, GraphAttributes& GA, int graphSize, Graph& subGraph, GraphAttributes& subGA) {
	node *nodes = new node[graphSize];

	for (node n : G.nodes) {
		if (GA.fillColor(n) == MAJOR_NODES || GA.fillColor(n) == SECONDARY_NODES) {
			node m = subGraph.newNode();
			subGA.label(m) = GA.label(n);
			subGA.fillColor(m) = GA.fillColor(n);
			nodes[n->index()] = m;
		}
	}

	for (edge e : G.edges) {
		node s = e->source();
		node t = e->target();
		if ((GA.fillColor(s) == MAJOR_NODES || GA.fillColor(s) == SECONDARY_NODES) && (GA.fillColor(t) == MAJOR_NODES || GA.fillColor(t) == SECONDARY_NODES)) {
			int sIndex = s->index();
			int tIndex = t->index();
			if (sIndex < graphSize && tIndex < graphSize) {
				edge w = subGraph.newEdge(nodes[sIndex], nodes[tIndex]);
			} else {
				cout << "wot? " << sIndex << " " << tIndex << endl;
			}			
		}
	}
}

// Model 1
void Model1(Graph& graph, GraphAttributes& GA) {
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
}

// Model 2
void Model2(Graph& graph, GraphAttributes& GA) {
	node Ref_Med_Type = graph.newNode();
	node Medications = graph.newNode();
	node Doctors = graph.newNode();
	node Addresses = graph.newNode();
	node Customers = graph.newNode();
	node Pharma_Companies = graph.newNode();
	node Prescriptions = graph.newNode();
	node Ref_Payment_Methods = graph.newNode();
	node Ref_Prescription_Status = graph.newNode();
	node Brand_Name_Med = graph.newNode();
	node Generic_Med = graph.newNode();
	node Gen_to_Brand_Name_Correspondence = graph.newNode();
	node Prescription_Items = graph.newNode();
	node Items_Ordered = graph.newNode();
	

	GA.label(Ref_Med_Type) = "Ref_Med_Type";
	GA.label(Medications) = "Medications";
	GA.label(Doctors) = "Doctors";
	GA.label(Addresses) = "Addresses";
	GA.label(Customers) = "Customers";
	GA.label(Pharma_Companies) = "Pharma_Companies";
	GA.label(Prescriptions) = "Prescriptions";
	GA.label(Ref_Payment_Methods) = "Ref_Payment_Methods";
	GA.label(Ref_Prescription_Status) = "Ref_Prescription_Status";
	GA.label(Brand_Name_Med) = "Brand_Name_Med";
	GA.label(Generic_Med) = "Generic_Med";
	GA.label(Gen_to_Brand_Name_Correspondence) = "Gen_to_Brand_Name_Correspondence";
	GA.label(Prescription_Items) = "Prescription_Items";
	GA.label(Items_Ordered) = "Items_Ordered";

	edge e1 = graph.newEdge(Ref_Med_Type, Medications);
	edge e2 = graph.newEdge(Pharma_Companies, Brand_Name_Med);
	edge e3 = graph.newEdge(Medications, Brand_Name_Med);
	edge e4 = graph.newEdge(Medications, Generic_Med);
	edge e5 = graph.newEdge(Brand_Name_Med, Gen_to_Brand_Name_Correspondence);
	edge e6 = graph.newEdge(Generic_Med, Gen_to_Brand_Name_Correspondence);
	edge e7 = graph.newEdge(Medications, Prescription_Items);
	edge e8 = graph.newEdge(Doctors, Addresses);
	edge e9 = graph.newEdge(Doctors, Prescriptions);
	edge e10 = graph.newEdge(Prescriptions, Prescription_Items);
	edge e11 = graph.newEdge(Prescription_Items, Items_Ordered);
	edge e12 = graph.newEdge(Addresses, Customers);
	edge e13 = graph.newEdge(Customers, Prescriptions);
	edge e14 = graph.newEdge(Prescriptions, Ref_Payment_Methods);
	edge e15 = graph.newEdge(Prescriptions, Ref_Prescription_Status);	
}

// Model 3
void Model3(Graph& graph, GraphAttributes& GA) {
	node Customers = graph.newNode();
	node Car_Features = graph.newNode();
	node Cars_for_Sale = graph.newNode();
	node Car_Models = graph.newNode();
	node Addresses = graph.newNode();
	node Customer_Preferences = graph.newNode();
	node Car_Manufacturers = graph.newNode();
	node features_on_Cars_for_Sale = graph.newNode();
	node Vehicle_Categories = graph.newNode();
	node Payment_Status = graph.newNode();
	node Cars_Sold = graph.newNode();
	node Customer_Payments = graph.newNode();
	node Car_Loans = graph.newNode();
	node Insurance_Policies = graph.newNode();
	node Finance_Companies = graph.newNode();
	node Insurance_Companies = graph.newNode();

	GA.label(Customers) = "Customers";
	GA.label(Car_Features) = "Car_Features";
	GA.label(Cars_for_Sale) = "Cars_for_Sale";
	GA.label(Car_Models) = "Car_Models";
	GA.label(Addresses) = "Addresses";
	GA.label(Customer_Preferences) = "Customer_Preferences";
	GA.label(Car_Manufacturers) = "Car_Manufacturers";
	GA.label(features_on_Cars_for_Sale) = "features_on_Cars_for_Sale";
	GA.label(Vehicle_Categories) = "Vehicle_Categories";
	GA.label(Payment_Status) = "Payment_Status";
	GA.label(Cars_Sold) = "Cars_Sold";
	GA.label(Customer_Payments) = "Customer_Payments";
	GA.label(Car_Loans) = "Car_Loans";
	GA.label(Insurance_Policies) = "Insurance_Policies";
	GA.label(Finance_Companies) = "Finance_Companies";
	GA.label(Insurance_Companies) = "Insurance_Companies";

	edge e1 = graph.newEdge(Customers, Addresses);
	edge e2 = graph.newEdge(Customers, Customer_Payments);
	edge e3 = graph.newEdge(Customers, Cars_Sold);
	edge e4 = graph.newEdge(Customers, Customer_Preferences);
	edge e5 = graph.newEdge(Car_Features, Customer_Preferences);
	edge e6 = graph.newEdge(Car_Features, features_on_Cars_for_Sale);
	edge e7 = graph.newEdge(Cars_for_Sale, features_on_Cars_for_Sale);
	edge e8 = graph.newEdge(Cars_for_Sale, Cars_Sold);
	edge e9 = graph.newEdge(Cars_for_Sale, Car_Models);
	edge e10 = graph.newEdge(Cars_for_Sale, Car_Manufacturers);
	edge e11 = graph.newEdge(Cars_for_Sale, Vehicle_Categories);
	edge e12 = graph.newEdge(Payment_Status, Customers);
	edge e13 = graph.newEdge(Customers, Customer_Payments);
	edge e14 = graph.newEdge(Cars_Sold, Customer_Payments);
	edge e15 = graph.newEdge(Cars_Sold, Car_Loans);
	edge e16 = graph.newEdge(Cars_Sold, Insurance_Policies);
	edge e17 = graph.newEdge(Car_Loans, Finance_Companies);
	edge e18 = graph.newEdge(Insurance_Policies, Insurance_Companies);
}

// Model 4
void Model4(Graph& graph, GraphAttributes& GA) {
	// Multiple
	node Addresses = graph.newNode();
	node Students = graph.newNode();
	GA.label(Addresses) = "Addresses";
	GA.label(Students) = "Students";

	// 4.1
	/*
	node Teachers = graph.newNode();	
	node Assessment_Notes = graph.newNode();
	node Detention = graph.newNode();
	node Students_in_Detention = graph.newNode();
	node Behaviour_Incident = graph.newNode();

	GA.label(Teachers) = "Teachers";	
	GA.label(Assessment_Notes) = "Assessment_Notes";
	GA.label(Detention) = "Medications";
	GA.label(Students_in_Detention) = "Students_in_Detention";
	GA.label(Behaviour_Incident) = "Behaviour_Incident";

	edge e1 = graph.newEdge(Teachers, Addresses);
	edge e2 = graph.newEdge(Teachers, Detention);
	edge e3 = graph.newEdge(Teachers, Assessment_Notes);
	edge e4 = graph.newEdge(Students, Assessment_Notes);
	edge e5 = graph.newEdge(Students, Students_in_Detention);
	edge e6 = graph.newEdge(Students, Behaviour_Incident);
	edge e7 = graph.newEdge(Detention, Students_in_Detention);
	edge e8 = graph.newEdge(Students_in_Detention, Behaviour_Incident);
	edge e9 = graph.newEdge(Students, Addresses);
	*/
	
	//4.2	
	node People = graph.newNode();
	node Ref_Gender = graph.newNode();
	node Ref_Subjects = graph.newNode();
	node Classes = graph.newNode();
	node Instructors = graph.newNode();
	node Ref_Retest_Status_Codes = graph.newNode();
	node Student_Classes = graph.newNode();
	node Instructors_Classes = graph.newNode();
	node Ref_Score_Status = graph.newNode();

	GA.label(People) = "People";
	GA.label(Ref_Gender) = "Ref_Gender";
	GA.label(Ref_Subjects) = "Ref_Subjects";
	GA.label(Classes) = "Classes";
	GA.label(Instructors) = "Instructors";
	GA.label(Ref_Retest_Status_Codes) = "Ref_Retest_Status_Codes";
	GA.label(Student_Classes) = "Student_Classes";
	GA.label(Instructors_Classes) = "Instructors_Classes";
	GA.label(Ref_Score_Status) = "Ref_Score_Status";

	edge ee1 = graph.newEdge(People, Addresses);
	edge ee2 = graph.newEdge(People, Ref_Gender);
	edge ee3 = graph.newEdge(People, Students);
	edge ee4 = graph.newEdge(People, Instructors);
	edge ee5 = graph.newEdge(Students, Student_Classes);
	edge ee6 = graph.newEdge(Students, Ref_Subjects);
	edge ee7 = graph.newEdge(Ref_Subjects, Classes);
	edge ee8 = graph.newEdge(Classes, Instructors_Classes);
	edge ee9 = graph.newEdge(Classes, Student_Classes);
	edge ee10 = graph.newEdge(Student_Classes, Ref_Retest_Status_Codes);
	edge ee11 = graph.newEdge(Instructors, Instructors_Classes);
	edge ee12 = graph.newEdge(Student_Classes, Ref_Score_Status);
	
	/*	
	//4.3	
	node Ref_Academic_Years = graph.newNode();
	node Ref_Salutations = graph.newNode();
	node Ref_Staff_Roles = graph.newNode();
	node Year_Groups = graph.newNode();
	node Activities = graph.newNode();
	node Staff = graph.newNode();
	node Forms = graph.newNode();
	node Scheduled_Activities = graph.newNode();
	node Activities_run_by_Staff = graph.newNode();
	node Student_Addresses = graph.newNode();
	node Ref_Attainment_Levels = graph.newNode();
	node Ref_Payment_Status = graph.newNode();
	node Student_Activities = graph.newNode();

	GA.label(Ref_Academic_Years) = "Ref_Academic_Years";
	GA.label(Ref_Salutations) = "Ref_Salutations";
	GA.label(Ref_Staff_Roles) = "Ref_Staff_Roles";
	GA.label(Year_Groups) = "Year_Groups";
	GA.label(Activities) = "Activities";
	GA.label(Staff) = "Staff";
	GA.label(Forms) = "Forms";
	GA.label(Scheduled_Activities) = "Scheduled_Activities";
	GA.label(Activities_run_by_Staff) = "Activities_run_by_Staff";
	GA.label(Student_Addresses) = "Student_Addresses";
	GA.label(Ref_Attainment_Levels) = "Ref_Attainment_Levels";
	GA.label(Ref_Payment_Status) = "Ref_Payment_Status";
	GA.label(Student_Activities) = "Student_Activities";

	edge eee1 = graph.newEdge(Ref_Academic_Years, Scheduled_Activities);
	edge eee2 = graph.newEdge(Activities, Scheduled_Activities);
	edge eee3 = graph.newEdge(Ref_Salutations, Staff);
	edge eee4 = graph.newEdge(Ref_Staff_Roles, Staff);
	edge eee5 = graph.newEdge(Year_Groups, Staff);
	edge eee6 = graph.newEdge(Forms, Staff);
	edge eee7 = graph.newEdge(Forms, Students);
	edge eee8 = graph.newEdge(Scheduled_Activities, Staff);
	edge eee9 = graph.newEdge(Addresses, Student_Addresses);
	edge eee10 = graph.newEdge(Students, Student_Addresses);
	edge eee11 = graph.newEdge(Scheduled_Activities, Activities_run_by_Staff);
	edge eee12 = graph.newEdge(Activities_run_by_Staff, Staff);
	edge eee13 = graph.newEdge(Students, Student_Activities);
	edge eee14 = graph.newEdge(Scheduled_Activities, Student_Activities);
	edge eee15 = graph.newEdge(Year_Groups, Forms);
	edge eee16 = graph.newEdge(Ref_Attainment_Levels, Student_Activities);
	edge eee17 = graph.newEdge(Ref_Payment_Status, Student_Activities);	
	*/
}

// Model 6
void Model6(Graph& graph, GraphAttributes& GA) {
	node Authors = graph.newNode();
	GA.label(Authors) = "Authors";
	node BookAuthor = graph.newNode();
	GA.label(BookAuthor) = "BookAuthor";
	node Books = graph.newNode();
	GA.label(Books) = "Books";
	node Employees = graph.newNode();
	GA.label(Employees) = "Employees";
	node Jobs = graph.newNode();
	GA.label(Jobs) = "Jobs";
	node Publishers = graph.newNode();
	GA.label(Publishers) = "Publishers";

	edge e0 = graph.newEdge(Jobs, Employees);
	edge e1 = graph.newEdge(Publishers, Books);
	edge e2 = graph.newEdge(Books, BookAuthor);
	edge e3 = graph.newEdge(Publishers, Employees);
	edge e4 = graph.newEdge(Authors, BookAuthor);
}

// Model 7
void Model7(Graph& graph, GraphAttributes& GA) {
	node Customers = graph.newNode();
	GA.label(Customers) = "Customers";
	node Employee = graph.newNode();
	GA.label(Employee) = "Employee";
	node Invoice = graph.newNode();
	GA.label(Invoice) = "Invoice";
	node Invoice_details = graph.newNode();
	GA.label(Invoice_details) = "Invoice details";
	node Parts = graph.newNode();
	GA.label(Parts) = "Parts";
	node Suppliers = graph.newNode();
	GA.label(Suppliers) = "Suppliers";
	node Vehicles = graph.newNode();
	GA.label(Vehicles) = "Vehicles";

	edge e0 = graph.newEdge(Customers, Vehicles);
	edge e1 = graph.newEdge(Employee, Invoice);
	edge e2 = graph.newEdge(Invoice, Invoice_details);
	edge e3 = graph.newEdge(Parts, Invoice_details);
	edge e4 = graph.newEdge(Suppliers, Parts);
	edge e5 = graph.newEdge(Vehicles, Invoice);
}

void Model8(Graph& graph, GraphAttributes& GA) {
	node Choice = graph.newNode();
	GA.label(Choice) = "Choice";
	node Choice_MultipleChoiceResponse = graph.newNode();
	GA.label(Choice_MultipleChoiceResponse) = "Choice_MultipleChoiceResponse";
	node DateResponse = graph.newNode();
	GA.label(DateResponse) = "DateResponse";
	node MultipleChoiceResponse = graph.newNode();
	GA.label(MultipleChoiceResponse) = "MultipleChoiceResponse";
	node NumericResponse = graph.newNode();
	GA.label(NumericResponse) = "NumericResponse";
	node Person = graph.newNode();
	GA.label(Person) = "Person";
	node Question = graph.newNode();
	GA.label(Question) = "Question";
	node QuestionType = graph.newNode();
	GA.label(QuestionType) = "QuestionType";
	node Reponse = graph.newNode();
	GA.label(Reponse) = "Reponse";
	node SingleChoiceResponse = graph.newNode();
	GA.label(SingleChoiceResponse) = "SingleChoiceResponse";
	node Survey = graph.newNode();
	GA.label(Survey) = "Survey";
	node SurveyParticipant = graph.newNode();
	GA.label(SurveyParticipant) = "SurveyParticipant";
	node TextResponse = graph.newNode();
	GA.label(TextResponse) = "TextResponse";
	node YesNoResponse = graph.newNode();
	GA.label(YesNoResponse) = "YesNoResponse";

	edge e0 = graph.newEdge(Choice, Choice_MultipleChoiceResponse);
	edge e1 = graph.newEdge(Choice, SingleChoiceResponse);
	edge e2 = graph.newEdge(MultipleChoiceResponse, Choice_MultipleChoiceResponse);
	edge e3 = graph.newEdge(Person, SurveyParticipant);
	edge e4 = graph.newEdge(QuestionType, Question);
	edge e5 = graph.newEdge(Question, Choice);
	edge e6 = graph.newEdge(Question, DateResponse);
	edge e7 = graph.newEdge(Question, MultipleChoiceResponse);
	edge e8 = graph.newEdge(Question, NumericResponse);
	edge e9 = graph.newEdge(Question, SingleChoiceResponse);
	edge e10 = graph.newEdge(Question, TextResponse);
	edge e11 = graph.newEdge(Question, YesNoResponse);
	edge e12 = graph.newEdge(Reponse, DateResponse);
	edge e13 = graph.newEdge(Reponse, MultipleChoiceResponse);
	edge e14 = graph.newEdge(Reponse, NumericResponse);
	edge e15 = graph.newEdge(Reponse, SingleChoiceResponse);
	edge e16 = graph.newEdge(Reponse, TextResponse);
	edge e17 = graph.newEdge(Reponse, YesNoResponse);
	edge e18 = graph.newEdge(SurveyParticipant, Reponse);
	edge e19 = graph.newEdge(Survey, Question);
	edge e20 = graph.newEdge(Survey, SurveyParticipant);
}

void Model9(Graph& graph, GraphAttributes& GA) {
	node Address = graph.newNode();
	GA.label(Address) = "Address";
	node Basket = graph.newNode();
	GA.label(Basket) = "Basket";
	node BasketProducts = graph.newNode();
	GA.label(BasketProducts) = "BasketProducts";
	node BasketStatus = graph.newNode();
	GA.label(BasketStatus) = "BasketStatus";
	node Category = graph.newNode();
	GA.label(Category) = "Category";
	node CategoryProduct = graph.newNode();
	GA.label(CategoryProduct) = "CategoryProduct";
	node City = graph.newNode();
	GA.label(City) = "City";
	node Country = graph.newNode();
	GA.label(Country) = "Country";
	node Customer = graph.newNode();
	GA.label(Customer) = "Customer";
	node Manufacturer = graph.newNode();
	GA.label(Manufacturer) = "Manufacturer";
	node Order = graph.newNode();
	GA.label(Order) = "Order";
	node OrderedProduct = graph.newNode();
	GA.label(OrderedProduct) = "OrderedProduct";
	node OrderStatus = graph.newNode();
	GA.label(OrderStatus) = "OrderStatus";
	node Product = graph.newNode();
	GA.label(Product) = "Product";
	node ProductImage = graph.newNode();
	GA.label(ProductImage) = "ProductImage";
	node ProductSpecial = graph.newNode();
	GA.label(ProductSpecial) = "ProductSpecial";
	node Review = graph.newNode();
	GA.label(Review) = "Review";
	node StateOrProvince = graph.newNode();
	GA.label(StateOrProvince) = "StateOrProvince";
	node TaxClass = graph.newNode();
	GA.label(TaxClass) = "TaxClass";
	node TaxRate = graph.newNode();
	GA.label(TaxRate) = "TaxRate";
	node ZipCode = graph.newNode();
	GA.label(ZipCode) = "ZipCode";

	edge e1 = graph.newEdge(Address, Order);
	edge e2 = graph.newEdge(BasketStatus, Basket);
	edge e3 = graph.newEdge(Basket, BasketProducts);
	//edge e4 = graph.newEdge(Category, Category);
	edge e5 = graph.newEdge(Category, CategoryProduct);
	edge e6 = graph.newEdge(City, ZipCode);
	edge e7 = graph.newEdge(Country, ZipCode);
	edge e8 = graph.newEdge(Customer, Address);
	edge e9 = graph.newEdge(Customer, Basket);
	edge e10 = graph.newEdge(Customer, Order);
	edge e11 = graph.newEdge(Customer, Review);
	edge e12 = graph.newEdge(Manufacturer, Product);
	edge e13 = graph.newEdge(OrderStatus, Order);
	edge e14 = graph.newEdge(Order, OrderedProduct);
	edge e15 = graph.newEdge(Product, BasketProducts);
	edge e16 = graph.newEdge(Product, CategoryProduct);
	edge e17 = graph.newEdge(Product, OrderedProduct);
	edge e18 = graph.newEdge(Product, ProductImage);
	edge e19 = graph.newEdge(Product, ProductSpecial);
	edge e20 = graph.newEdge(Product, Review);
	edge e21 = graph.newEdge(StateOrProvince, ZipCode);
	edge e22 = graph.newEdge(TaxClass, Product);
	edge e23 = graph.newEdge(TaxClass, TaxRate);
	edge e24 = graph.newEdge(ZipCode, Address);
}

void Model10(Graph& graph, GraphAttributes& GA) {
	node Address = graph.newNode();
	GA.label(Address) = "Address";
	node Agent = graph.newNode();
	GA.label(Agent) = "Agent";
	node City = graph.newNode();
	GA.label(City) = "City";
	node Client = graph.newNode();
	GA.label(Client) = "Client";
	node ClientWish = graph.newNode();
	GA.label(ClientWish) = "ClientWish";
	node Country = graph.newNode();
	GA.label(Country) = "Country";
	node Feature = graph.newNode();
	GA.label(Feature) = "Feature";
	node Owner = graph.newNode();
	GA.label(Owner) = "Owner";
	node Person = graph.newNode();
	GA.label(Person) = "Person";
	node Property = graph.newNode();
	GA.label(Property) = "Property";
	node PropertyFeature = graph.newNode();
	GA.label(PropertyFeature) = "PropertyFeature";
	node PropertyForRent = graph.newNode();
	GA.label(PropertyForRent) = "PropertyForRent";
	node PropertyForSale = graph.newNode();
	GA.label(PropertyForSale) = "PropertyForSale";
	node Rent = graph.newNode();
	GA.label(Rent) = "Rent";
	node RentalAgreement = graph.newNode();
	GA.label(RentalAgreement) = "RentalAgreement";
	node Sale = graph.newNode();
	GA.label(Sale) = "Sale";
	node StateOrProvince = graph.newNode();
	GA.label(StateOrProvince) = "StateOrProvince";
	node Visit = graph.newNode();
	GA.label(Visit) = "Visit";
	node ZipCode = graph.newNode();
	GA.label(ZipCode) = "ZipCode";

	edge e0 = graph.newEdge(Address, Property);
	edge e1 = graph.newEdge(Agent, Rent);
	edge e2 = graph.newEdge(Agent, Sale);
	edge e3 = graph.newEdge(Agent, Visit);
	edge e4 = graph.newEdge(City, ZipCode);
	edge e5 = graph.newEdge(Client, ClientWish);
	edge e6 = graph.newEdge(Client, Rent);
	edge e7 = graph.newEdge(Client, Sale);
	edge e8 = graph.newEdge(Client, Visit);
	edge e9 = graph.newEdge(Country, ZipCode);
	edge e10 = graph.newEdge(Feature, ClientWish);
	edge e11 = graph.newEdge(Feature, PropertyFeature);
	edge e12 = graph.newEdge(Person, Agent);
	edge e13 = graph.newEdge(Person, Client);
	edge e14 = graph.newEdge(Person, Owner);
	edge e15 = graph.newEdge(PropertyForRent, Rent);
	edge e16 = graph.newEdge(PropertyForSale, Sale);
	edge e17 = graph.newEdge(Property, Owner);
	edge e18 = graph.newEdge(Property, PropertyFeature);
	edge e19 = graph.newEdge(Property, PropertyForRent);
	edge e20 = graph.newEdge(Property, PropertyForSale);
	edge e21 = graph.newEdge(Property, Visit);
	edge e22 = graph.newEdge(RentalAgreement, Rent);
	edge e23 = graph.newEdge(StateOrProvince, ZipCode);
	edge e24 = graph.newEdge(ZipCode, Address);
}