/******************************************************************
* Minna Barry
* CS331 - Programming Assignment 1
*******************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <cstring>
#include <fstream>
#include <unordered_map>
#include <stack>
#include <queue>
#include <cstdlib>

typedef struct Node * Nodeptr;
typedef std::pair<int,Nodeptr> priotiyPair;

class Problem{
    public:
    std::string initState;
    std::string goalState;

    bool goalTest(std::string state){
        if(state == goalState)
            return true;
        return false;
    }

    int heuristic(std::string state){
        int current [6];    //represents current state
        int goal [6];       //represents goal state
        int result=0;
        int n=0;            //used for parsing strings
        //heuristic of the goal state is 0
        if(state != goalState)
            //calculate the least number of moves needed to move 1 cannibal or missionary to the goal state
            for(int i=0; i<3; i++){
                //parse state and goal strings
                current[i] = atoi(state.substr(n,state.find_first_of(",\n",n)).c_str());
                goal[i] = atoi(goalState.substr(n,goalState.find_first_of(",\n",n)).c_str());
                n=state.find_first_of(",\n",n)+1;

                if(i==2){
                    //if the boat was at the left bank add one more move to the result
                    //which corresponds to moving the boat back to the right
                    current[i]=!current[i];
                }
                result += goal[i]-current[i];
            }
        return result;
    }

    std::queue<std::string> * Successors(std::string state){
        //0 leftM, 1 leftC, 2 leftB, 3 rightM, 4 rightC, 5 rightB
        std::queue<std::string> * results  = new std::queue<std::string>;
        int current [6];
        int result [6];
        int n=0;
        bool manymv=0,who=0,mvtoleft=0;
        std::stringstream ss;
        //parse state string
        for(int i=0; i<6; i++){
            current[i] = atoi(state.substr(n,state.find_first_of(",\n",n)).c_str());
            n=state.find_first_of(",\n",n)+1;
        }
        mvtoleft = current[5];
        //each loop to apply one action
        for(int i=0; i<5; i++){
            //copy current state to result
            for(int j=0; j<6; j++){
                result[j] = current[j];
            }
            //apply action
            if(mvtoleft){
                if(i==3){
                    //move 1 m and 1 c to left
                    result[0]+=1;
                    result[0+3]-=1;
                    result[1]+=1;
                    result[1+3]-=1;
                }else{
                    result[who]+=manymv+1;
                    result[who+3]-=manymv+1;
                    if(manymv){
                        who = !who;
                    }
                    manymv=!manymv;
                }
                result[2]=1;
                result[5]=0;
            }else{
                if(i==3){
                    //move 1 m and 1 c to right
                    result[0]-=1;
                    result[0+3]+=1;
                    result[1]-=1;
                    result[1+3]+=1;
                }
                else{
                    result[who+3]+=manymv+1;
                    result[who]-=manymv+1;
                    if(manymv){
                        who = !who;
                    }
                    manymv=!manymv;
                }
                result[2]=0;
                result[5]=1;
            }
            //check if rules apply
            bool addresult=1;
            if((result[0] >= result[1] || result[0]==0) && (result[3] >= result[4] || result[3]==0))
                for(int j=0; j<6; j++){
                    if(result[j] < 0)
                        addresult=0;
                }
            else
                addresult=0;
            //if action passes add the resulting state to the result queue
            if(addresult){
                ss<< result[0]<<","<<result[1]<<","<<result[2]<<"\n"<<result[3]<<","<<result[4]<<","<<result[5]<<"\n";
                results->push(ss.str());
                ss.str("");
            }
            //update variables
        }
        return results;
    }
};

struct Node{
    std::string state;
    std::string action;
    struct Node * parent;
    float cost;
    int depth;
};

Nodeptr Node_init(std::string s,std::string a, Nodeptr p, float f, int d){
    Nodeptr n = new struct Node;
    n->state = s;
    n->action = a;
    n->parent = p;
    n->cost = f;
    n->depth = d;
    return n;
}

void Node_delete(Nodeptr n){
    Nodeptr p = n->parent;
    Nodeptr temp = n->parent;
    //delete(n);
    while(p!=0){
        temp = p;
        p = p->parent;
        delete(temp);
    }
    delete(n);
}

struct Solution{
    Nodeptr goal;
    int expanded;
};

void printSolution(struct Solution s){
    Nodeptr n = s.goal;
    if(n!=0){
        std::stack<std::string>l;
        while(n!=0){
            l.push(n->state);
            n=n->parent;
        }
        while(!l.empty()){
            std::cout << l.top() <<"\n";
            l.pop();
        }
        std::cout<< "Expanded: "<< s.expanded<<"\nPath length: "<<s.goal->depth<<"\n";
    }else{
        std::cout<<"No solution was found..\n";
    }
}

void writeSolution(char * fileName,struct Solution s){
    std::ofstream out;
    out.open(fileName);
    Nodeptr n = s.goal;
    if(n!=0){
        std::stack<std::string>l;
        while(n!=0){
            l.push(n->state);
            n=n->parent;
        }
        while(!l.empty())
        {
            out << l.top() <<"\n";
            l.pop();
        }
        out<< "Expanded: "<< s.expanded<<"\nPath length: "<<s.goal->depth<<"\n";
    }else{
        out <<"No solution was found..\n";
    }
}

std::queue<Nodeptr> * Expand(Nodeptr n, Problem p){

    std::queue<Nodeptr> * successors = new std::queue<Nodeptr>;
    std::queue<std::string> * results = p.Successors(n->state);

    while(!results->empty()){

        Nodeptr node = Node_init(results->front(),"",n,n->cost+ 1, n->depth+1);
        results->pop();
        successors->push(node);

    }
    return successors;
}

struct Solution GraphDepthSearch(Problem p){
    std::unordered_map<std::string, int> closed;
    std::stack<Nodeptr> fringe;
    struct Solution s;
    s.expanded=0;
    //add initial state node to fringe
    fringe.push(Node_init(p.initState,"",0,0,0));

    while(1){
        if(fringe.empty()){
            s.goal=0;
            return s;
        }
        Nodeptr n = fringe.top();
        fringe.pop();
        if(p.goalTest(n->state)){
            s.goal = n;
            return s;
        }
        if(closed.find(n->state)==closed.end()){
            //add it to closed
            closed[n->state]=1;
            //expand
            s.expanded++;
            std::queue<Nodeptr> * children = Expand(n,p);
            while(!children->empty()){
                fringe.push(children->front());
                children->pop();
            }
            delete(children);
        }
    }
}

struct Solution GraphBreadthSearch(Problem p){
    std::unordered_map<std::string, int> closed;
    std::queue<Nodeptr> fringe;
    struct Solution s;
    s.expanded=0;
    //add initial state node to fringe
    fringe.push(Node_init(p.initState,"",0,0,0));
    while(1){
        if(fringe.empty()){
            s.goal=0;
            return s;
        }
        Nodeptr n = fringe.front();
        fringe.pop();
        if(p.goalTest(n->state)){
            s.goal=n;
            return s;
        }
        if(closed.find(n->state)==closed.end()){
            closed[n->state]=1;
            //expand
            s.expanded++;
            std::queue<Nodeptr> * children = Expand(n,p);
            while(!children->empty()){
                fringe.push(children->front());
                children->pop();
            }
        }
    }

}

struct Solution GraphIDDFSearch(Problem p){
    std::unordered_map<std::string, int> closed;
    std::stack<Nodeptr> fringe;
    struct Solution s;
    s.expanded=0;
    int limit=0;
    //add initial state node to fringe
    fringe.push(Node_init(p.initState,"",0,0,0));
    Nodeptr n = fringe.top();
    for(;limit <2000;limit++){
        while(1){
            if(fringe.empty()){
                break;
            }
            Nodeptr n = fringe.top();
            fringe.pop();
            //test if current node is the goal
            if(p.goalTest(n->state)){
                s.goal=n;
                return s;
            }
            //if current state is not in closed
            if(closed.find(n->state)==closed.end()){
                //add it to closed
                closed[n->state]=1;
                //expand
                if(n->depth < limit){
                    s.expanded++;
                    std::queue<Nodeptr> * children = Expand(n,p);
                    while(!children->empty()){
                        fringe.push(children->front());
                        children->pop();
                    }
                }
            }
        }
        fringe.push(Node_init(p.initState,"",0,0,0));
        n = fringe.top();
        closed.erase(closed.begin(),closed.end());
        //s.expanded=0;
    }
    s.goal=0;
    return s;
}

struct Solution GraphAStarSearch(Problem p){
    std::unordered_map<std::string, int> closed;
    std::priority_queue<priotiyPair,std::vector<priotiyPair>,std::greater<priotiyPair>> fringe;
    struct Solution s;
    s.expanded=0;
    //add initial state node to fringe
    fringe.push(priotiyPair(p.heuristic(p.initState),Node_init(p.initState,"",0,0,0)));
    while(1){
        if(fringe.empty()){
            s.goal=0;
            return s;
        }
        Nodeptr n = fringe.top().second;
        fringe.pop();
        //test if current node is the goal
        if(p.goalTest(n->state)){
            s.goal=n;
            return s;
        }
        //if current state is not in closed
        if(closed.find(n->state)==closed.end()){
            //add it to closed
            closed[n->state]=1;
            //expand
            s.expanded++;
            std::queue<Nodeptr> * children = Expand(n,p);
            while(!children->empty()){
                fringe.push(priotiyPair(p.heuristic(children->front()->state)+children->front()->cost,children->front()));
                children->pop();
            }
        }
    }

}

std::string readFile(char * fileName){
    std::string s;
    std::stringstream ss;
    std::ifstream in;
    in.open(fileName);
    getline(in,s);
    ss << s.substr(0,s.find_first_of("\r\n"))<<"\n";
    getline(in,s);
    ss << s.substr(0,s.find_first_of("\r\n"))<<"\n";
    return ss.str();

}

int main(int argc, char* argv[]){

    if(argc < 4){
        std::cout<<"Usage: miscan <initial state file> <goal state file> <mode> <output file>\n";
        return -1;
    }
    char * startFile,* goalFile,* mode,* outfile;
    startFile = argv[1];
    goalFile = argv[2];
    mode = argv[3];
    outfile = argv[4];

    Problem p;

    p.initState = readFile(startFile);
    p.goalState = readFile(goalFile);
    struct Solution sol;

    if(strcmp(mode,"bfs")==0)
        sol = GraphBreadthSearch(p);
    if(strcmp(mode,"dfs")==0)
        sol = GraphDepthSearch(p);
    if(strcmp(mode,"iddfs")==0)
        sol = GraphIDDFSearch(p);
    if(strcmp(mode,"astar")==0)
        sol = GraphAStarSearch(p);

    printSolution(sol);
    writeSolution(outfile,sol);
    if(sol.goal!=0)
        Node_delete(sol.goal);

}

