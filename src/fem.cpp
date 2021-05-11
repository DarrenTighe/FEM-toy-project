#include "fem.hpp"

Fem::Fem(){

}

void Fem::addnode(node n){
    this->getOrCreateNode(n);
}

void Fem::addnode(int x, int y){
    this->getOrCreateNode(node{x=x, y=y});
}

node Fem::getOrCreateNode(node n){
    int x = n.x;
    int y = n.y;
    node res = {};
    bool found = false;
    for(std::vector<node>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        node n = *it;
        if (n.x == x && n.y == y){
            found = true;
            res = n;
            break;
        }
    };
    if (!found){
        res = n;
        nodes.push_back(n);
    }
    return res;
}

int Fem::getNodeIndex(node n){
    int index = -1;
    int cnt = 0;
    for(std::vector<node>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        node nd = *it;
        if (n.x == nd.x && n.y == nd.y){
            index = cnt;
            break;
        }
        cnt++;
    }
    return index;
}

connection Fem::getOrCreateConnection(connection c){
    connection res = {};
    bool found = false;
    for(std::vector<connection>::iterator it = connections.begin(); it != connections.end(); ++it){
        connection con = *it;
        if (con.inode.x == c.inode.x
            && con.inode.y == c.inode.y 
            && con.jnode.x == c.jnode.x 
            && con.jnode.y == c.jnode.y
        ){
            found = true;
            res = con;
            break;
        }
    }
    if(!found){
        res = c;
        connections.push_back(c);
    }
    return res;

    return c;
}

void Fem::connect(node i, node j){
    connection c = {};
    c.inode = i;
    c.jnode = j;
    getOrCreateConnection(c);
}

std::vector<node> Fem::getNodes(){
    return nodes;
}

std::vector<connection> Fem::getConnections(){
    return connections;
}

Eigen::MatrixXf Fem::getStiffness_1DOF(connection c){
    // spring material - super rough approach
    // k  = AE/L
    float L = sqrtf(powf((c.inode.x - c.jnode.x), 2)+ powf((c.inode.y - c.jnode.y), 2));
    float AE = 70000.0;
    Eigen::MatrixXf k(2,2);
    k(0,0) = AE/L;
    k(0,1) = -AE/L;
    k(1,0) = -AE/L;
    k(1,1) = AE/L;

    return k;
}

int Fem::getGlobalDegreesOfFreedomCount(){
    int dof = 0;
    for(std::vector<node>::iterator it = nodes.begin(); it != nodes.end(); ++it){
        dof++;
    }
    return dof;
}

std::vector<connection> Fem::getConnectionsForNode(node n){
    std::vector<connection> res;
    for (std::vector<connection>::iterator it = connections.begin(); it !=connections.end(); ++it){
        connection c = *it;
        if ((c.inode.x == n.x && c.inode.y == n.y) || (c.jnode.x == n.x && c.jnode.y == n.y)){
            res.push_back(c);
        }
    }
    return res;
}

Eigen::MatrixXf Fem::getGlobalStiffnessMatrix(){
    int totaldof = this->getGlobalDegreesOfFreedomCount();
    Eigen::MatrixXf k = Eigen::MatrixXf::Zero(totaldof, totaldof);
   
    std::vector<Eigen::MatrixXf> transformedStiffnessMatrices;
    for (std::vector<connection>::iterator it = connections.begin(); it != connections.end(); ++it){
        
        // translate local to global, e.g
        // | a, b | | u1 |   | a, 0, b | | u1 |
        // | c, d | | u3 | = | 0, 0, 0 | | u2 |
        //                   | c, 0, d | | u3 |
        connection c = *it;
        int ind_i = getNodeIndex(c.inode);
        int ind_j = getNodeIndex(c.jnode);
        Eigen::MatrixXf kc = Eigen::MatrixXf::Zero(totaldof, totaldof);
        Eigen::MatrixXf local = getStiffness_1DOF(c);

        kc(ind_i, ind_i) = local(0,0);
        kc(ind_i, ind_j) = local(0,1);
        kc(ind_j, ind_i) = local(1,0);
        kc(ind_j, ind_j) = local(0,1);

        // add this connection's stiffness to the global stiffness
        k = k + kc;
    }

    return k;
}



