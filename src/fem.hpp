#pragma once
#include <vector>
#include <Eigen/Dense>

struct degreesoffreedom{
    bool xfixed = false;
    bool yfixed = false;
    bool thetafixed = false;
};

struct node{

    int x;
    int y;
    degreesoffreedom dof = degreesoffreedom{};
};

struct connection{
    node inode;
    node jnode;
};

class Fem{

    std::vector<node> nodes;
    std::vector<connection> connections;

    public:
        Fem();
        void addnode(node n);
        void addnode(int x, int y);
        void connect(node i, node j);
        void connect(int x1, int y1, int x2, int y2);
        std::vector<node> getNodes();
        std::vector<connection> getConnections();
        Eigen::MatrixXf getStiffness_1DOF(connection c);
        Eigen::MatrixXf getGlobalStiffnessMatrix();
        Eigen::MatrixXf trimFixedDofFromMatrix(Eigen::MatrixXd global);

    protected:
        node getOrCreateNode(node n);
        connection getOrCreateConnection(connection c);
        int getGlobalDegreesOfFreedomCount();
        std::vector<connection> getConnectionsForNode(node n);
        int getNodeIndex(node n);
        
};