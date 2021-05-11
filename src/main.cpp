#include <raylib-cpp.hpp>
#include <iostream>
#include <Eigen/Dense>
#include "fem.hpp"
#include <vector>

using Eigen::MatrixXd;


static std::string toString(const Eigen::MatrixXf& mat){
    Eigen::IOFormat cleanFmt(4,0,", ","\n","|", "|");
    std::stringstream ss;
    ss << mat.format(cleanFmt);
    return ss.str();
}

int main() {
    
    // Initialization
    int screenWidth = 1600;
    int screenHeight = 900;

    raylib::Color textColor(DARKGRAY);
    raylib::Window w(screenWidth, screenHeight, "FEM");

    Fem fem = Fem();
    fem.addnode(node{100, 100, degreesoffreedom{true, true, true}}); //fixed endpoint
    fem.addnode(node{200,100});
    fem.addnode(node{400,100});
    fem.connect(node{100, 100}, node{200,100});
    fem.connect(node{200, 100}, node{400,100});

    SetTargetFPS(60);


    // Main game loop
    while (!w.ShouldClose()) // Detect window close button or ESC key
    {
        
        // one spring
        // k = F/d
        // 

        BeginDrawing();
        ClearBackground(RAYWHITE);
        std::vector<node> nodes = fem.getNodes();
        for(std::vector<node>::iterator it = nodes.begin(); it != nodes.end(); ++it){
            node n = *it;
            DrawCircle(n.x, n.y, 2, BLUE);
            
        }
         
        std::vector<connection> connections = fem.getConnections();
        for(std::vector<connection>::iterator it = connections.begin(); it != connections.end(); ++it){
            connection c = *it;
            DrawLine(c.inode.x, c.inode.y, c.jnode.x, c.jnode.y, GREEN);
            textColor.DrawText(
                toString(fem.getStiffness_1DOF(c)), 
                (c.inode.x + c.jnode.x)/2 - 20, 
                ((c.inode.y + c.jnode.y)/2) + 20, 
                10
            );
        }
        
        textColor.DrawText(
                toString(fem.getGlobalStiffnessMatrix()), 
                200, 
                400, 
                10
            );

        EndDrawing();
    }

    return 0;
}