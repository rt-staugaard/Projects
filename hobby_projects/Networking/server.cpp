#include "modules/physics.hpp"
#include "modules/network/server.hpp"

int main(){

    std::unique_ptr<Grid> current_grid = std::make_unique<Grid>();
    std::unique_ptr<Grid> previous_grid = std::make_unique<Grid>();
    std::unique_ptr<Grid> buffer = std::make_unique<Grid>();

    auto forceX = std::make_unique<float[]>(gridLength * gridHeight);
    auto forceY = std::make_unique<float[]>(gridLength * gridHeight);

    Interaction packet;
    
    std::unique_ptr<Server> myServer = std::make_unique<Server>(8080);

    std::thread listenerThread(&Server::startListening, myServer.get());
    listenerThread.detach();

    std::thread terminalThread([&](){
        std::string input;
        while(myServer->getStatus()){
            if(std::getline(std::cin, input)){
                if (input == "exit" || input == "quit") {
                    myServer->shutdown(); 
                }
            }
        }
    });
    terminalThread.detach();

    while (myServer->getStatus()){
        auto start = std::chrono::high_resolution_clock::now();

        if (myServer->getClientNumber() > 0){
            packet = myServer->recieveInput();
        }

        if (!packet.isDragging && packet.targetIdx != -1) {
            current_grid->Xs[packet.targetIdx] = packet.targetX;
            current_grid->Ys[packet.targetIdx] = packet.targetY;

            previous_grid->Xs[packet.targetIdx] = packet.targetX;
            previous_grid->Ys[packet.targetIdx] = packet.targetY;
            packet.targetIdx = -1;
        }

        physicsStep(current_grid.get(), previous_grid.get(), buffer.get(), forceX.get(), forceY.get(), 4);
        
        std::swap(previous_grid, current_grid); 
        std::swap(current_grid, buffer);

        myServer->broadcast(current_grid.get()->Xs, current_grid.get()->Ys, gridHeight * gridLength);

        std::this_thread::sleep_for(std::chrono::milliseconds(4));
    }
    myServer->shutdown();
}

