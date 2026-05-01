#pragma once
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <iostream>
#include "../graphics/InteractionPacket.hpp"
#include <atomic>

class Client{
private:
    int client_fd;
    bool isConnected = false;

public:

    Client(){
        client_fd = socket(AF_INET, SOCK_STREAM, 0);
        if (client_fd < 0){
            perror("socket");
            std::cout << "Could not create socket\n";
            return;
        }
    }

    void ConnectToServer(const char* address, int port){
        sockaddr_in server{};
        server.sin_family = AF_INET;
        server.sin_port = htons(port);
        
        if(inet_pton(AF_INET, address, &server.sin_addr) <= 0){
            std::cout << "Invalid address\n";
            return;        
        }

        if (connect(client_fd, (sockaddr*)&server, sizeof(server)) < 0){
            perror("connect");
            return;
        }
        isConnected = true;
    }

    bool receiveFullGrid(float* Xs, float* Ys, unsigned int count){
        size_t element_size = count * sizeof(float);

        if (!recvAll((char*)Xs, element_size)) return false;
        if (!recvAll((char*)Ys, element_size)) return false;

        return true;
    }

    bool recvAll(char* buffer, size_t total_size){
        ssize_t bytes_received = 0;
        
        while(total_size > bytes_received){
            ssize_t r = recv(client_fd, buffer + bytes_received, total_size - bytes_received, 0);

            if (r <= 0){
                return false;
            }
            bytes_received += r;
        }
        return true;
    }

    void sendPacket(Interaction inter){
        ssize_t r = send(client_fd, &inter, sizeof(inter), MSG_NOSIGNAL);
    }

    bool getConnectivity(){
        return isConnected;
    }

    ~Client(){
        close(client_fd);
    }
}; 

void NetworkWorker(Client &client, float* &Xs, float* &Ys, float* &XsBuffer, float* &YsBuffer, 
                    std::mutex &dataMtx, unsigned int count, bool& newFrameFlag){

    while(client.getConnectivity()){
        if (client.receiveFullGrid(XsBuffer, YsBuffer, count)){

            {
                std::lock_guard<std::mutex> lock(dataMtx);
                
                std::swap(Xs, XsBuffer);
                std::swap(Ys, YsBuffer);
                
                newFrameFlag = true; 
            }

        }
    }

}