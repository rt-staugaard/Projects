#pragma once
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <cstring>
#include <iostream>
#include <vector>
#include <mutex>
#include <atomic>
#include <fcntl.h>
#include "../graphics/InteractionPacket.hpp"

class Server{
private:
    int server_fd;
    std::vector<int> clients; 
    const int maxClients = 5;
    std::atomic<bool> running = true;
    std::mutex clientsMutex;

public:
    Server(int port){
        server_fd = socket(AF_INET, SOCK_STREAM, 0);
        if (server_fd < 0){
            perror("socket");
            std::cout << "Could not create server socket\n";
            return;
        }

        sockaddr_in addr{};
        addr.sin_family = AF_INET;
        addr.sin_addr.s_addr = INADDR_ANY;
        addr.sin_port = htons(port);

        if (bind(server_fd, (sockaddr*)&addr, sizeof(addr)) < 0){
            perror("bind");
            std::cout << "Could not bind to addresse\n";
            close(server_fd);
            return;
        }
    }

    void startListening() {
        if (::listen(server_fd, 5) < 0){
            perror("Too many clients trying to access server");
            return;
        }
        std::cout << "Server Online:\n";

        while (running) {
            int client_fd = accept(server_fd, nullptr, nullptr); 
        
            if (client_fd >= 0 && clients.size() < 5) {
                addClient(client_fd); 
                std::cout << "Client added.\n";
            }
        }
        running = false;
    }

    void addClient(int fd) {
        std::lock_guard<std::mutex> lock(clientsMutex);
        int flags = fcntl(fd, F_GETFL, 0);          
        fcntl(fd, F_SETFL, flags | O_NONBLOCK);
        clients.push_back(fd);
    }

    void broadcast(float* Xs, float* Ys, unsigned int count) {
        std::lock_guard<std::mutex> lock(clientsMutex); 
        for (auto it = clients.begin(); it != clients.end();) {
            
            ssize_t r1 = send(*it, Xs, count * sizeof(float), MSG_NOSIGNAL);
            ssize_t r2 = send(*it, Ys, count * sizeof(float), MSG_NOSIGNAL);

            if (r1 < 0 || r2 < 0){
                std::cout << "Client " << *it << " disconnected. Removing...\n";
                close(*it);
                it = clients.erase(it);
            }
            else{
                ++it;
            }
        }
    }

    bool isCompleteMessage(ssize_t bytesRead, unsigned int expectedSize) {
        return bytesRead == (ssize_t)expectedSize;
    }

    Interaction recieveInput(){
        Interaction packet; 
        size_t packet_size = sizeof(Interaction);

        std::lock_guard<std::mutex> lock(clientsMutex);
        for (auto it = clients.begin(); it != clients.end();){
            ssize_t r = recv(*it, &packet, packet_size, 0);
            
            if (r > 0) {
                if(!isCompleteMessage(r, packet_size)){
                    continue;
                }
                return packet;
            }
            else if (r == 0){
                std::cout << "Client " << *it << " closed connection.\n";
                close(*it);
                it = clients.erase(it);
                continue;
            }
            else {

                if (errno == EAGAIN || errno == EWOULDBLOCK){
                    ++it;
                }
                else{
                    std::cout << "Client " << *it << " error. Removing...\n";
                    close(*it);
                    it = clients.erase(it);
                }
            }
        }
        return packet;
    }

    int getClientNumber(){
        return clients.size();
    }

    bool getStatus(){
        if(running == true){
            return true;
        }
        return false;
    }

    void shutdown(){
        running = false;
    }

    ~Server() {
        for(int i = 0; i < clients.size(); i++) {
            close(clients[i]);
        }
        close(server_fd);
    }
};
