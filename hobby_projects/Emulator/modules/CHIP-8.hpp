#pragma once
#include <cstring>
#include <random>
#include <iostream>
#include <format>

const uint8_t CHIP8_FONTSET[80] = {
    0xF0, 0x90, 0x90, 0x90, 0xF0, // 0
    0x20, 0x60, 0x20, 0x20, 0x70, // 1
    0xF0, 0x10, 0xF0, 0x80, 0xF0, // 2
    0xF0, 0x10, 0xF0, 0x10, 0xF0, // 3
    0x90, 0x90, 0xF0, 0x10, 0x10, // 4
    0xF0, 0x80, 0xF0, 0x10, 0xF0, // 5
    0xF0, 0x80, 0xF0, 0x90, 0xF0, // 6
    0xF0, 0x10, 0x20, 0x40, 0x40, // 7
    0xF0, 0x90, 0xF0, 0x90, 0xF0, // 8
    0xF0, 0x90, 0xF0, 0x10, 0xF0, // 9
    0xF0, 0x90, 0xF0, 0x90, 0x90, // A
    0xE0, 0x90, 0xE0, 0x90, 0xE0, // B
    0xF0, 0x80, 0x80, 0x80, 0xF0, // C
    0xE0, 0x90, 0x90, 0x90, 0xE0, // D
    0xF0, 0x80, 0xF0, 0x80, 0xF0, // E
    0xF0, 0x80, 0xF0, 0x80, 0x80  // F
};

struct Stack{
    uint16_t stack[16];
    uint16_t *stack_pointer;

    Stack(){
        memset(this->stack, 0, sizeof(this->stack));
        stack_pointer = &this->stack[0]; 
    }

    void push(uint16_t address){
        if (stack_pointer == &stack[16]){
            std::cerr << "[ERROR] - Stack Overflow\n";
            exit(1);
        }
        *stack_pointer = address;
        stack_pointer += 1;
    }

    uint16_t pop(){
        if (stack_pointer == &stack[0]){
            std::cerr << "[ERROR] - Stack Underflow\n";
            exit(1);
        }
        stack_pointer -= 1;
        uint16_t result = *stack_pointer;
        *stack_pointer = 0;
        return result;
    }
};


void LogicalREG(uint8_t* REG, uint16_t remainder){
    uint8_t* Vx = REG + ((remainder & 0x0F00) >> 8);
    uint8_t* Vy = REG + ((remainder & 0x00F0) >> 4);
    uint8_t* VF = REG + 15;

    switch (remainder & (0x000F)){
        
        case 0x0:
            *Vx = *Vy;
            break;

        case 0x1:
            *Vx = *Vx | *Vy;
            break;
        
        case 0x2:
            *Vx = *Vx & *Vy;
            break;
        
        case 0x3:
            *Vx = *Vx xor *Vy;
            break;

        case 0x4:
            if (*Vx + *Vy > 255){
                *VF = 1;
            }
            else{
                *VF = 0;
            }
            *Vx += *Vy;
            break;

        case 0x5:
            *VF = (*Vx >= *Vy) ? 1 : 0;
            *Vx -= *Vy;
            break;

        case 0x6:
            *VF = (*Vx & 0x01) == 0x01 ? 1 : 0;
            *Vx = *Vx >> 1;
            break;

        case 0x7:
            *VF = (*Vy >= *Vx) ? 1 : 0;
            *Vx = *Vy - *Vx;
            break;

        case 0xE:
            *VF = (*Vx >> 7) == 0x01 ? 1 : 0;
            *Vx = *Vx << 1;
            break;

        default:
            break;
    }
}

class CHIP8{

    uint16_t OPCODE;

    uint8_t RAM[4096];

    uint8_t REG[16];

    uint16_t INDEX;
    uint16_t program_counter;

    uint8_t delay_timer = 0;
    uint8_t sound_timer = 0;

    Stack stack;

    public:
    uint8_t running;
    uint8_t Display[64 * 32];
    uint8_t KeyPad[16];
    bool drawFlag;

    CHIP8(){
        memset(this->RAM, 0, sizeof(RAM));
        memset(this->REG, 0, sizeof(REG));
        memset(this->Display, 0, sizeof(Display));
        memset(this->KeyPad, 0, sizeof(KeyPad));

        for (int i = 0; i < 80; i++){
            this->RAM[i] = CHIP8_FONTSET[i];
        }

        this->running = 1;
        this->program_counter = 0x200;
        this->OPCODE = 0;
        this->INDEX = 0;

        srand(time(NULL));
    }

    void LoadProgram(const char *filename){
        FILE *f = fopen(filename, "rb");
        if (f == NULL){
            std::cout << "[ERROR] - Could not open program\n";
            exit(1);
        }

        fseek(f, 0, SEEK_END);
        long fileSize = ftell(f);
        rewind(f);

        if (fileSize > (4096 - 0x200)){
            std::cout << "[ERROR] - File is too large (" << fileSize << " bytes)\n";
            fclose(f);
            exit(1);
        }

        fread(&this->RAM[0x200], 1, fileSize, f);
        std::cout << "[SUCCESS] - LOADED" << fileSize << " bytes from " << filename << std::endl;

        fclose(f);
    }

    void UpdateTimers(){
        if (delay_timer > 0) {
            delay_timer--;
        }
        if (sound_timer > 0) {
            sound_timer--;
        }
    }

    uint16_t Fetch(){
        uint16_t opcode = (RAM[program_counter] << 8) | RAM[program_counter + 1];
        //std::cout << std::format("Opcode: {:#06X}\n", opcode);
        program_counter += 2;

        return opcode;
    }

    void Decode(uint16_t opcode){
        uint16_t remainder = opcode & 0x0FFF;

        switch ((opcode & 0xF000) >> 12){

            case 0x0:{
                if (remainder == 0x00E0){
                    // CLEAR DISPLAY
                    memset(this->Display, 0, sizeof(Display));
                }
                else if (remainder == 0x00EE){
                    // RETURN FROM SYSTEM ROUTINE
                    program_counter = stack.pop();
                }
                else if ((remainder & 0x0FF0) == 0x00C0) { // 00CN: Scroll Down
                    uint8_t n = remainder & 0x000F;
                    for (int i = (64 * 32) - 1; i >= 64 * n; i--) {
                        Display[i] = Display[i - (64 * n)];
                    }
                    memset(Display, 0, 64 * n);
                }
                else{
                    // SPECIFIC ROUTINE ONLY USED ON OLDER PCs
                    // DO NOTHING
                }
                break;
            }

            case 0x1:{  // JUMP TO ADDRESSE NNN
                program_counter = remainder;
                break;
            }

            case 0x2:{ // PUSH TO STACK
                stack.push(program_counter);
                program_counter = remainder;
                break;
            }

            case 0x3:{ // SKIP INSTRUCTION IF REGx = NN
                if (REG[(0x0F00 & remainder) >> 8] == (0x00FF & remainder)){
                    program_counter += 2;
                }
                break;
            }

            case 0x4:{ // SKIP INSTRUCTION IF REGx != NN
                if (REG[(0x0F00 & remainder) >> 8] != (0x00FF & remainder)){
                    program_counter += 2;
                }
                break;
            }

            case 0x5:{ //SKIP INSTRUCTION IF REGx = REGy
                if (REG[(0x0F00 & remainder) >> 8] == REG[(0x00F0 & remainder) >> 4]){
                    program_counter += 2;
                }
                break;
            }

            case 0x6:{ // SET REGx = NN
                REG[(0x0F00 & remainder) >> 8] = (0x00FF & remainder);
                break;
            }

            case 0x7:{ // ADD NN TO REGx
                REG[(0x0F00 & remainder) >> 8] += (0x00FF & remainder);
                break;
            }

            case 0x8:{ // DO LOGICAL OPERATIONS ON REGISTERS
                LogicalREG(&REG[0], remainder);
                break;
            }

            case 0x9:{ // SKIP INSTRUCTION IF REGx != REGy
                if(REG[(remainder & 0x0F00) >> 8] != REG[(remainder & 0x00F0) >> 4]){
                    program_counter += 2;
                }
                break;
            }

            case 0xA:{ // SET INDEX TO NNN
                INDEX = remainder;
                break;
            }

            case 0xB:{ // JUMP TO NN + REG0
                program_counter = remainder + REG[0];
                break;
            }

            case 0xC:{ // SET REGx TO RANDOM + NN
                uint8_t random = rand() % 256;
                REG[(remainder & 0x0F00) >> 8] = random & (remainder & 0x00FF);
                break;
            }

            case 0xD:{ // DISPLAY N-BYTE SPRITE STARTING at (REGx, REGy), SET VF = COLLISION
                uint8_t Vx = REG[(remainder & 0x0F00) >> 8] % 64;
                uint8_t Vy = REG[(remainder & 0x00F0) >> 4] % 32;
                uint8_t height = remainder & 0x000F;
                REG[0xF] = 0;

                for (unsigned int row = 0; row < height; row++){
                    if (Vy + row >= 32) break;

                    uint8_t sprite_byte = RAM[INDEX + row];

                    if (INDEX + row > 4096){
                        std::cerr << "[ERORR] Trying to access out-of-bound memory\n";
                        break;
                    }
                        
                    for (unsigned bit = 0; bit < 8; bit++){
                        if (Vx + bit >= 64) break;

                        if((sprite_byte & (0x80 >> bit)) != 0){
                            uint16_t display_idx = (Vx + bit) + (Vy + row) * 64;

                            if(Display[display_idx] == 1){
                                REG[0xF] = 1;
                            }
                            Display[display_idx] ^= 1;
                            
                        }
                    }
                    
                }
                drawFlag = true;
                break;
            }

            case 0xE:{ // SKIP INSTRUCTION IF KEYx IS PRESSED / NOT PRESSED
                uint8_t Vx = REG[(remainder & 0x0F00) >> 8];

                if ((remainder & 0x00FF) == 0x9E){
                    if (KeyPad[Vx] == 1) program_counter += 2; 
                }
                else if ((remainder & 0x00FF) == 0xA1){
                    if (KeyPad[Vx] == 0) program_counter += 2; 
                }
                break;
            }

            case 0xF:{ 
                uint8_t idx = (remainder & 0x0F00) >> 8;
                uint8_t &Vx = REG[idx];
                uint8_t LastByte = remainder & 0x00FF;
                switch (LastByte){

                    case 0x07:{ // SET REGx = DELAY TIMER
                        Vx = delay_timer;
                        break;
                    }

                    case 0x0A:{ //WAIT FOR KEY PRESS (CONTINOUSLY LOOPS), THEN STORE IN REGx
                        bool keyPressed = false;
                        for (int i = 0; i < 16; i++) {
                            if (KeyPad[i] != 0) {
                                REG[(remainder & 0x0F00) >> 8] = i;
                                keyPressed = true;
                                break;
                            }
                        }

                        if (!keyPressed) {
                            program_counter -= 2; 
                        }
                        break;
                    }
                
                    case 0x15:{ // SET DELAY TIMER = REGx
                        delay_timer = Vx;
                        break;
                    }
                
                    case 0x18:{ // SET SOUND TIMER = REGx
                        sound_timer = Vx;
                        break;
                    }
                
                    case 0x1E:{ // ADD REGx TO INDEX 
                        if (INDEX + Vx > 0x0FFF){
                            REG[0xF] = 1;
                        }
                        else{
                            REG[0xF] = 0;
                        }
                        INDEX += Vx;
                        break;
                    }

                    case 0x29:{ // SET INDEX = LOCATION OF SPRITE FOR DIGIT REGx
                        INDEX = Vx * 5;
                        break;
                    }

                    case 0x33:{ // STORE BCD REPRESENTIATION IN REGx
                        RAM[INDEX]     = Vx / 100;          
                        RAM[INDEX + 1] = (Vx / 10) % 10;    
                        RAM[INDEX + 2] = Vx % 10;
                        break;
                    }

                    case 0x55:{ // STORE REG0 - REGx IN MEMORY STARTING AT INDEX 
                        uint8_t *reg_ptr = &REG[0]; 
                        for (int i = 0; i <= idx; i++){
                            RAM[INDEX + i] = *(reg_ptr + i);
                        }
                        break;
                    }

                    case 0x65:{ // READ MEMORY STARTING AT INDEX, INSERT INTO REG0 - REGx
                        uint8_t *reg_ptr = &REG[0]; 
                        for (int i = 0; i <= idx; i++){
                            *(reg_ptr + i) = RAM[INDEX + i];
                        }
                        break;
                    }
                
                    default:{
                        break;
                    }
                }
                break;
            }    
        }
    }

    ~CHIP8(){
    }
};


