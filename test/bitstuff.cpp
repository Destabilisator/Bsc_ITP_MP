#include <iostream>
#include <fstream>
#include <ctime>
#include <random>
#include <cmath>
#include <thread>
#include <mutex>
#include <string>
#include <bitset>

void printBits(int a) {
    std::bitset<8> x(a);
    std::cout << x << '\n';
}

void bitsStuff() {
    int a = 5, b = 9;
    std::cout << "a = " << a << " = ";
    printBits(a);
    std::cout << "b = " << b << " = ";
    printBits(b);

    std::cout << "\n";

    std::cout << "and & :\n";
    std::cout << "a & b = ";
    printBits(a & b);

    std::cout << "\n";


    std::cout << "or | :\n";
    std::cout << "a | b = ";
    printBits(a | b);

    std::cout << "\n";

    std::cout << "xor ^ :\n";
    std::cout << "a ^ b = ";
    printBits(a ^ b);

    std::cout << "\n";

    std::cout << "inverse ~ :\n";
    std::cout << "~(a) = ";
    printBits(~a);
    std::cout << "~(b) = ";
    printBits(~b);

    std::cout << "\n";


    std::cout << "bitshift left << :\n";
    std::cout << "a << 1" << " = ";
    printBits(a << 1);
    std::cout << "b << 1" << " = ";
    printBits(b << 1);

    std::cout << "\n";

    std::cout << "bitshift right >> :\n";
    std::cout << "a >> 1 " << "= ";
    printBits(a >> 1 );
    std::cout << "b >> 1 " << "= ";
    printBits(b >> 1 );

    std::cout << "\n";

    int c = 0;
    std::cout << "c = " << c << " = ";
    printBits(c);

    std::cout << "\n";

    std::cout << "set n th bit (n = 4):\n";
    c |= 1UL << 4;
    printBits(c);

    std::cout << "\n";

    std::cout << "clear n th bit (n = 4):\n";
    c &= ~(1UL << 4);
    printBits(c);

    std::cout << "\n";

    std::cout << "flip n th bit (n = 4):\n";
    c ^= 1UL << 4;
    printBits(c);
    c ^= 1UL << 4;
    printBits(c);
    c ^= 1UL << 4;
    printBits(c);

    std::cout << "\n";

    std::cout << "check n th bit (n = 4):\n";
    bool bit = (c >> 4) & 1UL;
    std::cout << bit << "\n";
}

int main() {

    bitsStuff();

    return 0;
}
