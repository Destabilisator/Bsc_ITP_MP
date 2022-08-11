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

void printBits2(int a) {
    for (int n = 8-1; n >= 0; n--) {
        //printBits((a >> n));
        std::cout << ( (a >> n) & 1 );
    } std::cout << "\n";
}

int invert(int s, int N) {
    //printBits(s);
    int a = pow(2,N)-1;
    //printBits(a);
    return s ^ a;
}

int invertBits(int s, int N) {
    int t = 0;
    for (int n = 0; n < N; n++) {
        t |= ( (s >> n) & 1 ^ 1) << n;
    }
    return t;
}

// cyclicly translate bits in s by n
int translate(int s, int n, int N) {
    for (int _ = 0; _ < n; _++) {
        int bit = s & 1;
        s = (s>>1) | (bit << (N-1));
    } return s;
}

int translate(int s, int N) {
    int bit = s & 1;
    s = (s>>1) | (bit << (N-1));
    return s;
}

// sum up all bits in s
int bitSum(int s, int N) {
    int sum = 0;
    for (int i = 0; i < N; i++) {
        sum += ((s >> i) & 1 );
    } return sum;
}

void bitsStuff() {
    int a = 5, b = 9;
    std::cout << "a = " << a << " = ";
    printBits(a);
    printBits2(a);
    std::cout << "b = " << b << " = ";
    printBits(b);
    printBits2(b);

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

    int N = 25;
    int size = std::pow(2, N);

    auto start1 = std::chrono::steady_clock::now();

    for (int i = 0; i < size; i++) {
        invert(i, N);
    }

    auto end1 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds1 = end1-start1;
    std::cout << "s^N: calculations done; this took: " << elapsed_seconds1.count() << "\n\n";

    auto start2 = std::chrono::steady_clock::now();

    for (int i = 0; i < size; i++) {
        invertBits(i, N);
    }

    auto end2 = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds2 = end2-start2;
    std::cout << "loop: calculations done; this took: " << elapsed_seconds2.count() << "\n\n";

    //bitsStuff();

//    printBits(invert(1, 8));
//    printBits2(invert(1, 8));

//    printBits(34);
//    printBits2(34);

    return 0;
}
