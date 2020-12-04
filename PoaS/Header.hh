#ifndef PoaSHeader
#define PoaSHeader

using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <list>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <unistd.h>

// Include the same RNG as used in Brem-cash
extern "C"{
#include "../Brem-cash/mersenne.h"
}

#define TRUE 1
#define FALSE 0
#endif
